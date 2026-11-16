"""Official Evo 2 adapter and causal, strand-aware context-offset experiments."""

from dataclasses import dataclass
import gc
import importlib.metadata
import logging
import math
from pathlib import Path
from typing import Protocol
import csv

from .inputs import Reference, Variant

LOG = logging.getLogger(__name__)
MODEL_SIZES = (("evo2_40b", 40_000_000_000), ("evo2_7b", 7_000_000_000),
               ("evo2_1b_base", 1_000_000_000))
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


@dataclass(frozen=True)
class Window:
    start: int
    end: int
    target: int
    strand: str
    offset: int


def windows(position: int, length: int, context: int,
            offsets: list[int]) -> list[Window]:
    if context < 2 or len(set(offsets)) != len(offsets):
        raise ValueError("Context must be >= 2 and offsets must be unique")
    result = []
    for offset in offsets:
        if not 1 <= offset < context:
            raise ValueError("Offsets must satisfy 1 <= offset < context length")
        for strand in ("+", "-"):
            target = offset if strand == "+" else context - 1 - offset
            start = position - 1 - target
            if start >= 0 and start + context <= length:
                result.append(Window(start, start + context, offset, strand, offset))
    return result


class Predictor(Protocol):
    def log_probs(self, sequence: str, target: int) -> dict[str, float]: ...


class EvoPredictor:
    def __init__(self, model: str, allow_rocm: bool, context: int):
        import torch

        if not 2 <= context <= 1_000_000:
            raise ValueError("Context must be between 2 and 1,000,000 tokens")
        self.torch = torch
        if not torch.cuda.is_available():
            raise RuntimeError("Evo 2 requires an available GPU backend; no CPU fallback")
        self.rocm = torch.version.hip is not None
        if self.rocm and not allow_rocm:
            raise RuntimeError("ROCm is experimental for official Evo 2; use --allow-rocm "
                               "only after validating your Vortex/kernel installation")
        from evo2 import Evo2

        candidates = MODEL_SIZES if model == "auto" else tuple(
            item for item in MODEL_SIZES if item[0] == model)
        if not candidates:
            raise ValueError(f"Unknown Evo 2 model: {model}")
        self.attempts: list[dict[str, str]] = []
        self.selected = ""
        self.model = None
        free = sum(torch.cuda.mem_get_info(i)[0] for i in range(torch.cuda.device_count()))
        for name, parameters in candidates:
            # Conservative BF16 weights + reserve; successful inference is the actual fit test.
            if self.rocm and name != "evo2_7b":
                reason = "Upstream requires NVIDIA FP8/Transformer Engine for this checkpoint"
            elif name == "evo2_1b_base" and context > 8192:
                reason = "evo2_1b_base supports at most 8192 context tokens"
            elif name != "evo2_7b" and any(
                    torch.cuda.get_device_capability(i)[0] < 9
                    for i in range(torch.cuda.device_count())):
                reason = "Upstream FP8 checkpoint requires Hopper-class or newer hardware"
            elif parameters * 2 + 3 * 1024**3 > free:
                reason = "Insufficient aggregate free memory for BF16 weights plus 3 GiB reserve"
            else:
                try:
                    self.model = Evo2(name)
                    self.model.model.eval()
                    self.selected = name
                    probe = ("ACGT" * ((context + 3) // 4))[:context]
                    self.log_probs(probe, context - 1)
                    LOG.info("Loaded %s and passed a %d-token context probe", name, context)
                    return
                except torch.OutOfMemoryError as error:
                    self.model = None
                    gc.collect()
                    torch.cuda.empty_cache()
                    reason = f"Model load ran out of memory: {error}"
            LOG.warning("Skipping %s: %s", name, reason)
            self.attempts.append({"model": name, "reason": reason})
        raise RuntimeError(f"No eligible Evo 2 checkpoint fits: {self.attempts}")

    def log_probs(self, sequence: str, target: int) -> dict[str, float]:
        if self.model is None:
            raise RuntimeError("Evo 2 model not loaded")
        tokens = self.model.tokenizer.tokenize(sequence)
        if len(tokens) != len(sequence) or target < 1:
            raise ValueError("Expected single-base tokenization and a preceding context token")
        inputs = self.torch.tensor([tokens], dtype=self.torch.long, device="cuda:0")
        with self.torch.inference_mode():
            outputs, _ = self.model(inputs)
            logits = outputs[0]
            # Autoregressive output at t-1 predicts the token at t, not at t+1.
            probabilities = logits[0, target - 1].float().log_softmax(dim=-1)
            result = {base: float(probabilities[self.model.tokenizer.tokenize(base)[0]].item())
                      for base in "ACGT"}
        if any(not math.isfinite(value) for value in result.values()):
            raise RuntimeError("Evo 2 produced nonfinite log probabilities")
        return result

    def metadata(self) -> dict:
        return {"model": self.selected, "evo2_version": importlib.metadata.version("evo2"),
                "torch_version": self.torch.__version__, "rocm": self.rocm,
                "devices": [self.torch.cuda.get_device_name(i)
                            for i in range(self.torch.cuda.device_count())],
                "skipped_models": self.attempts}


def score_sites(sites: list[Variant], reference: Reference, predictor: Predictor,
                output: Path, context: int, offsets: list[int]) -> int:
    written = 0
    with output.open("w") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["chrom", "position", "strand", "window_start0", "window_end0",
                         "offset", "context_length", "ref", "alt", "logp_ref", "logp_alt",
                         "delta_logp", "p_A", "p_C", "p_G", "p_T", "status"])
        for site in sites:
            reference.check(site)
            planned = windows(site.position, reference.length, context, offsets)
            for window in planned:
                sequence = reference.fetch(window.start, window.end)
                if window.strand == "-":
                    sequence = sequence.translate(COMPLEMENT)[::-1]
                if any(base not in "ACGT" for base in sequence):
                    writer.writerow([site.chrom, site.position, window.strand, window.start,
                                     window.end, window.offset, context, site.ref, site.alt,
                                     "", "", "", "", "", "", "", "ambiguous_context"])
                    continue
                probs = predictor.log_probs(sequence, window.target)
                ref = site.ref.translate(COMPLEMENT) if window.strand == "-" else site.ref
                alt = site.alt.translate(COMPLEMENT) if window.strand == "-" else site.alt
                # Expose probabilities in genomic (+) orientation on both strands.
                base_probs = [math.exp(probs[b.translate(COMPLEMENT) if window.strand == "-" else b])
                              for b in "ACGT"]
                writer.writerow([site.chrom, site.position, window.strand, window.start,
                                 window.end, window.offset, context, site.ref, site.alt,
                                 probs[ref], probs[alt], probs[alt] - probs[ref],
                                 *base_probs, "ok"])
                written += 1
            if not planned:
                writer.writerow([site.chrom, site.position, "", "", "", "", context,
                                 site.ref, site.alt, "", "", "", "", "", "", "",
                                 "contig_edge"])
    if written == 0:
        raise ValueError("No Evo scores: check target positions, context length, and ambiguous bases")
    return written
