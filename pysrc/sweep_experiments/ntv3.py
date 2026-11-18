"""NTv3 masked-base scoring and strand-aware context-offset experiments."""

from dataclasses import dataclass
import importlib.metadata
import logging
import math
from pathlib import Path
from typing import Protocol
import csv

from .inputs import Reference, Variant

LOG = logging.getLogger(__name__)
MODEL = "InstaDeepAI/NTv3_pretrained"
REVISION = "fd58207244d828957e633c591f9babd8bb2b970e"
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def validate_context(context: int) -> None:
    if not 128 <= context <= 1_048_576 or context % 128:
        raise ValueError("NTv3 context must be a multiple of 128 between 128 and 1,048,576 bases")


@dataclass(frozen=True)
class Window:
    start: int
    end: int
    target: int
    strand: str
    offset: int


def windows(position: int, length: int, context: int,
            offsets: list[int]) -> list[Window]:
    if context < 2 or not offsets or len(set(offsets)) != len(offsets):
        raise ValueError("Context must be >= 2 and offsets must be nonempty and unique")
    result = []
    for offset in offsets:
        if not 0 <= offset < context:
            raise ValueError("Offsets must satisfy 0 <= offset < context length")
        for strand in ("+", "-"):
            target = offset if strand == "+" else context - 1 - offset
            start = position - 1 - target
            if start >= 0 and start + context <= length:
                result.append(Window(start, start + context, offset, strand, offset))
    return result


class Predictor(Protocol):
    def log_probs(self, sequence: str, target: int) -> dict[str, float]: ...


class NTv3Predictor:
    def __init__(self, model: str, revision: str, context: int):
        validate_context(context)
        if model != MODEL:
            raise ValueError(f"Expected the NTv3 masked-language checkpoint {MODEL}")
        if not revision:
            raise ValueError("NTv3 model revision must not be empty")
        import torch

        self.torch = torch
        if not torch.cuda.is_available():
            raise RuntimeError("NTv3 requires an available CUDA/ROCm GPU; no CPU fallback")
        self.rocm = torch.version.hip is not None
        self.device = "cuda:0"
        self.context = context
        self.selected = model
        self.revision = revision
        from transformers import AutoModelForMaskedLM, AutoTokenizer

        LOG.info("Loading %s at revision %s (custom Hub code, float32)", model, revision)
        self.tokenizer = AutoTokenizer.from_pretrained(
            model, revision=revision, trust_remote_code=True)
        self.base_ids = {base: self.tokenizer.convert_tokens_to_ids(base) for base in "ACGT"}
        self.mask_id = self.tokenizer.mask_token_id
        token_ids = [*self.base_ids.values(), self.mask_id]
        if (any(value is None or value == self.tokenizer.unk_token_id for value in token_ids)
                or len(set(token_ids)) != 5):
            raise ValueError("NTv3 tokenizer must expose distinct A/C/G/T and mask tokens")
        self.model = AutoModelForMaskedLM.from_pretrained(
            model, revision=revision, trust_remote_code=True, dtype=torch.float32)
        self.model.to(self.device).eval()
        self.log_probs(("ACGT" * ((context + 3) // 4))[:context], context // 2)
        LOG.info("Loaded %s and passed a %d-base masked context probe", model, context)

    def log_probs(self, sequence: str, target: int) -> dict[str, float]:
        if len(sequence) != self.context or not 0 <= target < len(sequence):
            raise ValueError("Expected the configured NTv3 context length and an in-range target")
        if any(base not in "ACGT" for base in sequence):
            raise ValueError("NTv3 scoring requires uppercase A/C/G/T sequence")
        inputs = self.tokenizer(
            [sequence], return_tensors="pt", add_special_tokens=False)["input_ids"]
        if tuple(inputs.shape) != (1, len(sequence)):
            raise ValueError("Expected NTv3 single-base tokenization without special tokens")
        inputs = inputs.to(self.device)
        inputs[0, target] = self.mask_id
        with self.torch.inference_mode():
            logits = self.model(inputs)["logits"]
            if (logits.ndim != 3 or tuple(logits.shape[:2]) != (1, len(sequence))
                    or logits.shape[2] != len(self.tokenizer)):
                raise ValueError("Expected NTv3 logits with shape (1, context_length, vocabulary)")
            probabilities = logits[0, target].float().log_softmax(dim=-1)
            result = {base: float(probabilities[index].item())
                      for base, index in self.base_ids.items()}
        if any(not math.isfinite(value) for value in result.values()):
            raise RuntimeError("NTv3 produced nonfinite log probabilities")
        return result

    def metadata(self) -> dict:
        return {"model": self.selected, "revision": self.revision,
                "transformers_version": importlib.metadata.version("transformers"),
                "torch_version": self.torch.__version__, "rocm": self.rocm,
                "device": self.device, "device_name": self.torch.cuda.get_device_name(0),
                "hip_version": self.torch.version.hip, "cuda_version": self.torch.version.cuda,
                "precision": "float32", "context_length": self.context,
                "scoring": "masked_target_bidirectional_logp_alt_minus_logp_ref",
                "probability_normalization": "full_vocabulary"}


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
        raise ValueError("No NTv3 scores: check target positions, context length, and ambiguous bases")
    return written
