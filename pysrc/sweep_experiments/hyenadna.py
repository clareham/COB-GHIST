"""Causal single-base scoring with the public HyenaDNA large checkpoint."""

import importlib.metadata
import logging
import math

LOG = logging.getLogger(__name__)
MODEL = "LongSafari/hyenadna-large-1m-seqlen-hf"
REVISION = "0a629abf9c7f85b4ec9aa6a1aefa3adcf1907446"
MAX_CONTEXT = 1_000_000


def validate_context(context: int, offsets: list[int]) -> None:
    if not 2 <= context <= MAX_CONTEXT:
        raise ValueError(f"HyenaDNA context must be between 2 and {MAX_CONTEXT:,} bases")
    if not offsets or any(not 1 <= offset < context for offset in offsets):
        raise ValueError("HyenaDNA offsets must satisfy 1 <= offset < context length")


class HyenaDNAPredictor:
    def __init__(self, revision: str, context: int, offsets: list[int]):
        validate_context(context, offsets)
        if not revision:
            raise ValueError("HyenaDNA model revision must not be empty")

        import torch

        if not torch.cuda.is_available():
            raise RuntimeError("HyenaDNA requires an available CUDA/ROCm GPU; no CPU fallback")
        from transformers import AutoModelForCausalLM, AutoTokenizer

        self.torch = torch
        self.device = "cuda:0"
        self.context = context
        self.max_prefix = max(offsets)
        self.revision = revision
        LOG.info("Loading %s at revision %s (custom Hub code, float32)", MODEL, revision)
        self.tokenizer = AutoTokenizer.from_pretrained(
            MODEL, revision=revision, trust_remote_code=True)
        self.base_ids = {base: self.tokenizer.convert_tokens_to_ids(base) for base in "ACGT"}
        if (any(index is None or index < 0 or index == self.tokenizer.unk_token_id
                for index in self.base_ids.values()) or len(set(self.base_ids.values())) != 4):
            raise ValueError("HyenaDNA tokenizer must expose distinct A/C/G/T tokens")
        self.model = AutoModelForCausalLM.from_pretrained(
            MODEL, revision=revision, trust_remote_code=True, dtype=torch.float32)
        self.model.to(self.device).eval()
        self.log_probs(("ACGT" * ((context + 3) // 4))[:context], self.max_prefix)
        LOG.info("Loaded %s and passed a %d-base causal prefix probe", MODEL, self.max_prefix)

    def log_probs(self, sequence: str, target: int) -> dict[str, float]:
        if len(sequence) != self.context or not 1 <= target <= self.max_prefix:
            raise ValueError("Expected the configured HyenaDNA window length and a probed target offset")
        prefix = sequence[:target]
        if any(base not in "ACGT" for base in prefix):
            raise ValueError("HyenaDNA scoring requires an uppercase A/C/G/T prefix")
        # A causal pass needs only the prefix. Neither the target nor its suffix enters the model.
        inputs = self.tokenizer(prefix, return_tensors="pt", add_special_tokens=False)["input_ids"]
        if tuple(inputs.shape) != (1, target):
            raise ValueError("Expected HyenaDNA single-base tokenization without special tokens")
        with self.torch.inference_mode():
            logits = self.model(input_ids=inputs.to(self.device)).logits
            if (logits.ndim != 3 or tuple(logits.shape[:2]) != (1, target)
                    or logits.shape[2] <= max(self.base_ids.values())):
                raise ValueError("Expected HyenaDNA logits with shape (1, prefix_length, vocabulary)")
            # The LM head includes padded vocabulary channels beyond the tokenizer's vocabulary.
            probabilities = logits[0, target - 1].float().log_softmax(dim=-1)
            result = {base: float(probabilities[index].item())
                      for base, index in self.base_ids.items()}
        if any(not math.isfinite(value) for value in result.values()):
            raise RuntimeError("HyenaDNA produced nonfinite log probabilities")
        return result

    def metadata(self) -> dict:
        return {"model": MODEL, "revision": self.revision,
                "transformers_version": importlib.metadata.version("transformers"),
                "torch_version": self.torch.__version__, "rocm": self.torch.version.hip is not None,
                "device": self.device, "device_name": self.torch.cuda.get_device_name(0),
                "hip_version": self.torch.version.hip, "cuda_version": self.torch.version.cuda,
                "precision": "float32", "context_length": self.context,
                "max_prefix_length": self.max_prefix,
                "scoring": "causal_prefix_next_base_logp_alt_minus_logp_ref",
                "probability_normalization": "full_lm_head_including_padding"}
