"""GPU-free adapter tests using real CPU tensors and a tiny model fixture."""

import importlib.util
import math
import os
import sys
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from pysrc.sweep_experiments import ntv3


@unittest.skipUnless(importlib.util.find_spec("torch"), "PyTorch not installed")
class NTv3Tests(unittest.TestCase):
    def setUp(self):
        import torch

        self.torch = torch
        self.tokenizer_calls = []
        calls = self.tokenizer_calls

        class Tokenizer:
            mask_token_id = 2
            unk_token_id = 0
            vocab = {"A": 6, "T": 7, "C": 8, "G": 9}

            def __len__(self):
                return 11

            def convert_tokens_to_ids(self, base):
                return self.vocab[base]

            def __call__(self, sequences, **kwargs):
                calls.append((sequences, kwargs))
                return {"input_ids": torch.tensor([
                    [self.vocab[base] for base in sequence] for sequence in sequences])}

        self.tokenizer = Tokenizer()
        self.seen = []

        def forward(tokens):
            self.seen.append(tokens.clone())
            positions = torch.arange(tokens.shape[1], dtype=torch.float32) + 1
            vocabulary = torch.arange(11, dtype=torch.float32)
            return {"logits": (positions[:, None] * vocabulary[None, :] / 100)[None]}

        self.model = Mock(side_effect=forward)
        self.model.to.return_value = self.model
        self.model.eval.return_value = self.model
        self.predictor = ntv3.NTv3Predictor.__new__(ntv3.NTv3Predictor)
        self.predictor.torch = torch
        self.predictor.model = self.model
        self.predictor.tokenizer = self.tokenizer
        self.predictor.context = 128
        self.predictor.device = "cpu"
        self.predictor.base_ids = {base: self.tokenizer.vocab[base] for base in "ACGT"}
        self.predictor.mask_id = self.tokenizer.mask_token_id

    def test_masking_same_position_logits_and_full_vocabulary_normalization(self):
        for target in (0, 37, 127):
            with self.subTest(target=target):
                sequence = "ACGT" * 32
                result = self.predictor.log_probs(sequence, target)
                expected_inputs = self.torch.tensor([[self.tokenizer.vocab[b] for b in sequence]])
                expected_inputs[0, target] = self.tokenizer.mask_token_id
                self.assertTrue(self.torch.equal(self.seen[-1], expected_inputs))
                self.assertEqual(self.tokenizer_calls[-1],
                                 ([sequence], {"return_tensors": "pt", "add_special_tokens": False}))
                expected = (self.torch.arange(11) * ((target + 1) / 100)).log_softmax(-1)
                for base in "ACGT":
                    self.assertAlmostEqual(result[base], expected[self.tokenizer.vocab[base]].item(),
                                           delta=2e-6)
                self.assertLess(sum(math.exp(p) for p in result.values()), 1)

    def test_target_identity_does_not_leak_into_masked_input(self):
        for base in "ACGT":
            sequence = list("ACGT" * 32)
            sequence[37] = base
            self.predictor.log_probs("".join(sequence), 37)
        for tokens in self.seen[1:]:
            self.assertTrue(self.torch.equal(tokens, self.seen[0]))

    def test_invalid_inputs_fail_before_model_execution(self):
        for sequence, target in (("ACGT", 1), ("ACGT" * 32, -1), ("ACGT" * 32, 128),
                                 ("acgt" * 32, 1), ("N" + "A" * 127, 1)):
            with self.subTest(target=target, sequence=sequence[:4]):
                with self.assertRaises(ValueError):
                    self.predictor.log_probs(sequence, target)
        self.model.assert_not_called()

    def test_bad_tokenization_and_output_shapes_fail_explicitly(self):
        with patch.object(type(self.tokenizer), "__call__", return_value={
                "input_ids": self.torch.zeros((1, 129), dtype=self.torch.long)}):
            with self.assertRaisesRegex(ValueError, "single-base"):
                self.predictor.log_probs("ACGT" * 32, 37)
        self.model.side_effect = None
        for shape in ((1, 127, 11), (1, 128, 10), (128, 11)):
            self.model.return_value = {"logits": self.torch.zeros(shape)}
            with self.subTest(shape=shape), self.assertRaisesRegex(ValueError, "logits"):
                self.predictor.log_probs("ACGT" * 32, 37)

    def test_nonfinite_and_backend_errors_propagate(self):
        self.model.side_effect = None
        self.model.return_value = {"logits": self.torch.full((1, 128, 11), float("nan"))}
        with self.assertRaisesRegex(RuntimeError, "nonfinite"):
            self.predictor.log_probs("ACGT" * 32, 37)
        for error in (ImportError("missing kernel"), self.torch.OutOfMemoryError("test OOM")):
            self.model.side_effect = error
            with self.subTest(error=type(error)), self.assertRaises(type(error)):
                self.predictor.log_probs("ACGT" * 32, 37)

    def loaders(self):
        return SimpleNamespace(
            AutoTokenizer=SimpleNamespace(from_pretrained=Mock(return_value=self.tokenizer)),
            AutoModelForMaskedLM=SimpleNamespace(from_pretrained=Mock(return_value=self.model)))

    def test_rocm_loading_pins_model_and_tokenizer_and_probes_requested_context(self):
        loaders = self.loaders()
        with patch.dict(sys.modules, {"transformers": loaders}), \
                patch.object(self.torch.cuda, "is_available", return_value=True), \
                patch.object(self.torch.version, "hip", "7.2"), \
                patch.object(ntv3.NTv3Predictor, "log_probs") as probe:
            predictor = ntv3.NTv3Predictor(ntv3.MODEL, ntv3.REVISION, 8192)
        loaders.AutoTokenizer.from_pretrained.assert_called_once_with(
            ntv3.MODEL, revision=ntv3.REVISION, trust_remote_code=True)
        loaders.AutoModelForMaskedLM.from_pretrained.assert_called_once_with(
            ntv3.MODEL, revision=ntv3.REVISION, trust_remote_code=True, dtype=self.torch.float32)
        self.model.to.assert_called_once_with("cuda:0")
        self.model.eval.assert_called_once_with()
        probe.assert_called_once_with("ACGT" * 2048, 4096)
        self.assertTrue(predictor.rocm)
        with patch.object(ntv3.importlib.metadata, "version", return_value="test"), \
                patch.object(self.torch.cuda, "get_device_name", return_value="AMD fixture"):
            metadata = predictor.metadata()
        self.assertEqual(metadata["model"], ntv3.MODEL)
        self.assertEqual(metadata["revision"], ntv3.REVISION)
        self.assertEqual(metadata["context_length"], 8192)
        self.assertEqual(metadata["precision"], "float32")
        self.assertEqual(metadata["probability_normalization"], "full_vocabulary")
        self.assertEqual(metadata["scoring"], "masked_target_bidirectional_logp_alt_minus_logp_ref")
        self.assertEqual(metadata["device"], "cuda:0")

    def test_loader_errors_do_not_downgrade_or_fallback(self):
        for error in (OSError("gated repository"), ImportError("missing dependency"),
                      self.torch.OutOfMemoryError("test OOM")):
            loaders = self.loaders()
            loaders.AutoModelForMaskedLM.from_pretrained.side_effect = error
            with patch.dict(sys.modules, {"transformers": loaders}), \
                    patch.object(self.torch.cuda, "is_available", return_value=True), \
                    self.subTest(error=type(error)), self.assertRaises(type(error)):
                ntv3.NTv3Predictor(ntv3.MODEL, ntv3.REVISION, 128)
            loaders.AutoModelForMaskedLM.from_pretrained.assert_called_once()

    def test_constructor_rejects_bad_context_model_and_tokenizer_and_missing_gpu(self):
        for context in (127, 129, 1_048_704):
            with self.subTest(context=context), self.assertRaisesRegex(ValueError, "multiple of 128"):
                ntv3.NTv3Predictor(ntv3.MODEL, ntv3.REVISION, context)
        for context in (128, 8192, 1_048_576):
            ntv3.validate_context(context)
        with self.assertRaisesRegex(ValueError, "checkpoint"):
            ntv3.NTv3Predictor("InstaDeepAI/NTv3_tracks", ntv3.REVISION, 128)
        with self.assertRaisesRegex(ValueError, "revision"):
            ntv3.NTv3Predictor(ntv3.MODEL, "", 128)
        with patch.object(self.torch.cuda, "is_available", return_value=False):
            with self.assertRaisesRegex(RuntimeError, "no CPU fallback"):
                ntv3.NTv3Predictor(ntv3.MODEL, ntv3.REVISION, 128)
        for mask_id in (None, 0, 6):
            self.tokenizer.mask_token_id = mask_id
            with patch.dict(sys.modules, {"transformers": self.loaders()}), \
                    patch.object(self.torch.cuda, "is_available", return_value=True), \
                    self.subTest(mask_id=mask_id), self.assertRaisesRegex(ValueError, "tokenizer"):
                ntv3.NTv3Predictor(ntv3.MODEL, ntv3.REVISION, 128)


@unittest.skipUnless(os.environ.get("NTV3_GPU_TEST") == "1",
                     "set NTV3_GPU_TEST=1 for gated-checkpoint GPU validation")
class NTv3GPUTests(unittest.TestCase):
    def test_real_checkpoint_masked_scoring(self):
        context = int(os.environ.get("NTV3_TEST_CONTEXT", "8192"))
        predictor = ntv3.NTv3Predictor(ntv3.MODEL, ntv3.REVISION, context)
        sequence = "ACGT" * (context // 4)
        target = context // 2
        first = predictor.log_probs(sequence, target)
        changed = sequence[:target] + "T" + sequence[target + 1:]
        second = predictor.log_probs(changed, target)
        self.assertEqual(set(first), set("ACGT"))
        self.assertTrue(all(math.isfinite(p) and p <= 0 for p in first.values()))
        self.assertLessEqual(sum(math.exp(p) for p in first.values()), 1 + 1e-6)
        for base in "ACGT":
            self.assertAlmostEqual(first[base], second[base], delta=1e-5)


if __name__ == "__main__":
    unittest.main()
