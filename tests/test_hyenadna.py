"""GPU-free tests of the HyenaDNA causal scoring contract."""

import importlib.util
import sys
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from pysrc.sweep_experiments import hyenadna


@unittest.skipUnless(importlib.util.find_spec("torch"), "PyTorch not installed")
class HyenaDNATests(unittest.TestCase):
    def test_pinned_load_and_probe_largest_prefix_on_gpu(self):
        import torch

        base_ids = dict(zip("ACGT", (7, 8, 9, 10)))
        tokens = Mock(shape=(1, 5))
        tokenizer = Mock(return_value={"input_ids": tokens}, unk_token_id=6)
        tokenizer.convert_tokens_to_ids.side_effect = base_ids.__getitem__
        logits = torch.arange(80, dtype=torch.float32).reshape(1, 5, 16) / 100
        logits[0, 4, 15] = 4
        model = Mock(return_value=SimpleNamespace(logits=logits))
        transformers = SimpleNamespace(
            AutoTokenizer=SimpleNamespace(from_pretrained=Mock(return_value=tokenizer)),
            AutoModelForCausalLM=SimpleNamespace(from_pretrained=Mock(return_value=model)))
        with patch.dict(sys.modules, {"transformers": transformers}), \
                patch.object(torch.cuda, "is_available", return_value=True):
            predictor = hyenadna.HyenaDNAPredictor(hyenadna.REVISION, 8, [5, 2])
        opts = {"revision": hyenadna.REVISION, "trust_remote_code": True}
        transformers.AutoTokenizer.from_pretrained.assert_called_once_with(hyenadna.MODEL, **opts)
        transformers.AutoModelForCausalLM.from_pretrained.assert_called_once_with(
            hyenadna.MODEL, dtype=torch.float32, **opts)
        model.to.assert_called_once_with("cuda:0")
        tokens.to.assert_called_once_with("cuda:0")
        tokenizer.assert_called_once_with("ACGTA", return_tensors="pt", add_special_tokens=False)
        model.assert_called_once_with(input_ids=tokens.to.return_value)

        expected = logits[0, 4].log_softmax(-1)
        for sequence in ("ACGTAAAA", "ACGTAGGG"):
            result = predictor.log_probs(sequence, 5)
            tokenizer.assert_called_with("ACGTA", return_tensors="pt", add_special_tokens=False)
            self.assertEqual(result, {base: expected[i].item() for base, i in base_ids.items()})

        with self.assertRaisesRegex(ValueError, "probed target"):
            predictor.log_probs("ACGTAAAA", 0)
        with self.assertRaisesRegex(ValueError, "probed target"):
            predictor.log_probs("ACGTAAAA", 6)
        logits[0, 4, 15] = torch.nan
        with self.assertRaisesRegex(RuntimeError, "nonfinite"):
            predictor.log_probs("ACGTAAAA", 5)

    def test_no_gpu_fails_before_load(self):
        import torch

        with patch.object(torch.cuda, "is_available", return_value=False), \
                self.assertRaisesRegex(RuntimeError, "no CPU fallback"):
            hyenadna.HyenaDNAPredictor(hyenadna.REVISION, 8, [1])
