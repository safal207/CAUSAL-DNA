import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest

try:
    import numpy as np
except ModuleNotFoundError:  # lightweight causal validator intentionally omits scientific stack
    np = None

ROOT = Path(__file__).resolve().parents[1]
MOD = None
if np is not None:
    SPEC = importlib.util.spec_from_file_location(
        "gap001_codebook_motif_switch",
        ROOT / "analysis" / "gap001_codebook_motif_switch.py",
    )
    MOD = importlib.util.module_from_spec(SPEC)
    assert SPEC.loader is not None
    # dataclasses resolves postponed annotations through sys.modules on Python 3.12.
    # Register the dynamically loaded analysis module before executing it.
    sys.modules[SPEC.name] = MOD
    SPEC.loader.exec_module(MOD)


@unittest.skipUnless(np is not None, "scientific stack is validated in the dedicated Codebook workflow")
class CodebookMotifSwitchTests(unittest.TestCase):
    def test_same_placement_detects_gain_on_c(self):
        # Three-base motif strongly prefers C at the middle position.
        pwm = np.array([
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 2.0, 0.0, -1.0],
            [1.0, 0.0, 0.0, 0.0],
        ])
        ref = "AATAG"
        alt = "AACAG"
        result = MOD.best_variant_overlapping_switch("TFX", pwm, ref, alt, 2)
        self.assertEqual(result.direction, "gain_on_C")
        self.assertGreater(result.delta_normalized, 0)

    def test_only_variant_overlapping_windows_are_considered(self):
        pwm = np.array([
            [3.0, 0.0, 0.0, 0.0],
            [3.0, 0.0, 0.0, 0.0],
        ])
        ref = "AATAA"
        alt = "AACAA"
        result = MOD.best_variant_overlapping_switch("TFX", pwm, ref, alt, 2)
        # An AA motif exists away from the variant, but it must not create a fake allele effect.
        self.assertEqual(result.abs_delta_normalized, 0.0)

    def test_reverse_complement(self):
        self.assertEqual(MOD.reverse_complement("ACGTTA"), "TAACGT")

    def test_pwm_parser_requires_four_columns(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "bad.pwm"
            path.write_text(">BAD\n1 2 3\n4 5 6\n")
            with self.assertRaises(ValueError):
                MOD.load_pwm(path)

    def test_mouse_symbol_candidates_for_znf(self):
        self.assertEqual(MOD.mouse_symbol_candidates("ZNF445"), ["Znf445", "Zfp445"])


if __name__ == "__main__":
    unittest.main()
