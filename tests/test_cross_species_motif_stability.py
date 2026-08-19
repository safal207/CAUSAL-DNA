import csv
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "gap001_cross_species_motif_stability",
    ROOT / "analysis" / "gap001_cross_species_motif_stability.py",
)
MOD = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MOD
SPEC.loader.exec_module(MOD)


def row(tf, rid, delta, max_site=0.8, crossing=False):
    return {
        "tf": tf,
        "motif_record_id": rid,
        "direction": "gain_on_C" if delta > 0 else "loss_on_C" if delta < 0 else "no_change",
        "delta_normalized": str(delta),
        "abs_delta_normalized": str(abs(delta)),
        "max_site_strength": str(max_site),
        "threshold_075_crossing": str(crossing),
    }


class CrossSpeciesMotifStabilityTests(unittest.TestCase):
    def test_same_record_direction_is_conserved(self):
        human = [row("CXXC4", "CXXC4__rep1", 0.30)]
        mouse = [row("CXXC4", "CXXC4__rep1", 0.21)]
        rows, summary = MOD.compare_species(human, mouse, 0.10)
        self.assertTrue(rows[0]["direction_conserved"])
        self.assertTrue(rows[0]["robust_cross_species_delta"])
        self.assertEqual(summary["robust_cross_species_tfs"], 1)

    def test_reversal_demotes_candidate(self):
        human = [row("TFX", "TFX__rep1", 0.30)]
        mouse = [row("TFX", "TFX__rep1", -0.25)]
        rows, _ = MOD.compare_species(human, mouse, 0.10)
        self.assertFalse(rows[0]["direction_conserved"])
        self.assertFalse(rows[0]["robust_cross_species_delta"])

    def test_human_representative_is_frozen_before_mouse_evaluation(self):
        human = [
            row("TFX", "TFX__rep1", 0.35),
            row("TFX", "TFX__rep2", 0.15),
        ]
        mouse = [
            row("TFX", "TFX__rep1", 0.02),
            row("TFX", "TFX__rep2", 0.40),
        ]
        rows, _ = MOD.compare_species(human, mouse, 0.10)
        self.assertEqual(rows[0]["motif_record_id"], "TFX__rep1")
        self.assertFalse(rows[0]["robust_cross_species_delta"])

    def test_record_sets_must_match(self):
        human = [row("TFX", "TFX__rep1", 0.2)]
        mouse = [row("TFX", "TFX__rep2", 0.2)]
        with self.assertRaises(ValueError):
            MOD.compare_species(human, mouse, 0.10)

    def test_tf_identity_must_match_for_same_record(self):
        human = [row("TFX", "shared", 0.2)]
        mouse = [row("TFY", "shared", 0.2)]
        with self.assertRaises(ValueError):
            MOD.compare_species(human, mouse, 0.10)

    def test_crossing_requires_direction_conservation(self):
        human = [row("TFX", "TFX__rep1", 0.2, crossing=True)]
        mouse = [row("TFX", "TFX__rep1", -0.2, crossing=True)]
        rows, summary = MOD.compare_species(human, mouse, 0.10)
        self.assertFalse(rows[0]["crossing_conserved"])
        self.assertEqual(summary["conserved_crossing_tfs"], 0)

    def test_outputs_do_not_claim_causality(self):
        human = [row("TFX", "TFX__rep1", 0.2)]
        mouse = [row("TFX", "TFX__rep1", 0.2)]
        rows, summary = MOD.compare_species(human, mouse, 0.10)
        forbidden = {"cause_found", "causal_status", "materialized", "verified"}
        self.assertFalse(forbidden & set(summary))
        self.assertFalse(forbidden & set(rows[0]))


if __name__ == "__main__":
    unittest.main()
