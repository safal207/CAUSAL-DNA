import importlib.util
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "gap001_motif_source_stability",
    ROOT / "analysis" / "gap001_motif_source_stability.py",
)
MOD = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MOD
assert SPEC.loader is not None
SPEC.loader.exec_module(MOD)


def row(tf, delta, site=0.8, crossing=False, record=""):
    return {
        "tf": tf,
        "direction": "gain_on_C" if delta > 0 else "loss_on_C" if delta < 0 else "no_change",
        "delta_normalized": str(delta),
        "abs_delta_normalized": str(abs(delta)),
        "max_site_strength": str(site),
        "threshold_075_crossing": str(crossing).lower(),
        "motif_record_id": record or tf,
    }


class MotifSourceStabilityTests(unittest.TestCase):
    def test_best_representative_is_selected_per_tf(self):
        best = MOD.best_by_tf([
            row("TFX", 0.1, record="rep1"),
            row("TFX", 0.3, record="rep2"),
        ])
        self.assertEqual(best["TFX"]["motif_record_id"], "rep2")
        self.assertEqual(best["TFX"]["representative_count"], 2)

    def test_same_direction_and_min_delta_define_robust_candidate(self):
        comparisons, summary = MOD.compare_sources(
            [row("TFX", 0.2)], [row("TFX", 0.15)], "a", "b", 0.1,
        )
        self.assertTrue(comparisons[0]["direction_consistent"])
        self.assertTrue(comparisons[0]["robust_delta_direction"])
        self.assertEqual(summary["robust_delta_direction_tfs"], 1)

    def test_opposite_direction_is_not_stable_even_if_large(self):
        comparisons, summary = MOD.compare_sources(
            [row("TFX", 0.5, crossing=True)], [row("TFX", -0.5, crossing=True)], "a", "b", 0.1,
        )
        self.assertFalse(comparisons[0]["direction_consistent"])
        self.assertFalse(comparisons[0]["crossing_replicated"])
        self.assertEqual(summary["crossing_replicated_tfs"], 0)

    def test_crossing_mismatch_is_explicit_source_sensitivity(self):
        comparisons, summary = MOD.compare_sources(
            [row("TFX", 0.2, crossing=True)], [row("TFX", 0.2, crossing=False)], "a", "b", 0.1,
        )
        self.assertTrue(comparisons[0]["source_sensitive_crossing"])
        self.assertIn("TFX", summary["source_a_only_crossings"])
        self.assertNotIn("TFX", summary["replicated_crossings"])

    def test_replicated_crossing_requires_same_direction(self):
        comparisons, summary = MOD.compare_sources(
            [row("TFX", -0.2, crossing=True)], [row("TFX", -0.3, crossing=True)], "a", "b", 0.1,
        )
        self.assertTrue(comparisons[0]["crossing_replicated"])
        self.assertEqual(summary["replicated_crossings"], ["TFX"])

    def test_nonoverlap_does_not_inflate_agreement(self):
        _, summary = MOD.compare_sources(
            [row("A", 0.2), row("B", 0.2)],
            [row("B", 0.2), row("C", 0.2)],
            "a", "b", 0.1,
        )
        self.assertEqual(summary["overlap_tfs"], 1)
        self.assertEqual(summary["direction_consistent_tfs"], 1)

    def test_result_has_no_causal_promotion_fields(self):
        comparisons, summary = MOD.compare_sources(
            [row("TFX", 0.2)], [row("TFX", 0.2)], "a", "b", 0.1,
        )
        forbidden = {"cause_found", "causal_status", "materialized", "verified"}
        self.assertFalse(forbidden & set(summary))
        self.assertFalse(forbidden & set(comparisons[0]))


if __name__ == "__main__":
    unittest.main()
