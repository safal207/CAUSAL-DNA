import importlib.util
from pathlib import Path
import sys
import unittest

try:
    import numpy as np
    import pandas as pd
except ModuleNotFoundError as exc:
    raise unittest.SkipTest(
        "scientific stack is validated in the dedicated Irx3 candidate-context workflow"
    ) from exc

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "gap001_irx3_candidate_context",
    ROOT / "analysis" / "gap001_irx3_candidate_context.py",
)
MOD = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MOD
SPEC.loader.exec_module(MOD)


class Irx3CandidateContextTests(unittest.TestCase):
    def test_candidate_set_is_exact_decisive_shortlist(self):
        self.assertEqual(
            list(MOD.CANDIDATES),
            ["ARID5B", "CUX1", "TET1", "KDM2B", "DNMT1", "MECP2", "CXXC1"],
        )

    def test_summary_measures_irx3_context_without_causal_claim(self):
        irx = np.array([1, 1, 0, 0, 1, 0])
        cand = np.array([1, 0, 1, 0, 1, 0])
        row = MOD.summarize_candidate("S1", "male", "TET1", "Tet1", irx, cand)
        self.assertEqual(row["irx3_positive"], 3)
        self.assertEqual(row["irx3_candidate_double"], 2)
        self.assertAlmostEqual(row["candidate_given_irx3_pct"], 200 / 3)
        self.assertAlmostEqual(row["candidate_given_non_irx3_pct"], 100 / 3)
        self.assertGreater(row["irx3_context_delta_pp"], 0)
        self.assertFalse(MOD.FORBIDDEN_CAUSAL_FIELDS & set(row))

    def test_depletion_is_preserved_not_clipped(self):
        irx = np.array([1, 1, 0, 0])
        cand = np.array([0, 0, 1, 1])
        row = MOD.summarize_candidate("S1", "female", "CUX1", "Cux1", irx, cand)
        self.assertEqual(row["candidate_given_irx3_pct"], 0.0)
        self.assertEqual(row["candidate_given_non_irx3_pct"], 100.0)
        self.assertEqual(row["irx3_context_delta_pp"], -100.0)

    def test_replicate_shortlist_prioritizes_coverage_before_mean(self):
        rows = []
        # Candidate A: weak in every replicate -> coverage 4.
        # Candidate B: very high in only two replicates -> coverage 2.
        for sample, sex in MOD.SAMPLES.items():
            rows.append({
                "sample": sample, "sex": sex, "candidate_tf": "A",
                "mouse_gene": "A", "irx3_candidate_double": 1,
                "candidate_given_irx3_pct": 5.0,
            })
        for i, (sample, sex) in enumerate(MOD.SAMPLES.items()):
            rows.append({
                "sample": sample, "sex": sex, "candidate_tf": "B",
                "mouse_gene": "B", "irx3_candidate_double": 10 if i < 2 else 0,
                "candidate_given_irx3_pct": 90.0 if i < 2 else 0.0,
            })
        sample_df = pd.DataFrame(rows)
        original = MOD.CANDIDATES
        try:
            MOD.CANDIDATES = {"A": "A", "B": "B"}
            pooled_df = pd.DataFrame([
                {"candidate_tf": "A", "candidate_given_irx3_pct": 5.0,
                 "candidate_given_non_irx3_pct": 3.0, "irx3_context_delta_pp": 2.0},
                {"candidate_tf": "B", "candidate_given_irx3_pct": 45.0,
                 "candidate_given_non_irx3_pct": 30.0, "irx3_context_delta_pp": 15.0},
            ])
            ranked = MOD.replicate_summary(sample_df, pooled_df)
            self.assertEqual(ranked.iloc[0]["candidate_tf"], "A")
            self.assertEqual(ranked.iloc[0]["replicates_with_irx3_codetection"], 4)
        finally:
            MOD.CANDIDATES = original

    def test_summary_field_contract_excludes_promotion_terms(self):
        allowed = {
            "dataset", "samples", "candidate_count", "candidates", "detection_rule",
            "ranking_rule", "top_candidate_order", "claim_scope",
        }
        self.assertFalse(MOD.FORBIDDEN_CAUSAL_FIELDS & allowed)


if __name__ == "__main__":
    unittest.main()
