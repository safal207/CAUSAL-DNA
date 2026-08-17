import copy
import json
import unittest
from pathlib import Path

from causal_dna.temporal_lattice import TemporalLattice, TemporalLatticeError

ROOT = Path(__file__).resolve().parents[1]
PATH = ROOT / "cases" / "CDNA-001-rs1421085.temporal-lattice.json"


class TemporalLatticeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.document = json.loads(PATH.read_text(encoding="utf-8"))

    def lattice(self):
        return TemporalLattice(copy.deepcopy(self.document))

    def test_reference_lattice_is_valid(self):
        lattice = self.lattice()
        lattice.assert_valid()
        self.assertEqual([], lattice.errors())

    def test_future_material_fact_is_forbidden(self):
        doc = copy.deepcopy(self.document)
        state = next(s for s in doc["states"] if s["id"] == "T2_DECISIVE_EXPERIMENT")
        state["space"] = "material"
        lattice = TemporalLattice(doc)
        self.assertTrue(any("future state cannot be material fact" in e for e in lattice.errors()))

    def test_projective_cannot_bypass_bardo(self):
        doc = copy.deepcopy(self.document)
        doc["transitions"].append({
            "id": "TL_BAD_BYPASS",
            "from": "T1_GAP_MODEL",
            "to": "T0_IRX3_PH_UP",
            "relation": "bounds",
            "evidence_refs": ["fake:evidence"]
        })
        lattice = TemporalLattice(doc)
        self.assertTrue(any("projective→material bypass forbidden" in e for e in lattice.errors()))

    def test_materialization_requires_bardo_and_provenance(self):
        doc = copy.deepcopy(self.document)
        doc["transitions"].append({
            "id": "TL_BAD_MATERIALIZE",
            "from": "T1_GAP_MODEL",
            "to": "T0_IRX3_PH_UP",
            "relation": "materializes",
            "evidence_refs": []
        })
        lattice = TemporalLattice(doc)
        errors = lattice.errors()
        self.assertTrue(any("materializes must be Bardo→Material" in e for e in errors))
        self.assertTrue(any("materialization requires evidence_refs" in e for e in errors))

    def test_future_frontier_contains_decisive_experiment_and_gate(self):
        ids = {s["id"] for s in self.lattice().future_frontier()}
        self.assertIn("T2_DECISIVE_EXPERIMENT", ids)
        self.assertIn("T2_RESOLUTION_GATE", ids)

    def test_materialized_history_is_evidence_backed(self):
        history = self.lattice().materialized_history()
        self.assertGreaterEqual(len(history), 3)
        self.assertTrue(all(s["evidence_refs"] for s in history))

    def test_observer_views_are_separated(self):
        views = self.lattice().observer_views()
        self.assertIn("model", views)
        self.assertIn("experiment", views)
        self.assertIn("independent_verifier", views)


if __name__ == "__main__":
    unittest.main()
