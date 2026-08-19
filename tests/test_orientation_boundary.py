import copy
import json
import unittest
from pathlib import Path

from causal_dna.orientation_boundary import OrientationBoundaryError, OrientationBoundaryGate

ROOT = Path(__file__).resolve().parents[1]
CASE_PATH = ROOT / "cases" / "CDNA-001.orientation-boundary.json"


class OrientationBoundaryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reference = json.loads(CASE_PATH.read_text(encoding="utf-8"))

    def boundary(self, doc, boundary_id="COB_A_OCCUPANCY"):
        return next(b for b in doc["boundaries"] if b["id"] == boundary_id)

    def test_reference_case_is_valid_and_projective(self):
        gate = OrientationBoundaryGate(copy.deepcopy(self.reference))
        gate.assert_valid()
        self.assertFalse(gate.document["boundary_found"])
        self.assertTrue(all(gate.certificate(b["id"]).status == "PROJECTIVE" for b in gate.boundaries))

    def test_temporal_evidence_is_not_enough_for_causal_boundary(self):
        doc = copy.deepcopy(self.reference)
        b = self.boundary(doc)
        b["evidence"]["temporal_precedence"] = ["EXP:T1"]
        b["declared_status"] = "TEMPORALLY_SUPPORTED"
        gate = OrientationBoundaryGate(doc)
        gate.assert_valid()
        self.assertEqual("TEMPORALLY_SUPPORTED", gate.certificate(b["id"]).status)

    def test_intervention_requires_temporal_precedence(self):
        doc = copy.deepcopy(self.reference)
        b = self.boundary(doc)
        b["evidence"]["boundary_crossing_intervention"] = ["EXP:I1"]
        b["declared_status"] = "PROJECTIVE"
        gate = OrientationBoundaryGate(doc)
        self.assertTrue(any("requires temporal precedence" in e for e in gate.errors()))

    def test_crossing_intervention_promotes_only_after_temporal_support(self):
        doc = copy.deepcopy(self.reference)
        b = self.boundary(doc)
        b["evidence"]["temporal_precedence"] = ["EXP:T2"]
        b["evidence"]["boundary_crossing_intervention"] = ["EXP:I2"]
        b["declared_status"] = "INTERVENTION_SUPPORTED"
        gate = OrientationBoundaryGate(doc)
        gate.assert_valid()
        self.assertEqual("INTERVENTION_SUPPORTED", gate.certificate(b["id"]).status)

    def test_reversal_requires_crossing_intervention(self):
        doc = copy.deepcopy(self.reference)
        b = self.boundary(doc)
        b["evidence"]["temporal_precedence"] = ["EXP:T3"]
        b["evidence"]["reversal_rescue"] = ["EXP:R3"]
        b["declared_status"] = "TEMPORALLY_SUPPORTED"
        gate = OrientationBoundaryGate(doc)
        self.assertTrue(any("requires crossing intervention" in e for e in gate.errors()))

    def test_full_chain_is_verified_boundary(self):
        doc = copy.deepcopy(self.reference)
        b = self.boundary(doc)
        b["evidence"]["temporal_precedence"] = ["EXP:T4"]
        b["evidence"]["boundary_crossing_intervention"] = ["EXP:I4"]
        b["evidence"]["reversal_rescue"] = ["EXP:R4"]
        b["evidence"]["independent_verification"] = ["VERIFY:V4"]
        b["declared_status"] = "VERIFIED"
        doc["gap_status"] = "CLOSED"
        doc["boundary_found"] = True
        gate = OrientationBoundaryGate(doc)
        gate.assert_valid()
        cert = gate.certificate(b["id"])
        self.assertEqual("VERIFIED", cert.status)
        self.assertTrue(cert.independent_verification)

    def test_boundary_found_requires_verified_boundary(self):
        doc = copy.deepcopy(self.reference)
        doc["gap_status"] = "CLOSED"
        doc["boundary_found"] = True
        gate = OrientationBoundaryGate(doc)
        self.assertTrue(any("requires at least one VERIFIED boundary" in e for e in gate.errors()))

    def test_open_gap_cannot_claim_boundary_found(self):
        doc = copy.deepcopy(self.reference)
        doc["boundary_found"] = True
        gate = OrientationBoundaryGate(doc)
        self.assertTrue(any("gap_status=OPEN forbids boundary_found=true" in e for e in gate.errors()))

    def test_declared_status_cannot_override_evidence(self):
        doc = copy.deepcopy(self.reference)
        b = self.boundary(doc)
        b["declared_status"] = "VERIFIED"
        gate = OrientationBoundaryGate(doc)
        with self.assertRaises(OrientationBoundaryError):
            gate.assert_valid()


if __name__ == "__main__":
    unittest.main()
