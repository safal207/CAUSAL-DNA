import copy
import json
import unittest
from pathlib import Path

from jsonschema import Draft202012Validator

from causal_dna.mediation_gate import MediationGate, MediationGateError

ROOT = Path(__file__).resolve().parents[1]
CASE_PATH = ROOT / "cases" / "CDNA-001.mediation-gate.json"
SCHEMA_PATH = ROOT / "schemas" / "causal-mediation-gate.schema.json"


class MediationGateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.document = json.loads(CASE_PATH.read_text(encoding="utf-8"))
        cls.schema = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))

    def gate(self, doc=None):
        return MediationGate(copy.deepcopy(doc if doc is not None else self.document))

    def candidate(self, doc, mediator_id="M_ARID5B"):
        return next(m for m in doc["mediators"] if m["id"] == mediator_id)

    def test_reference_document_matches_schema(self):
        errors = list(Draft202012Validator(self.schema).iter_errors(self.document))
        self.assertEqual([], errors)

    def test_reference_gate_is_valid_and_open(self):
        gate = self.gate()
        gate.assert_valid()
        self.assertFalse(gate.summary()["cause_found"])
        self.assertEqual([], gate.summary()["materialized"])

    def test_current_candidates_are_not_silently_promoted(self):
        gate = self.gate()
        for cert in gate.certificates():
            self.assertEqual("PROJECTIVE", cert.status)
            self.assertEqual("projective", cert.space)
            self.assertFalse(cert.materialized)

    def test_association_only_stays_in_bardo(self):
        doc = copy.deepcopy(self.document)
        c = self.candidate(doc)
        c["evidence"]["genotype_to_mediator"] = ["evidence:allele-mediator"]
        c["declared_status"] = "BARDO_ASSOCIATED"
        c["declared_space"] = "bardo"
        gate = self.gate(doc)
        gate.assert_valid()
        self.assertEqual("BARDO_ASSOCIATED", gate.certificate("M_ARID5B").status)

    def test_intervention_without_rescue_stays_in_bardo(self):
        doc = copy.deepcopy(self.document)
        c = self.candidate(doc)
        c["evidence"]["genotype_to_mediator"] = ["evidence:allele-mediator"]
        c["evidence"]["intervention_on_mediator"] = ["evidence:intervention"]
        c["declared_status"] = "BARDO_INTERVENTION_SUPPORTED"
        c["declared_space"] = "bardo"
        gate = self.gate(doc)
        gate.assert_valid()
        cert = gate.certificate("M_ARID5B")
        self.assertEqual("BARDO_INTERVENTION_SUPPORTED", cert.status)
        self.assertFalse(cert.materialized)

    def test_intervention_and_rescue_are_not_final_without_independent_verification(self):
        doc = copy.deepcopy(self.document)
        c = self.candidate(doc)
        c["evidence"]["genotype_to_mediator"] = ["evidence:allele-mediator"]
        c["evidence"]["intervention_on_mediator"] = ["evidence:intervention"]
        c["evidence"]["rescue"] = ["evidence:rescue"]
        c["declared_status"] = "BARDO_RESCUE_SUPPORTED"
        c["declared_space"] = "bardo"
        gate = self.gate(doc)
        gate.assert_valid()
        cert = gate.certificate("M_ARID5B")
        self.assertEqual("BARDO_RESCUE_SUPPORTED", cert.status)
        self.assertFalse(cert.materialized)
        self.assertIn("independent_verification", cert.missing_requirements)

    def test_full_chain_materializes_only_after_independent_verification(self):
        doc = copy.deepcopy(self.document)
        doc["gap_status"] = "RESOLVED"
        c = self.candidate(doc)
        c["evidence"]["genotype_to_mediator"] = ["evidence:allele-mediator"]
        c["evidence"]["intervention_on_mediator"] = ["evidence:intervention"]
        c["evidence"]["rescue"] = ["evidence:rescue"]
        c["evidence"]["independent_verification"] = ["evidence:verification"]
        c["declared_status"] = "MATERIAL"
        c["declared_space"] = "material"
        gate = self.gate(doc)
        gate.assert_valid()
        cert = gate.certificate("M_ARID5B")
        self.assertTrue(cert.materialized)
        self.assertEqual("MATERIAL", cert.status)

    def test_rescue_without_intervention_is_invalid(self):
        doc = copy.deepcopy(self.document)
        c = self.candidate(doc)
        c["evidence"]["rescue"] = ["evidence:rescue"]
        gate = self.gate(doc)
        errors = gate.errors()
        self.assertTrue(any("rescue requires intervention evidence" in e for e in errors))

    def test_manual_material_promotion_is_rejected(self):
        doc = copy.deepcopy(self.document)
        c = self.candidate(doc)
        c["declared_status"] = "MATERIAL"
        c["declared_space"] = "material"
        gate = self.gate(doc)
        errors = gate.errors()
        self.assertTrue(any("evidence computes PROJECTIVE" in e for e in errors))
        self.assertTrue(any("evidence computes 'projective'" in e for e in errors))

    def test_open_gap_forbids_material_certificate(self):
        doc = copy.deepcopy(self.document)
        c = self.candidate(doc)
        c["evidence"]["genotype_to_mediator"] = ["evidence:allele-mediator"]
        c["evidence"]["intervention_on_mediator"] = ["evidence:intervention"]
        c["evidence"]["rescue"] = ["evidence:rescue"]
        c["evidence"]["independent_verification"] = ["evidence:verification"]
        c["declared_status"] = "MATERIAL"
        c["declared_space"] = "material"
        gate = self.gate(doc)
        errors = gate.errors()
        self.assertTrue(any("gap_status=OPEN forbids a MATERIAL mediator certificate" in e for e in errors))

    def test_cause_found_cannot_be_claimed_without_material_mediator(self):
        doc = copy.deepcopy(self.document)
        doc["gap_status"] = "RESOLVED"
        doc["cause_found"] = True
        gate = self.gate(doc)
        with self.assertRaises(MediationGateError):
            gate.assert_valid()


if __name__ == "__main__":
    unittest.main()
