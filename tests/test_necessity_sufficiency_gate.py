import copy
import json
import unittest
from pathlib import Path

from causal_dna.necessity_sufficiency_gate import (
    NecessitySufficiencyError,
    NecessitySufficiencyGate,
)

ROOT = Path(__file__).resolve().parents[1]
CASE_PATH = ROOT / "cases" / "CDNA-001.necessity-sufficiency-gate.json"


class NecessitySufficiencyGateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.document = json.loads(CASE_PATH.read_text(encoding="utf-8"))

    def gate(self):
        return NecessitySufficiencyGate(copy.deepcopy(self.document))

    def candidate(self, doc, mediator_id="M_ARID5B"):
        return next(m for m in doc["mediators"] if m["id"] == mediator_id)

    def test_reference_case_is_valid_and_unresolved(self):
        gate = self.gate()
        gate.assert_valid()
        self.assertEqual([], gate.errors())
        self.assertTrue(all(c.role == "UNRESOLVED" for c in gate.certificates()))
        self.assertFalse(gate.document["causal_role_found"])

    def test_abolished_plus_recapitulated_is_necessary_and_sufficient(self):
        doc = copy.deepcopy(self.document)
        m = self.candidate(doc)
        m["necessity_test"]["outcome"] = "abolished"
        m["necessity_test"]["evidence_refs"] = ["EXP:N1"]
        m["sufficiency_test"]["outcome"] = "recapitulated"
        m["sufficiency_test"]["evidence_refs"] = ["EXP:S1"]
        m["residual_genotype_effect_after_block"] = "absent"
        m["declared_role"] = "NECESSARY_AND_SUFFICIENT"
        gate = NecessitySufficiencyGate(doc)
        gate.assert_valid()
        self.assertEqual("NECESSARY_AND_SUFFICIENT", gate.certificate("M_ARID5B").role)
        self.assertFalse(gate.certificate("M_ARID5B").independently_verified)

    def test_attenuation_is_partial_mediation(self):
        doc = copy.deepcopy(self.document)
        m = self.candidate(doc)
        m["necessity_test"]["outcome"] = "attenuated"
        m["necessity_test"]["evidence_refs"] = ["EXP:N2"]
        m["residual_genotype_effect_after_block"] = "present"
        m["declared_role"] = "PARTIAL_MEDIATOR"
        gate = NecessitySufficiencyGate(doc)
        gate.assert_valid()
        self.assertEqual("PARTIAL_MEDIATOR", gate.certificate("M_ARID5B").role)

    def test_sufficient_but_not_necessary_detects_redundant_or_parallel_route(self):
        doc = copy.deepcopy(self.document)
        m = self.candidate(doc)
        m["necessity_test"]["outcome"] = "no_change"
        m["necessity_test"]["evidence_refs"] = ["EXP:N3"]
        m["sufficiency_test"]["outcome"] = "recapitulated"
        m["sufficiency_test"]["evidence_refs"] = ["EXP:S3"]
        m["residual_genotype_effect_after_block"] = "present"
        m["declared_role"] = "SUFFICIENT_NOT_NECESSARY"
        gate = NecessitySufficiencyGate(doc)
        gate.assert_valid()
        self.assertEqual("SUFFICIENT_NOT_NECESSARY", gate.certificate("M_ARID5B").role)

    def test_no_change_in_both_tests_is_non_load_bearing(self):
        doc = copy.deepcopy(self.document)
        m = self.candidate(doc)
        m["necessity_test"]["outcome"] = "no_change"
        m["necessity_test"]["evidence_refs"] = ["EXP:N4"]
        m["sufficiency_test"]["outcome"] = "no_change"
        m["sufficiency_test"]["evidence_refs"] = ["EXP:S4"]
        m["residual_genotype_effect_after_block"] = "present"
        m["declared_role"] = "NON_LOAD_BEARING"
        gate = NecessitySufficiencyGate(doc)
        gate.assert_valid()
        self.assertEqual("NON_LOAD_BEARING", gate.certificate("M_ARID5B").role)

    def test_tested_outcome_requires_provenance(self):
        doc = copy.deepcopy(self.document)
        m = self.candidate(doc)
        m["necessity_test"]["outcome"] = "abolished"
        m["declared_role"] = "NECESSITY_SUPPORTED"
        gate = NecessitySufficiencyGate(doc)
        errors = gate.errors()
        self.assertTrue(any("tested necessity_test requires evidence_refs" in e for e in errors))

    def test_residual_effect_requires_necessity_test(self):
        doc = copy.deepcopy(self.document)
        m = self.candidate(doc)
        m["residual_genotype_effect_after_block"] = "present"
        gate = NecessitySufficiencyGate(doc)
        self.assertTrue(any("residual genotype effect requires" in e for e in gate.errors()))

    def test_abolished_cannot_leave_residual_genotype_effect(self):
        doc = copy.deepcopy(self.document)
        m = self.candidate(doc)
        m["necessity_test"]["outcome"] = "abolished"
        m["necessity_test"]["evidence_refs"] = ["EXP:N5"]
        m["residual_genotype_effect_after_block"] = "present"
        m["declared_role"] = "NECESSITY_SUPPORTED"
        gate = NecessitySufficiencyGate(doc)
        self.assertTrue(any("abolished conflicts" in e for e in gate.errors()))

    def test_declared_role_cannot_override_evidence(self):
        doc = copy.deepcopy(self.document)
        m = self.candidate(doc)
        m["declared_role"] = "NECESSARY_AND_SUFFICIENT"
        gate = NecessitySufficiencyGate(doc)
        with self.assertRaises(NecessitySufficiencyError):
            gate.assert_valid()

    def test_independent_verification_upgrades_certificate_not_role_logic(self):
        doc = copy.deepcopy(self.document)
        m = self.candidate(doc)
        m["necessity_test"]["outcome"] = "abolished"
        m["necessity_test"]["evidence_refs"] = ["EXP:N6"]
        m["sufficiency_test"]["outcome"] = "recapitulated"
        m["sufficiency_test"]["evidence_refs"] = ["EXP:S6"]
        m["residual_genotype_effect_after_block"] = "absent"
        m["declared_role"] = "NECESSARY_AND_SUFFICIENT"
        m["independent_verification"]["evidence_refs"] = ["VERIFY:V1"]
        gate = NecessitySufficiencyGate(doc)
        gate.assert_valid()
        cert = gate.certificate("M_ARID5B")
        self.assertTrue(cert.independently_verified)
        self.assertEqual("VERIFIED_NECESSARY_AND_SUFFICIENT", cert.verified_role)

    def test_open_gap_cannot_claim_causal_role_found(self):
        doc = copy.deepcopy(self.document)
        doc["causal_role_found"] = True
        gate = NecessitySufficiencyGate(doc)
        self.assertTrue(any("gap_status=OPEN forbids causal_role_found=true" in e for e in gate.errors()))


if __name__ == "__main__":
    unittest.main()
