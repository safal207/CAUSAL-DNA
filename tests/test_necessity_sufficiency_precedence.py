import copy
import json
import unittest
from pathlib import Path

from causal_dna.necessity_sufficiency_gate import NecessitySufficiencyGate

ROOT = Path(__file__).resolve().parents[1]
CASE_PATH = ROOT / "cases" / "CDNA-001.necessity-sufficiency-gate.json"


class NecessitySufficiencyPrecedenceTests(unittest.TestCase):
    def test_attenuation_beats_no_sufficiency_and_remains_partial(self):
        doc = json.loads(CASE_PATH.read_text(encoding="utf-8"))
        doc = copy.deepcopy(doc)
        mediator = next(m for m in doc["mediators"] if m["id"] == "M_ARID5B")
        mediator["necessity_test"]["outcome"] = "attenuated"
        mediator["necessity_test"]["evidence_refs"] = ["EXP:N-PARTIAL"]
        mediator["sufficiency_test"]["outcome"] = "no_change"
        mediator["sufficiency_test"]["evidence_refs"] = ["EXP:S-NOCHANGE"]
        mediator["residual_genotype_effect_after_block"] = "present"
        mediator["declared_role"] = "PARTIAL_MEDIATOR"

        gate = NecessitySufficiencyGate(doc)
        gate.assert_valid()
        self.assertEqual("PARTIAL_MEDIATOR", gate.certificate("M_ARID5B").role)


if __name__ == "__main__":
    unittest.main()
