import copy
import json
import unittest
from pathlib import Path

from causal_dna.processor import CausalProcessor, CausalProcessorError

ROOT = Path(__file__).resolve().parents[1]
EVENTS = ROOT / "cases" / "CDNA-001.events.json"


class CausalProcessorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.events = json.loads(EVENTS.read_text(encoding="utf-8"))

    def test_reference_stream_is_valid(self):
        p = CausalProcessor(self.events)
        self.assertEqual(1, p.generation)

    def test_projection_is_rebuildable(self):
        a = CausalProcessor(self.events).project()
        b = CausalProcessor(copy.deepcopy(self.events)).project()
        self.assertEqual(a, b)
        self.assertIn("H3_ACCESSIBILITY", a.open_hypotheses)
        self.assertIn("LHX6_BROAD_ADULT_PH", a.rejected_hypotheses)

    def test_existing_event_cannot_be_rewritten(self):
        p = CausalProcessor(self.events)
        with self.assertRaises(CausalProcessorError):
            p.append(copy.deepcopy(self.events[0]))

    def test_historical_generation_cannot_be_appended(self):
        p = CausalProcessor(self.events)
        with self.assertRaises(CausalProcessorError):
            p.append({"event_id":"EV-X","generation":0,"event_type":"model_update","subject":"x","observer":"model","evidence_refs":[]})

    def test_verification_requires_independent_verifier(self):
        bad = copy.deepcopy(self.events)
        bad.append({"event_id":"EV-X","generation":1,"event_type":"verification","subject":"H3_ACCESSIBILITY","observer":"model","evidence_refs":["x"]})
        with self.assertRaises(CausalProcessorError):
            CausalProcessor(bad)

    def test_experiment_outcome_requires_experiment_observer(self):
        bad = copy.deepcopy(self.events)
        bad.append({
            "event_id":"EV-X",
            "generation":2,
            "event_type":"experiment_outcome",
            "subject":"E1_ARID5B_OCCUPANCY",
            "observer":"model",
            "evidence_refs":["test:measurement"],
            "payload":{"experiment_id":"E1_ARID5B_OCCUPANCY","outcome_name":"no_allele_specific_occupancy"}
        })
        with self.assertRaises(CausalProcessorError):
            CausalProcessor(bad)

    def test_experiment_outcome_requires_provenance(self):
        bad = copy.deepcopy(self.events)
        bad.append({
            "event_id":"EV-X",
            "generation":2,
            "event_type":"experiment_outcome",
            "subject":"E1_ARID5B_OCCUPANCY",
            "observer":"experiment",
            "evidence_refs":[],
            "payload":{"experiment_id":"E1_ARID5B_OCCUPANCY","outcome_name":"no_allele_specific_occupancy"}
        })
        with self.assertRaises(CausalProcessorError):
            CausalProcessor(bad)

    def test_next_generation_does_not_rewrite_history(self):
        p = CausalProcessor(self.events)
        old = copy.deepcopy(p.events)
        gen = p.next_generation([
            {"event_id":"EV-0009","generation":2,"event_type":"hypothesis_superseded","subject":"H1_ARID5B_REUSE","observer":"model","evidence_refs":["future:test"],"payload":{}},
            {"event_id":"EV-0010","generation":2,"event_type":"model_update","subject":"GAP-001","observer":"model","evidence_refs":["future:test"],"payload":{"cause_found":False}}
        ])
        self.assertEqual(2, gen)
        self.assertEqual(old, p.events[:len(old)])
        self.assertIn("H1_ARID5B_REUSE", p.project().superseded_hypotheses)


if __name__ == "__main__":
    unittest.main()
