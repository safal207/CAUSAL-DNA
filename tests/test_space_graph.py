import copy
import json
import unittest
from pathlib import Path

from causal_dna.space_graph import SpaceGraph, SpaceGraphError

ROOT = Path(__file__).resolve().parents[1]
GRAPH_PATH = ROOT / "cases" / "CDNA-001-rs1421085.space-graph.json"


class SpaceGraphTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.document = json.loads(GRAPH_PATH.read_text(encoding="utf-8"))

    def graph(self):
        return SpaceGraph(copy.deepcopy(self.document))

    def test_reference_graph_is_valid(self):
        graph = self.graph()
        graph.assert_valid()
        self.assertEqual([], graph.errors())

    def test_all_three_spaces_are_populated(self):
        graph = self.graph()
        self.assertGreater(len(graph.nodes_in("projective")), 0)
        self.assertGreater(len(graph.nodes_in("bardo")), 0)
        self.assertGreater(len(graph.nodes_in("material")), 0)

    def test_candidate_hypothesis_has_full_contour(self):
        graph = self.graph()
        contours = graph.contours("P_H3_ACCESS", "M_IRX3_UP")
        self.assertGreaterEqual(len(contours), 1)
        self.assertTrue(
            all(
                path.compressed_spaces == ("projective", "bardo", "material")
                for path in contours
            )
        )

    def test_material_nodes_require_evidence(self):
        doc = copy.deepcopy(self.document)
        node = next(n for n in doc["nodes"] if n["id"] == "M_IRX3_UP")
        node["evidence_refs"] = []
        graph = SpaceGraph(doc)
        with self.assertRaises(SpaceGraphError):
            graph.assert_valid()

    def test_projective_cannot_bypass_bardo(self):
        doc = copy.deepcopy(self.document)
        doc["edges"].append(
            {
                "id": "SG_TEST_BYPASS",
                "from": "P_H3_ACCESS",
                "to": "M_IRX3_UP",
                "relation": "supports",
                "status": "supported",
                "evidence_refs": ["fake:evidence"],
            }
        )
        graph = SpaceGraph(doc)
        errors = graph.errors()
        self.assertTrue(any("projective→material bypass forbidden" in e for e in errors))

    def test_edge_entering_material_requires_provenance(self):
        doc = copy.deepcopy(self.document)
        edge = next(e for e in doc["edges"] if e["id"] == "SG_E040")
        edge["evidence_refs"] = []
        graph = SpaceGraph(doc)
        errors = graph.errors()
        self.assertTrue(any("entering material space requires evidence_refs" in e for e in errors))

    def test_materialization_relation_only_bardo_to_material(self):
        doc = copy.deepcopy(self.document)
        edge = next(e for e in doc["edges"] if e["id"] == "SG_E018")
        edge["relation"] = "materializes_as"
        graph = SpaceGraph(doc)
        errors = graph.errors()
        self.assertTrue(any("materializes_as must connect Bardo→material" in e for e in errors))

    def test_open_bardo_is_preserved(self):
        graph = self.graph()
        ids = {node["id"] for node in graph.open_bardo()}
        self.assertIn("B_H1_ARID5B", ids)
        self.assertIn("B_H4_CONTACT", ids)
        self.assertIn("B_INTEGRATION_GATE", ids)

    def test_materialization_frontier_contains_decisive_gate(self):
        graph = self.graph()
        ids = {node["id"] for node in graph.materialization_frontier()}
        self.assertIn("B_INTEGRATION_GATE", ids)

    def test_material_backbone_remains_queryable(self):
        graph = self.graph()
        paths = list(graph.paths("M_RS1421085_EDIT", "M_WEIGHT_UP"))
        self.assertGreaterEqual(len(paths), 1)
        self.assertTrue(any("M_IRX3_UP" in p.node_ids for p in paths))

    def test_mermaid_is_generated_from_same_graph(self):
        rendered = self.graph().to_mermaid()
        self.assertIn("flowchart LR", rendered)
        self.assertIn('subgraph PROJECTIVE["Projective space"]', rendered)
        self.assertIn('subgraph BARDO["Bardo transition space"]', rendered)
        self.assertIn('subgraph MATERIAL["Material space"]', rendered)
        self.assertIn("P_H3_ACCESS", rendered)
        self.assertIn("B_INTEGRATION_GATE", rendered)
        self.assertIn("M_IRX3_UP", rendered)
        self.assertIn('SG_E040', json.dumps(self.document))


if __name__ == "__main__":
    unittest.main()
