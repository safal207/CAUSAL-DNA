import copy
import json
import unittest
from pathlib import Path

from causal_dna.orientation_kernel import OrientationKernelError, OrientationKernelGraph

ROOT = Path(__file__).resolve().parents[1]
CASE_PATH = ROOT / "cases" / "CDNA-001.exchange-process-graph.json"


def node(node_id):
    return {"id": node_id, "class": "CHROMATIN_STATE", "label": node_id}


def edge(start, end, status="INTERVENTION_SUPPORTED", relation="CAUSES"):
    return {
        "from": start,
        "to": end,
        "relation": relation,
        "status": status,
        "evidence_refs": [f"EXP:{start}->{end}"] if status != "HYPOTHESIS" else [],
    }


def synthetic(nodes, edges, source="S", target="T", declared_kernel=None, kernel_found=False, gap_status="CLOSED"):
    return {
        "graph_id": "SYNTH",
        "case_id": "TEST",
        "gap_id": "GAP",
        "gap_status": gap_status,
        "source": source,
        "target": target,
        "kernel_found": kernel_found,
        "declared_kernel": declared_kernel or [],
        "nodes": [node(n) for n in nodes],
        "edges": edges,
    }


class OrientationKernelTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reference = json.loads(CASE_PATH.read_text(encoding="utf-8"))

    def test_reference_graph_is_valid_but_has_no_kernel(self):
        graph = OrientationKernelGraph(copy.deepcopy(self.reference))
        graph.assert_valid()
        self.assertEqual([], graph.minimal_kernels())
        self.assertIsNone(graph.certificate())
        self.assertFalse(graph.document["kernel_found"])

    def test_observational_edges_do_not_create_causal_path(self):
        doc = synthetic(
            ["S", "M", "T"],
            [edge("S", "M", status="OBSERVATIONAL"), edge("M", "T", status="OBSERVATIONAL")],
        )
        graph = OrientationKernelGraph(doc)
        self.assertEqual([], graph.causal_paths())
        self.assertEqual([], graph.minimal_kernels())

    def test_single_intermediate_node_is_minimal_kernel(self):
        doc = synthetic(
            ["S", "M", "T"],
            [edge("S", "M"), edge("M", "T")],
            declared_kernel=["M"],
            kernel_found=True,
        )
        graph = OrientationKernelGraph(doc)
        graph.assert_valid()
        self.assertEqual([("M",)], graph.minimal_kernels())
        cert = graph.certificate()
        self.assertEqual(("M",), cert.kernel)
        self.assertEqual("INTERVENTION_SUPPORTED_KERNEL", cert.status)

    def test_parallel_paths_require_combined_kernel(self):
        doc = synthetic(
            ["S", "M1", "M2", "T"],
            [
                edge("S", "M1"), edge("M1", "T"),
                edge("S", "M2"), edge("M2", "T"),
            ],
            declared_kernel=["M1", "M2"],
            kernel_found=True,
        )
        graph = OrientationKernelGraph(doc)
        graph.assert_valid()
        self.assertEqual([("M1", "M2")], graph.minimal_kernels())

    def test_direct_causal_bypass_prevents_intermediate_kernel(self):
        doc = synthetic(
            ["S", "M", "T"],
            [edge("S", "M"), edge("M", "T"), edge("S", "T")],
        )
        graph = OrientationKernelGraph(doc)
        self.assertEqual([], graph.minimal_kernels())

    def test_feedback_edge_is_not_used_as_forward_causal_path(self):
        doc = synthetic(
            ["S", "M", "T"],
            [edge("S", "M"), edge("M", "T"), edge("T", "M", relation="FEEDBACK")],
            declared_kernel=["M"],
            kernel_found=True,
        )
        graph = OrientationKernelGraph(doc)
        graph.assert_valid()
        self.assertEqual([("S", "M", "T")], graph.causal_paths())

    def test_all_verified_path_edges_upgrade_certificate(self):
        doc = synthetic(
            ["S", "M", "T"],
            [edge("S", "M", status="VERIFIED"), edge("M", "T", status="VERIFIED")],
            declared_kernel=["M"],
            kernel_found=True,
        )
        graph = OrientationKernelGraph(doc)
        graph.assert_valid()
        cert = graph.certificate()
        self.assertTrue(cert.independently_verified)
        self.assertEqual("VERIFIED_KERNEL", cert.status)

    def test_declared_kernel_cannot_override_evidence(self):
        doc = synthetic(
            ["S", "M", "T"],
            [edge("S", "M", status="OBSERVATIONAL"), edge("M", "T", status="OBSERVATIONAL")],
            declared_kernel=["M"],
            kernel_found=True,
        )
        graph = OrientationKernelGraph(doc)
        with self.assertRaises(OrientationKernelError):
            graph.assert_valid()

    def test_open_gap_cannot_claim_kernel_found(self):
        doc = synthetic(
            ["S", "M", "T"],
            [edge("S", "M"), edge("M", "T")],
            declared_kernel=["M"],
            kernel_found=True,
            gap_status="OPEN",
        )
        graph = OrientationKernelGraph(doc)
        self.assertTrue(any("gap_status=OPEN forbids kernel_found=true" in e for e in graph.errors()))


if __name__ == "__main__":
    unittest.main()
