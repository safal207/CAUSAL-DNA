import importlib.util
from pathlib import Path
import gzip
import sys
import tempfile
import unittest

try:
    import numpy as np
except ModuleNotFoundError as exc:
    raise unittest.SkipTest("scientific stack is validated in the dedicated GSE226277 workflow") from exc

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "gap001_gse226277_same_nucleus_locus",
    ROOT / "analysis" / "gap001_gse226277_same_nucleus_locus.py",
)
MOD = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MOD
SPEC.loader.exec_module(MOD)


class SameNucleusLocusTests(unittest.TestCase):
    def test_stream_extracts_only_overlapping_locus_barcodes(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "fragments.tsv.gz"
            rows = [
                "chr7\t10\t30\tX-1\t1\n",
                "chr8\t91373000\t91373100\tA-1\t1\n",
                "chr8\t91373990\t91374010\tA-1\t2\n",
                "chr8\t91374100\t91374200\tB-1\t1\n",
                "chr8\t91374490\t91374510\tB-1\t3\n",
                "chr8\t91374500\t91374600\tC-1\t1\n",
                "chr9\t1\t10\tD-1\t1\n",
            ]
            with gzip.open(path, "wt") as handle:
                handle.writelines(rows)
            locus, meta = MOD.stream_locus_fragments(path.as_uri(), "chr8", 91374000, 91374500)
            self.assertEqual(set(locus), {"A-1", "B-1"})
            self.assertEqual(locus["B-1"]["read_support"], 4)
            self.assertTrue(meta["target_chromosome_seen"])
            self.assertTrue(meta["early_stop_after_locus"])

    def test_join_requires_barcode_intersection(self):
        barcodes = np.asarray(["A-1", "B-1", "C-1"])
        expr = {label: np.zeros(3, dtype=int) for label in MOD.GENES}
        expr["IRX3"] = np.asarray([1, 0, 1])
        expr["ARID5B"] = np.asarray([2, 1, 0])
        locus = {
            "A-1": {"fragment_records": 1, "read_support": 2},
            "B-1": {"fragment_records": 2, "read_support": 2},
            "OUTSIDE-1": {"fragment_records": 1, "read_support": 1},
        }
        joined, summary = MOD.build_join(barcodes, expr, locus)
        self.assertEqual({r["barcode"] for r in joined}, {"A-1", "B-1"})
        self.assertEqual(summary["raw_locus_fragment_barcodes"], 3)
        self.assertEqual(summary["locus_barcodes_in_filtered_rna"], 2)
        self.assertEqual(summary["locus_irx3_positive"], 1)
        self.assertEqual(summary["candidate_same_nucleus"]["ARID5B"]["irx3_candidate_positive_among_locus_barcodes"], 1)

    def test_zero_irx3_locus_is_preserved_not_promoted(self):
        barcodes = np.asarray(["A-1"])
        expr = {label: np.zeros(1, dtype=int) for label in MOD.GENES}
        locus = {"A-1": {"fragment_records": 1, "read_support": 1}}
        _, summary = MOD.build_join(barcodes, expr, locus)
        self.assertEqual(summary["locus_irx3_positive"], 0)
        self.assertFalse(MOD.FORBIDDEN_CAUSAL_FIELDS & set(summary))

    def test_invalid_interval_is_rejected(self):
        with self.assertRaises(ValueError):
            MOD.stream_locus_fragments("file:///does/not/matter", "chr8", 10, 10)


if __name__ == "__main__":
    unittest.main()
