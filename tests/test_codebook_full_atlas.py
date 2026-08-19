import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


FETCH = load("gap001_fetch_codebook_atlas", ROOT / "analysis" / "gap001_fetch_codebook_atlas.py")
PREP = load("gap001_prepare_codebook_atlas", ROOT / "analysis" / "gap001_prepare_codebook_atlas.py")


class CodebookAtlasFetchTests(unittest.TestCase):
    def test_resolves_explicit_supplementary_data_1_anchor(self):
        html = '''<html><body>
        <a href="/other.zip">Supplementary Tables</a>
        <a href="https://media.springernature.com/codebook-data1.zip">Supplementary Data 1 (download ZIP)</a>
        </body></html>'''
        got = FETCH.resolve_supplementary_data_1("https://www.nature.com/articles/x", html)
        self.assertEqual(got, "https://media.springernature.com/codebook-data1.zip")

    def test_does_not_guess_when_label_is_missing(self):
        html = '<a href="x.zip">Supplementary Data 2</a>'
        with self.assertRaises(RuntimeError):
            FETCH.resolve_supplementary_data_1("https://www.nature.com/articles/x", html)

    def test_zip_validation_rejects_non_zip(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "x.zip"
            p.write_text("not a zip")
            with self.assertRaises(RuntimeError):
                FETCH.validate_zip(p)

    def test_springer_prefix_is_derived_from_exact_doi(self):
        base, stem = FETCH.springer_media_prefix("10.1038/s41586-026-10798-9")
        self.assertIn("art%3A10.1038%2Fs41586-026-10798-9", base)
        self.assertEqual(stem, "41586_2026_10798")

    def test_motif_likeness_requires_motif_signal(self):
        generic = [{"name": "table1.xlsx", "size": 10}, {"name": "source.csv", "size": 10}]
        self.assertFalse(FETCH.motif_likeness(generic)[2])
        motifs = [{"name": f"TF{i}.pwm", "size": 10} for i in range(100)]
        self.assertTrue(FETCH.motif_likeness(motifs)[2])

    def test_select_motif_zip_requires_unique_candidate(self):
        candidates = [
            {"looks_motif_like": False, "index": 1, "url": "a", "final_url": "a", "sha256": "a", "size_bytes": 1, "member_count": 1, "motif_named_members": 0, "motif_extension_members": 0},
            {"looks_motif_like": True, "index": 2, "url": "b", "final_url": "b", "sha256": "b", "size_bytes": 2, "member_count": 2, "motif_named_members": 2, "motif_extension_members": 2, "data": b"x", "members": []},
        ]
        self.assertEqual(FETCH.select_motif_zip(candidates)["index"], 2)
        candidates.append({"looks_motif_like": True, "index": 3, "url": "c", "final_url": "c", "sha256": "c", "size_bytes": 3, "member_count": 3, "motif_named_members": 3, "motif_extension_members": 3})
        with self.assertRaises(RuntimeError):
            FETCH.select_motif_zip(candidates)


class CodebookAtlasParserTests(unittest.TestCase):
    def test_official_codebook_tsv(self):
        text = (
            b"TF\tARID5B\n"
            b"Motif\tM01234_1.94d\n"
            b"Pos\tA\tC\tG\tT\n"
            b"1\t0.1\t0.2\t0.3\t0.4\n"
            b"2\t0.4\t0.3\t0.2\t0.1\n"
        )
        motifs = PREP.parse_member("SupplementaryData1/ARID5B.txt", text)
        self.assertEqual(motifs, [("ARID5B", [[0.1, 0.2, 0.3, 0.4], [0.4, 0.3, 0.2, 0.1]])])

    def test_codebook_tsv_requires_monotonic_positions(self):
        text = (
            "TF\tTFX\nMotif\tM1\nPos\tA\tC\tG\tT\n"
            "1\t0.1\t0.2\t0.3\t0.4\n3\t0.4\t0.3\t0.2\t0.1\n"
        )
        self.assertEqual(PREP.parse_codebook_tsv(text), [])

    def test_simple_lx4(self):
        motifs = PREP.parse_member("TFX.pwm", b">TFX\n1 2 3 4\n4 3 2 1\n")
        self.assertEqual(motifs[0][0], "TFX")
        self.assertEqual(len(motifs[0][1]), 2)

    def test_jaspar_rows(self):
        text = b">MA0001.1 TFX\nA [1 2 3]\nC [2 3 4]\nG [3 4 5]\nT [4 5 6]\n"
        motifs = PREP.parse_member("x.pfm", text)
        self.assertEqual(motifs[0][0], "TFX")
        self.assertEqual(motifs[0][1][0], [1.0, 2.0, 3.0, 4.0])

    def test_meme(self):
        text = b'''MEME version 4\n\nALPHABET= ACGT\nMOTIF X TFX\nletter-probability matrix: alength= 4 w= 2\n0.1 0.2 0.3 0.4\n0.4 0.3 0.2 0.1\n'''
        motifs = PREP.parse_member("motifs.meme", text)
        self.assertEqual(motifs[0][0], "TFX")
        self.assertEqual(len(motifs[0][1]), 2)

    def test_transfac(self):
        text = b'''ID TFX\nP0 A C G T\n01 1 2 3 4\n02 4 3 2 1\n//\n'''
        motifs = PREP.parse_member("motifs.transfac", text)
        self.assertEqual(motifs[0][0], "TFX")
        self.assertEqual(motifs[0][1][1], [4.0, 3.0, 2.0, 1.0])

    def test_duplicate_matrices_can_remain_distinct_representatives(self):
        a = [[1.0, 2.0, 3.0, 4.0], [4.0, 3.0, 2.0, 1.0]]
        b = [[1.0, 2.0, 3.0, 4.0], [4.0, 3.0, 2.0, 2.0]]
        self.assertNotEqual(PREP.matrix_signature(a), PREP.matrix_signature(b))
        # Different representatives for one TF are valid published records;
        # matrix inequality itself must not be treated as an atlas conflict.


if __name__ == "__main__":
    unittest.main()
