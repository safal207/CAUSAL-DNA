import importlib.util
import io
from pathlib import Path
import sys
import tempfile
import unittest
import zipfile

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


class CodebookAtlasParserTests(unittest.TestCase):
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

    def test_duplicate_conflict_signature_differs(self):
        a = [[1.0, 2.0, 3.0, 4.0], [4.0, 3.0, 2.0, 1.0]]
        b = [[1.0, 2.0, 3.0, 4.0], [4.0, 3.0, 2.0, 2.0]]
        self.assertNotEqual(PREP.matrix_signature(a), PREP.matrix_signature(b))


if __name__ == "__main__":
    unittest.main()
