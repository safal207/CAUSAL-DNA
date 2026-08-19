import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest
import zipfile

ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "analysis"
if str(ANALYSIS) not in sys.path:
    sys.path.insert(0, str(ANALYSIS))
SPEC = importlib.util.spec_from_file_location(
    "gap001_mesmileseq_methylation_triage",
    ANALYSIS / "gap001_mesmileseq_methylation_triage.py",
)
MOD = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MOD
assert SPEC.loader is not None
SPEC.loader.exec_module(MOD)


class MeSmileMethylationTriageTests(unittest.TestCase):
    def build_xlsx(self, path: Path) -> None:
        strings = [
            "TF", "DBD", "classification", "misc",
            "CXXC4", "CXXC", "methyl minus", "",
            "ZNF445", "C2H2 ZF", "weak methyl plus", "note",
        ]
        shared_items = "".join(f"<si><t>{s}</t></si>" for s in strings)
        workbook = '''<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets><sheet name="Classification by methylation" sheetId="1" r:id="rId1"/></sheets></workbook>'''
        rels = '''<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="worksheet" Target="worksheets/sheet1.xml"/></Relationships>'''
        shared = f'''<sst xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">{shared_items}</sst>'''
        # Header 0..3; CXXC4 4..7; ZNF445 8..11.
        sheet = '''<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheetData>
        <row r="1"><c r="A1" t="s"><v>0</v></c><c r="B1" t="s"><v>1</v></c><c r="C1" t="s"><v>2</v></c><c r="D1" t="s"><v>3</v></c></row>
        <row r="2"><c r="A2" t="s"><v>4</v></c><c r="B2" t="s"><v>5</v></c><c r="C2" t="s"><v>6</v></c><c r="D2" t="s"><v>7</v></c></row>
        <row r="3"><c r="A3" t="s"><v>8</v></c><c r="B3" t="s"><v>9</v></c><c r="C3" t="s"><v>10</v></c><c r="D3" t="s"><v>11</v></c></row>
        </sheetData></worksheet>'''
        with zipfile.ZipFile(path, "w") as zf:
            zf.writestr("xl/workbook.xml", workbook)
            zf.writestr("xl/_rels/workbook.xml.rels", rels)
            zf.writestr("xl/sharedStrings.xml", shared)
            zf.writestr("xl/worksheets/sheet1.xml", sheet)

    def test_normalize_preserves_weakness(self):
        self.assertEqual(MOD.normalize_class("methyl plus"), ("methyl_plus", "strong"))
        self.assertEqual(MOD.normalize_class("weak methyl minus"), ("methyl_minus", "weak"))
        self.assertEqual(MOD.normalize_class("no CG"), ("little_effect_no_cg", "no_cg"))
        self.assertEqual(MOD.normalize_class("*"), ("unknown", "*"))

    def test_join_unknown_is_not_negative(self):
        classes = [{
            "tf": "CXXC4", "dbd": "CXXC", "published_classification": "methyl minus",
            "mesmile_class": "methyl_minus", "class_strength": "strong", "misc": "",
            "mesmile_profiled": True, "source_sheet": MOD.CLASS_SHEET,
        }]
        motifs = [
            {"tf": "CXXC4", "threshold_075_crossing": "false", "max_site_strength": "0.62", "abs_delta_normalized": "0.34"},
            {"tf": "NEWTF", "threshold_075_crossing": "true", "max_site_strength": "0.80", "abs_delta_normalized": "0.20"},
        ]
        joined = MOD.join_triage(motifs, classes, [])
        by_tf = {x["tf"]: x for x in joined}
        self.assertEqual(by_tf["NEWTF"]["mesmile_class"], "unknown")
        self.assertFalse(by_tf["NEWTF"]["mesmile_profiled"])
        self.assertEqual(by_tf["CXXC4"]["mesmile_class"], "methyl_minus")
        self.assertTrue(by_tf["CXXC4"]["methylation_relevant_class"])

    def test_threshold_crossing_has_priority_without_composite_score(self):
        classes = []
        motifs = [
            {"tf": "A", "threshold_075_crossing": "false", "max_site_strength": "0.95", "abs_delta_normalized": "0.90"},
            {"tf": "B", "threshold_075_crossing": "true", "max_site_strength": "0.76", "abs_delta_normalized": "0.10"},
        ]
        joined = MOD.join_triage(motifs, classes, [])
        self.assertEqual(joined[0]["tf"], "B")
        self.assertNotIn("composite_score", joined[0])

    def test_published_sheet_minimum_guard(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "small.xlsx"
            self.build_xlsx(path)
            with self.assertRaises(ValueError):
                MOD.published_classification_rows(path)


if __name__ == "__main__":
    unittest.main()
