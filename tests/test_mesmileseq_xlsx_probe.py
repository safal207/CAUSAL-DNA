import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest
import zipfile

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "gap001_probe_mesmileseq_xlsx",
    ROOT / "analysis" / "gap001_probe_mesmileseq_xlsx.py",
)
MOD = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MOD
assert SPEC.loader is not None
SPEC.loader.exec_module(MOD)


class MeSmileXlsxProbeTests(unittest.TestCase):
    def build_xlsx(self, path: Path) -> None:
        workbook = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
        <workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"
          xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">
          <sheets><sheet name="TF classes" sheetId="1" r:id="rId1"/></sheets>
        </workbook>'''
        rels = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
        <Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
          <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/>
        </Relationships>'''
        shared = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
        <sst xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" count="3" uniqueCount="3">
          <si><t>TF</t></si><si><t>class</t></si><si><t>CXXC4</t></si>
        </sst>'''
        sheet = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
        <worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
          <sheetData>
            <row r="1"><c r="A1" t="s"><v>0</v></c><c r="B1" t="s"><v>1</v></c></row>
            <row r="2"><c r="A2" t="s"><v>2</v></c><c r="B2" t="inlineStr"><is><t>methyl minus</t></is></c><c r="C2"><v>1.5</v></c></row>
          </sheetData>
        </worksheet>'''
        with zipfile.ZipFile(path, "w") as zf:
            zf.writestr("xl/workbook.xml", workbook)
            zf.writestr("xl/_rels/workbook.xml.rels", rels)
            zf.writestr("xl/sharedStrings.xml", shared)
            zf.writestr("xl/worksheets/sheet1.xml", sheet)

    def test_column_index(self):
        self.assertEqual(MOD.column_index("A1"), 0)
        self.assertEqual(MOD.column_index("Z1"), 25)
        self.assertEqual(MOD.column_index("AA9"), 26)

    def test_inspect_reads_shared_inline_and_numeric_cells(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "supp.xlsx"
            self.build_xlsx(path)
            out = MOD.inspect_xlsx(path, max_rows=10)
            self.assertEqual(out["sheet_count"], 1)
            self.assertEqual(out["sheets"][0]["title"], "TF classes")
            self.assertEqual(out["sheets"][0]["sample_rows"][0], ["TF", "class"])
            self.assertEqual(out["sheets"][0]["sample_rows"][1], ["CXXC4", "methyl minus", 1.5])

    def test_non_xlsx_is_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "bad.xlsx"
            path.write_text("not xlsx")
            with self.assertRaises(ValueError):
                MOD.inspect_xlsx(path)


if __name__ == "__main__":
    unittest.main()
