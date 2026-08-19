from pathlib import Path
import importlib.util
import sys

MODULE = Path(__file__).resolve().parents[1] / "analysis" / "gap001_mecp2_chip_locus.py"
spec = importlib.util.spec_from_file_location("gap001_mecp2_chip_locus", MODULE)
mod = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = mod
spec.loader.exec_module(mod)


def write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return p


def test_fixedstep_parser(tmp_path):
    p = write(tmp_path, "x.wig", "fixedStep chrom=chr8 start=101 step=10 span=5\n1\n2\n")
    rows = list(mod.iter_wig(p))
    assert [(r.start0, r.end0, r.value) for r in rows] == [(100, 105, 1.0), (110, 115, 2.0)]


def test_variablestep_parser(tmp_path):
    p = write(tmp_path, "x.wig", "variableStep chrom=chr8 span=3\n101 2.5\n")
    rows = list(mod.iter_wig(p))
    assert rows[0].start0 == 100 and rows[0].end0 == 103 and rows[0].value == 2.5


def test_bedgraph_parser(tmp_path):
    p = write(tmp_path, "x.wig", "chr8\t100\t110\t3.0\n")
    rows = list(mod.iter_wig(p))
    assert rows[0] == mod.WigPoint("chr8", 100, 110, 3.0)


def test_summarize_ratio():
    chip = [mod.WigPoint("chr8", 90, 111, 4.0)]
    inp = [mod.WigPoint("chr8", 90, 111, 2.0)]
    out = mod.summarize(chip, inp, 100)
    assert out["exact_base_signal"]["chip_over_input"] == 2.0
    assert out["windows"]["250"]["chip_over_input"] == 2.0


def test_no_causal_promotion_fields_in_summary():
    out = mod.summarize([mod.WigPoint("chr8", 0, 10, 1.0)], [mod.WigPoint("chr8", 0, 10, 1.0)], 5)
    assert "causal_status" not in out
    assert "materialized" not in out
