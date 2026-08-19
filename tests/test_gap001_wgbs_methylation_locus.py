from pathlib import Path
import importlib.util
import sys

MODULE = Path(__file__).resolve().parents[1] / "analysis" / "gap001_wgbs_methylation_locus.py"
spec = importlib.util.spec_from_file_location("gap001_wgbs_methylation_locus", MODULE)
mod = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = mod
spec.loader.exec_module(mod)


def test_forward_contexts():
    seq="AACGCACAAA"
    assert mod.classify_context(seq,100,102)=="CG"
    assert mod.classify_context(seq,100,104)=="CAC"


def test_reverse_contexts():
    # Reverse complement of CG is CG: terminal G at index 3 corresponds to reverse-strand CG.
    seq="AACGAAA"
    assert mod.classify_context(seq,100,103)=="CG"
    # TGT reverse-complements to ACA, so use GTG: reverse complement CAC.
    seq2="AAGTGAA"
    assert mod.classify_context(seq2,100,104)=="CAC"


def test_scan_cov_and_min_coverage(tmp_path):
    fasta=tmp_path/'r.fa'; fasta.write_text('>chr8\nAACGCACAAA\n')
    cov=tmp_path/'x.cov'; cov.write_text('chr8 103 103 80 8 2\nchr8 105 105 50 2 2\n')
    _,seq=mod.read_fasta(fasta)
    calls=mod.scan_cov(cov,'chr8',102,10,seq,100,5)
    assert len(calls)==1
    assert calls[0].context=='CG'
    assert calls[0].fraction==0.8


def test_aggregate_is_coverage_weighted():
    calls=[mod.Call(99,9,1,'CG'),mod.Call(101,1,9,'CG')]
    out=mod.aggregate(calls,100,10,'CG')
    assert out['sites']==2
    assert out['coverage_weighted_methylation']==0.5


def test_no_causal_promotion_in_aggregate():
    out=mod.aggregate([mod.Call(1,1,1,'CG')],1,10)
    assert 'cause_found' not in out
    assert 'causal_status' not in out
