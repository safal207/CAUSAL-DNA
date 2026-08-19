import importlib.util
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "gap001_gse226277_multiome_manifest",
    ROOT / "analysis" / "gap001_gse226277_multiome_manifest.py",
)
MOD = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MOD
SPEC.loader.exec_module(MOD)


def sample_block(acc, title, modality, condition, replicate):
    suffix = "matrix.h5" if modality == "RNA" else "atac_fragments.tsv.gz"
    return f"""^SAMPLE = {acc}
!Sample_title = {title}
!Sample_characteristics_ch1 = strain: C57BL6N
!Sample_characteristics_ch1 = genotype: WT
!Sample_characteristics_ch1 = phenotype: diabetic (HFD+STZ)
!Sample_characteristics_ch1 = Sex: male
!Sample_characteristics_ch1 = tissue: Hypothalamus
!Sample_supplementary_file_1 = ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSMx/{acc}/suppl/{acc}_{suffix}
"""


class MultiomeManifestTests(unittest.TestCase):
    def fixture(self):
        rows = []
        n = 100
        for condition in ("RH", "AH"):
            for replicate in (1, 2):
                for modality in ("RNA", "ATAC"):
                    acc = f"GSM{n}"
                    n += 1
                    title = f"Hypothalamus {condition}, replicate {replicate}, Multiome ({modality})"
                    rows.append(sample_block(acc, title, modality, condition, replicate))
        return "\n".join(rows)

    def test_exact_four_pairs_are_built(self):
        selected = MOD.select_hypothalamus_multiome(MOD.parse_soft(self.fixture()))
        pairs = MOD.pair_samples(selected)
        self.assertEqual(len(pairs), 4)
        self.assertEqual({p["pair_id"] for p in pairs}, {"HYPO-AH-R1", "HYPO-AH-R2", "HYPO-RH-R1", "HYPO-RH-R2"})

    def test_each_pair_has_rna_and_atac(self):
        selected = MOD.select_hypothalamus_multiome(MOD.parse_soft(self.fixture()))
        for pair in MOD.pair_samples(selected):
            self.assertTrue(pair["rna_accession"].startswith("GSM"))
            self.assertTrue(pair["atac_accession"].startswith("GSM"))
            self.assertEqual(pair["sex"].lower(), "male")
            self.assertEqual(pair["genotype"].upper(), "WT")

    def test_ftp_urls_are_normalized_to_https(self):
        samples = MOD.parse_soft(self.fixture())
        urls = [u for s in samples for u in s["supplementary_files"]]
        self.assertTrue(urls)
        self.assertTrue(all(u.startswith("https://ftp.ncbi.nlm.nih.gov/") for u in urls))

    def test_missing_modality_is_rejected(self):
        selected = MOD.select_hypothalamus_multiome(MOD.parse_soft(self.fixture()))
        selected = [r for r in selected if not (r["condition"] == "RH" and r["replicate"] == 1 and r["modality"] == "ATAC")]
        with self.assertRaises(ValueError):
            MOD.pair_samples(selected)

    def test_wrong_sex_is_rejected(self):
        text = self.fixture().replace("!Sample_characteristics_ch1 = Sex: male", "!Sample_characteristics_ch1 = Sex: female", 1)
        with self.assertRaises(ValueError):
            MOD.select_hypothalamus_multiome(MOD.parse_soft(text))


if __name__ == "__main__":
    unittest.main()
