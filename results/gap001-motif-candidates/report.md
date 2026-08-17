# GAP-001 motif-candidate presence screen

Status: **PASS — LHX6/NKX6-3 direct adult-VPH routes demoted; ARID5B/CUX1 remain cell-feasible**

GitHub Actions run: `31988331718`
Artifact: `gap001-motif-candidates`
Artifact digest: `sha256:9c4d4e547333c2d5b350c33d27700d5a0b887c677a675d19619271a64eb07f33`

## Why this screen

The original rs1421085 motif analysis reported that the T>C substitution disrupts predicted motifs for ARID-family regulators, LHX6 and NKX6-3. ARID5B was prioritized experimentally in adipocytes because it was the best adipocyte-context regulator, not because the other motif candidates were universal impossibilities.

LHX6 is independently known to mark neurons in the posterior hypothalamic region, making it a biologically relevant brain-context candidate worth testing.

This analysis asks only:

> Are these candidate TF transcripts detected in the same wild-type VPH cells as `Irx3`?

## Pooled VPH cells

Across 44,149 GSE146692 filtered cells:

| candidate | candidate+ cells | Irx3+/candidate+ | candidate+ among Irx3+ |
|---|---:|---:|---:|
| Arid5b | 7,132 | 509 | **28.420%** |
| Cux1 | 9,279 | 632 | **35.288%** |
| Lhx6 | 282 | 11 | **0.614%** |
| Nkx6-3 | 200 | 5 | **0.279%** |

All four target genes were present as unique features in all four Cell Ranger matrices.

## By sex (descriptive only)

| sex | candidate | candidate+ among Irx3+ |
|---|---|---:|
| female | Arid5b | 28.015% |
| male | Arid5b | 28.763% |
| female | Cux1 | 35.932% |
| male | Cux1 | 34.742% |
| female | Lhx6 | 0.853% |
| male | Lhx6 | 0.412% |
| female | Nkx6-3 | 0.365% |
| male | Nkx6-3 | 0.206% |

Sex comparisons are not biological-replicate inference because sex is confounded with harvest/10x chemistry in GSE146692.

## Interpretation

### LHX6

**Status: `DIRECT_ADULT_CELL_MECHANISM_DEMOTED; DEVELOPMENTAL_ROUTE_OPEN`.**

LHX6 is attractive on sequence/anatomical grounds: the rs1421085 substitution was reported to disrupt an LHX6 motif, and Lhx6+ neurons exist in posterior hypothalamic regions. However, only 11 of 1,791 `Irx3+` cells co-detected `Lhx6` transcript in this VPH dataset.

That makes a broad direct model such as:

```text
adult PH LHX6 binds rs1421085 element
  -> directly regulates Irx3 in most Irx3+ cells
```

poorly supported by current single-cell feasibility evidence.

It does **not** fully exclude LHX6 because:

- scRNA-seq dropout can miss low TF transcripts;
- protein can persist after transcript falls;
- LHX6 could act in a rare causal subpopulation;
- a developmental/imprinting mechanism could act before the sampled P30–P35 state.

### NKX6-3

**Status: `DIRECT_ADULT_CELL_MECHANISM_STRONGLY_DEMOTED`.**

Only 5 `Irx3+` cells co-detected `Nkx6-3`. There is currently even less cellular support for NKX6-3 as the adult PH rs1421085→Irx3 interpreter.

### ARID5B and CUX1

Both remain substantially more cell-feasible than LHX6/NKX6-3 by transcript co-detection:

```text
P(Arid5b+ | Irx3+) = 28.4%
P(Cux1+   | Irx3+) = 35.3%
```

But neither has been shown to occupy the rs1421085 element in the relevant adult PH cell type, and the CUX1 neuronal literature currently points to an FTO/RPGRIP1L target program rather than a general IRX3 program.

## Updated ranking for the missing adult PH edge

1. **H3/H4 — allele-dependent chromatin state/contact**: rises, because universal single-TF explanations keep failing across contexts.
2. **H1 — ARID5B reuse**: cell-feasible but under negative tension from DEL82 bulk hypothalamus.
3. **H2a — CUX1 combinatorial/parallel role**: cell-feasible and exact-site-sensitive, but target-gene conflict remains.
4. **H9 — temporal/developmental gating**: remains strong architectural explanation.
5. **LHX6 direct adult route**: demoted.
6. **NKX6-3 direct adult route**: strongly demoted.

## Next discriminator

The highest-value test is shifting from `which TF is present?` to:

```text
exact T vs C allele
  x male vs female
  x Irx3+ neuronal/astrocyte subtype
    -> accessibility
    -> TF occupancy
    -> enhancer-promoter contact
    -> Irx3 transcription
```

This screen reduces candidate space; it does not close GAP-001.
