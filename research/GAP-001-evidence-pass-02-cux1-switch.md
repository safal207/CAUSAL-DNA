# GAP-001 Evidence Pass 02 — exact-site CUX1 switch candidate

Date: 2026-08-17
Status: competing mechanism identified; not causal in posterior hypothalamus

## Why this pass matters

Evidence pass 01 exposed a discordance:

```text
DEL82: delete the rs1421085 cis-regulatory module / ARID motif
  -> no detectable Irx3/Irx5 change in bulk hypothalamus

Exact T>C substitution
  -> Irx3 increases in male posterior hypothalamus
```

A deletion and a single-nucleotide substitution need not be mechanistically equivalent if the nucleotide participates in **multiple factor-binding grammars**.

## Primary evidence: rs1421085 also alters CUX1 binding

Stratigopoulos et al. (J Clin Invest, 2016; PMID 27064284; PMCID PMC4855930) tested the rs1421085 site in the context of CUX1.

Reported observations include:

- rs1421085 is embedded within a CUX1 binding site;
- in gel-shift experiments, the CUX1 P200 form showed preference for the obesity-risk **C** allele;
- the CUX1 P110 form showed higher affinity for the protective **T** allele;
- allele/form-dependent promoter effects were observed in reporter assays involving FTO/RPGRIP1L regulatory constructs.

Primary source:

- https://pmc.ncbi.nlm.nih.gov/articles/PMC4855930/
- PMID 27064284

## Relation to the established adipocyte mechanism

Claussnitzer et al. (NEJM, 2015; PMID 26287746) showed in adipocyte progenitors that T>C disrupts a conserved ARID5B repressor motif and causes loss of repression with increased IRX3/IRX5.

Therefore at least two exact-site regulatory observations exist in different experimental contexts:

```text
rs1421085 T>C
  ├─ adipocyte context: ARID5B repressor binding lost
  └─ CUX1 assays: allele-dependent P110/P200 binding preference changes
```

These observations should **not** be collapsed into one universal mechanism.

## H2a — CUX1 factor/isoform switch

New concrete sub-hypothesis:

> In a relevant posterior-hypothalamic cell state, rs1421085 T>C changes occupancy or functional balance of CUX1 and/or ARID-family factors at the intact element, and this altered factor grammar changes long-range control of Irx3.

Status: `PLAUSIBLE_EXACT_SITE_CANDIDATE`, not yet supported in PH.

## Why H2a explains deletion-vs-substitution discordance

An 82-bp deletion removes the local sequence context entirely. A point substitution preserves the element while changing one base.

Thus:

```text
point mutation:
intact enhancer + altered motif grammar
  -> possible factor switch / competition / contact reweighting

deletion:
local enhancer grammar removed
  -> no opportunity for the same allele-specific switch
```

This is a mechanistic explanation for why `DEL82` need not phenocopy `T>C` in CNS.

It is still only a hypothesis because the 2016 CUX1 experiments did not establish a `CUX1 -> IRX3` causal edge in posterior-hypothalamic cells.

## Immediate computational prerequisites

### P1 — Cux1 presence

Test whether `Cux1` is detectable in the VPH populations that contain `Irx3` and `Arid5b`.

Possible outcomes:

- absent/near-absent in relevant cells -> H2a sharply weakened;
- present but spatially disjoint -> H2a weakened;
- enriched in relevant double-positive cells -> H2a survives and becomes a high-priority wet-lab target.

Standard 3' scRNA-seq may not distinguish functional CUX1 protein forms; RNA presence is only a prerequisite.

### P2 — exact-site motif occupancy

In WT-T vs exact-C PH cells, stratified by sex where possible:

- ARID5B CUT&RUN/ChIP at the orthologous rs1421085 element;
- CUX1 CUT&RUN/ChIP at the same element;
- allele-specific occupancy in heterozygous cells if feasible;
- ATAC footprinting as supporting, not definitive, evidence.

### P3 — perturb the candidate factor

In exact-C PH cells:

```text
CUX1 perturbation
  -> measure Irx3
  -> measure local accessibility/contact
  -> measure neuronal excitability
```

A causal H2a result requires an effect on the SNP-to-Irx3 path, not merely CUX1 binding.

## Updated causal picture

```text
                    ┌─ ARID5B loss? ───────────┐
rs1421085 T>C ──────┼─ CUX1 state switch? ─────┼─ ??? regulatory state
                    ├─ accessibility change? ──┤
                    └─ 3D contact change? ──────┘
                                                   ↓
                                                Irx3 ↑
                                                   ↓
                                      PH neuronal excitability ↓
                                                   ↓
                                           food intake / weight ↑
```

The next objective is not to select one story prematurely. It is to determine which of these candidate edges coexist in the **same PH cell population** and then design the smallest intervention that separates them.
