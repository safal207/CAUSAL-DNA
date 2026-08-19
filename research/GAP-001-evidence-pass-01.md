# GAP-001 Evidence Pass 01 — rs1421085 → posterior hypothalamus Irx3

Date: 2026-08-17
Status: evidence triage; no causal gap claimed resolved

## Question

What mechanism connects the obesity-associated `rs1421085 T>C` substitution to increased `Irx3` expression in the posterior hypothalamus (PH)?

```text
rs1421085 T>C
  -> ???
  -> Irx3 ↑ in posterior hypothalamus
  -> PH IRX3+ neuron excitability ↓
  -> food intake ↑
  -> body weight ↑
```

## Finding 1 — the FTO obesity interval physically contacts Irx3 in adult brain

Smemo et al. (Nature, 2014; PMID 24646999; PMCID PMC4113484) used 4C/3C and reported long-range interactions between the obesity-associated interval in `Fto` and the `Irx3` promoter in adult mouse brain. They also reported brain eQTL evidence connecting obesity-associated FTO variants to `IRX3` expression in human brain.

Interpretation:

- supports a **brain-level long-range regulatory route** from the FTO obesity interval to `IRX3`;
- does **not** establish that `rs1421085` itself drives that contact;
- does **not** localize the interaction to posterior-hypothalamic cell types.

H4 (`3D contact`) status after this pass: `SUPPORTED_AT_LOCUS_LEVEL`, not yet causal at variant/PH-cell resolution.

Primary source: https://pmc.ncbi.nlm.nih.gov/articles/PMC4113484/

## Finding 2 — deleting the rs1421085 CRM did not change Irx3/Irx5 in bulk hypothalamus

Laber et al. (Science Advances, 2021; PMID 34290091; PMCID PMC8294759) generated `rs1421085-DEL82`, an 82-nt deletion spanning the conserved rs1421085 cis-regulatory module and the ARID motif.

The deletion changed `Irx3/Irx5` expression in adipocyte progenitors, but the authors reported **no detectable effect on Irx3 or Irx5 expression in the hypothalamus of male or female mice** at 6–8 weeks.

This is important negative evidence against a simple context-independent model:

```text
same rs1421085 CRM / ARID mechanism
  -> Irx3 regulation in every tissue
```

However, it is not a clean falsification of CNS ARID5B involvement because:

1. the perturbation was an 82-bp deletion, not the exact T>C substitution;
2. expression was measured in bulk hypothalamus, which can dilute a posterior-hypothalamus/cell-type-specific effect;
3. physiological context differs from later obesogenic-condition experiments.

H1 (`ARID5B reuse`) status after this pass: `INSUFFICIENT_DATA_WITH_NEGATIVE_TENSION`.

Primary source: https://pmc.ncbi.nlm.nih.gov/articles/PMC8294759/

## Finding 3 — the exact T>C substitution produces a region- and sex-dependent CNS effect

Sullivan et al. (Molecular Metabolism, 2025; PMID 40835181; PMCID PMC12419104) generated mice with the exact homologous substitution corresponding to human `rs1421085 T>C` using CRISPR/Cas9.

They found:

- genotype-dependent increases of `Irx3` across brain regions;
- a notable dose-dependent increase of `Irx3` in the **posterior hypothalamus of male mice**;
- this PH increase was **not observed in female mice**;
- increasing IRX3 specifically in PH IRX3+ neurons was sufficient to increase food intake and body weight;
- increased IRX3 reduced PH IRX3+ neuronal excitability;
- silencing PH IRX3+ neurons reproduced increased food intake/body-weight effects.

Primary source: https://pmc.ncbi.nlm.nih.gov/articles/PMC12419104/

## The key discordance

Two perturbations of nearly the same locus yield different CNS observations:

```text
2021: delete 82-bp CRM incl. ARID motif
      -> no detectable Irx3/Irx5 change in bulk hypothalamus

2025: exact T>C point substitution
      -> Irx3 ↑ in male posterior hypothalamus
```

This discordance is not noise to discard. It is now a **load-bearing clue**.

### New hypothesis H8 — substitution ≠ deletion

`rs1421085 C` may create/reweight a regulatory state that requires the surrounding element to remain intact. Deleting the element can therefore fail to phenocopy the point substitution in CNS.

Candidate mechanisms:

- altered binding competition between ARID5B and another TF;
- gain/loss of a motif whose effect requires neighboring sequence;
- allele-dependent chromatin accessibility;
- allele-dependent enhancer-promoter contact;
- developmental or sex-hormone-dependent cofactor occupancy.

Status: `SUPPORTED_AS_DISCRIMINATING_HYPOTHESIS`, not proven.

## Finding 4 — immediate mouse single-cell dataset is available

The PH single-cell dataset used by Sullivan et al. includes GEO `GSE146692`, a ventral posterior hypothalamus scRNA-seq dataset with >16,000 cells from four juvenile C57BL/6 mouse batches (two male, two female).

This dataset can directly test a prerequisite for H1:

```text
Do Arid5b and Irx3 occur in the same posterior-hypothalamic cells / cell classes?
```

It cannot test rs1421085 genotype effects because these are wild-type reference animals, but it can eliminate mechanisms that require impossible cell-type co-localization.

Dataset: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146692

## Finding 5 — human translation dataset is available

Tadross et al. (Nature, 2025; PMID 39910307; PMCID PMC11922758) created HYPOMAP: 433,369 human hypothalamic nuclei with spatial transcriptomics and 452 hypothalamic cell types.

Use for GAP-001:

- map `IRX3`-expressing human hypothalamic cell types;
- test whether `ARID5B` is co-expressed in candidate cells;
- identify candidate transcriptional cofactors;
- spatially prioritize posterior regions for later genotype-aware experiments.

HYPOMAP itself does not provide a causal rs1421085 perturbation.

Primary source: https://pmc.ncbi.nlm.nih.gov/articles/PMC11922758/
GEO spatial dataset: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278848

## Updated hypothesis ledger

| Hypothesis | Status | Evidence after pass 01 |
|---|---|---|
| H1 ARID5B reuse in PH | `INSUFFICIENT_DATA_WITH_NEGATIVE_TENSION` | DEL82/ARID-motif deletion had no bulk-hypothalamus Irx3/Irx5 effect; PH-specific occupancy remains untested |
| H2 alternative TF | `INSUFFICIENT_DATA` | now more plausible because intact point substitution differs from deletion |
| H3 allele-dependent accessibility | `INSUFFICIENT_DATA` | direct PH ATAC/caQTL evidence not yet identified |
| H4 3D contact | `SUPPORTED_AT_LOCUS_LEVEL` | adult-brain 3C/4C supports Fto obesity interval ↔ Irx3 promoter interaction |
| H5 indirect mediator | `INSUFFICIENT_DATA` | Sullivan found Fto/Rpgrip1l changes inconsistent across regions, but mediator space remains open |
| H6 cell-composition effect | `INSUFFICIENT_DATA` | bulk PH expression cannot distinguish per-cell regulation vs abundance/state shift |
| H7 sex/context interaction | `SUPPORTED_AS_EFFECT_MODIFIER` | PH Irx3 genotype effect seen in males but not females; diet/temperature/context differences remain relevant |
| H8 point-substitution-specific regulatory state | `SUPPORTED_AS_DISCRIMINATING_HYPOTHESIS` | exact T>C and DEL82 perturbations diverge in CNS observations |

## Highest-value next tests

### T1 — PH cell-type feasibility

Use `GSE146692` to quantify:

- fraction of `Irx3+` cells also expressing `Arid5b`;
- neuronal vs astrocytic/tanycytic co-expression;
- male/female differences;
- candidate TFs enriched in `Irx3+` PH populations.

Outcome can demote H1 before any wet-lab experiment.

### T2 — variant-context chromatin test

In WT vs exact-T>C PH tissue/cells, ideally separated by sex:

- ATAC-seq or targeted accessibility assay at the orthologous rs1421085 element;
- H3K27ac/H3K4me1 profiling;
- ARID5B CUT&RUN/ChIP-qPCR;
- candidate-TF occupancy if H1 fails.

### T3 — PH-specific 3D contact

Compare WT vs T>C:

- 3C/4C/Capture-C or H3K27ac HiChIP;
- bait at the rs1421085-containing element;
- target `Irx3` promoter and alternative candidate promoters.

### T4 — decisive intervention

A strong causal discriminator would edit `T↔C` in a relevant PH cell population while leaving the surrounding enhancer intact, then measure:

1. `Irx3` expression;
2. local chromatin accessibility/TF occupancy;
3. neuronal excitability;
4. feeding/weight phenotype where appropriate.

## Current conclusion

GAP-001 is **not resolved**.

But the search space has narrowed substantially. The strongest new clue is that **an intact single-base substitution and deletion of the surrounding CRM are not equivalent in the CNS**. The next question is therefore not simply “does the FTO interval regulate Irx3?” — adult-brain evidence says it can — but:

> What allele-dependent regulatory state exists in posterior-hypothalamic cells only when the rs1421085 element remains intact?
