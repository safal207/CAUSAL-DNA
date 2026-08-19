# GAP-001 Evidence Pass 03 — developmental gating and multi-branch CNS architecture

Date: 2026-08-17
Status: major search-space refinement; adult exact-SNP PH edge remains unresolved

## Executive finding

The CNS evidence does not support one universal `rs1421085 -> one factor -> one target` mechanism.

Instead, the FTO obesity locus behaves like a **context-gated regulatory system** with at least two experimentally distinguishable CNS programs:

```text
FTO obesity-associated interval
  ├─ developmental hypothalamus / POMC-lineage program
  │    -> long-range IRX3/IRX5 regulation
  │
  └─ adult neuronal program(s)
       ├─ CUX1 -> FTO/RPGRIP1L branch
       └─ exact rs1421085 -> male posterior-hypothalamus Irx3 branch (GAP-001)
```

The last branch remains the target of GAP-001.

---

## Evidence A — long-range IRX3/IRX5 contacts exist in human arcuate-like neurons

Sobreira et al. (Science, 2021; PMID 34083488; PMCID PMC8386003) performed promoter-capture Hi-C in human iPSC-derived hypothalamic arcuate-like neurons.

They reported long-range interactions between the FTO obesity-associated interval and promoters of `IRX3` and `IRX5`. In their maps, the interval did not show the same target pattern toward `FTO/RPGRIP1L`.

Primary source:

- https://pmc.ncbi.nlm.nih.gov/articles/PMC8386003/
- DOI: 10.1126/science.abf1008

Interpretation for GAP-001:

- upgrades H4 (`3D contact`) from generic adult-brain support to **human hypothalamic-neuron-context support**;
- still does not prove that `rs1421085` specifically changes that contact;
- still does not localize the effect to adult posterior-hypothalamic IRX3+ cells.

H4 status: `SUPPORTED_IN_HYPOTHALAMIC_NEURON_CONTEXT_AT_LOCUS_LEVEL`.

---

## Evidence B — deleting the larger obesity interval changes Irx3/Irx5 only during hypothalamic development

Sobreira et al. engineered `mmFtoΔ20`, deleting ~20.2 kb spanning the orthologous obesity-associated interval.

In mouse hypothalamus:

- at embryonic day E17, deletion reduced `Irx3` and `Irx5` expression;
- in adult hypothalamus (~10 weeks), the same deletion did not alter `Irx3/Irx5` expression;
- expression of Irx3/Irx5 was highest during development and fell strongly with maturation.

This is strong evidence for a **temporal gate**.

It also helps reconcile apparent contradictions:

```text
large interval deletion during development
  -> Irx3/Irx5 down

82-bp rs1421085 CRM deletion in young/adult bulk hypothalamus
  -> no detectable Irx3/Irx5 change

exact rs1421085 T>C in adult male posterior hypothalamus
  -> Irx3 up
```

These observations do not need to represent the same regulatory state.

---

## Evidence C — human deletion effect is cell-subtype restricted

Sobreira et al. deleted ~36.1 kb spanning the human FTO obesity-associated locus in iPSCs and differentiated them into hypothalamic arcuate-like neurons.

Single-cell RNA-seq at the neuronal progenitor stage (~91,825 cells) identified multiple developmental populations. The authors reported that reduced `IRX3/IRX5` expression after deletion was detected specifically in a `POMC`-expressing cellular subgroup, rather than uniformly across all cell populations.

Interpretation:

> `region -> target gene` is conditional on **cell identity and developmental state**.

This is directly aligned with the CAUSAL-DNA rule that context belongs inside the causal claim rather than being treated as metadata.

---

## Evidence D — CUX1 exact-site evidence points to a parallel neuronal target program

Stratigopoulos et al. (JCI, 2016; PMID 27064284; PMCID PMC4855930) reported allele-dependent CUX1 P110/P200 binding at `rs1421085` and genotype-associated effects on `FTO`, `RPGRIP1L`, and `AKTIP` in human iPSC-derived neurons.

Crucially, in that neuronal model they **did not observe genotype-associated IRX3 or IRX5 expression**.

Primary source:

- https://pmc.ncbi.nlm.nih.gov/articles/PMC4855930/
- https://www.jci.org/articles/view/85526

This is important negative evidence.

### H2a status correction

Previous status:

`CELLULARLY_FEASIBLE_AND_ENRICHED`

Refined status:

`EXACT_SITE_AND_CELLULARLY_FEASIBLE_BUT_TARGET_GENE_CONFLICT`

Why:

- exact-site CUX1 allele-binding difference: supported;
- Cux1 presence/enrichment in VPH Irx3-containing cells: supported by CAUSAL-DNA T1/P1;
- CUX1 -> adult PH Irx3 causal edge: **not supported**;
- a neuronal model instead connects CUX1/risk alleles to `FTO/RPGRIP1L`, not `IRX3/IRX5`.

Therefore CUX1 may be:

1. a parallel CNS branch;
2. one component of a PH-specific combinatorial grammar that requires other cofactors;
3. irrelevant to the adult PH Irx3 branch despite occupying the same locus.

We must not promote option 2 merely because it is attractive.

---

## New hypothesis H9 — temporal regulatory gating

> The FTO obesity locus changes target-gene routing across development. A developmental enhancer state regulates `IRX3/IRX5` in hypothalamic progenitor/POMC lineages, whereas adult posterior-hypothalamic regulation uses a distinct regulatory grammar that can be reactivated or reweighted by the exact rs1421085 C allele.

Status: `SUPPORTED_AS_ARCHITECTURAL_HYPOTHESIS`.

What is supported:

- developmental vs adult differential sensitivity to interval deletion;
- cell-subtype-specific IRX3/IRX5 response in POMC-expressing precursor cells;
- adult exact-SNP PH Irx3 response is sex- and region-specific.

What is not supported:

- that adult rs1421085 reactivates the developmental POMC enhancer program;
- that the same cell population mediates the adult PH phenotype;
- that CUX1 or ARID5B is the temporal gate.

---

## New hypothesis H10 — multi-target enhancer grammar

> The intact rs1421085-containing interval can route regulatory output to different target genes depending on cell state, developmental stage, sex, and factor occupancy.

Candidate routing states:

```text
State A — adipocyte progenitor
rs1421085 C -> ARID5B repression lost -> IRX3/IRX5 up

State B — neuronal CUX1-sensitive context
rs1421085 C -> CUX1 P200/P110 binding balance -> FTO/RPGRIP1L down

State C — developmental hypothalamic POMC lineage
FTO obesity interval -> 3D contact -> IRX3/IRX5 expression

State D — adult male posterior hypothalamus
rs1421085 C -> ??? -> Irx3 up
```

The scientific target is now to identify the transition rule that makes **State D** distinct.

---

## Highest-value discriminating measurements

The adult PH experiment should no longer measure only one candidate TF.

For the same sorted / spatially defined PH cell population, compare WT-T vs exact-C, male and female:

1. `Irx3` expression at single-cell or spatial resolution;
2. accessibility at the rs1421085-containing element;
3. ARID5B occupancy;
4. CUX1 occupancy / protein-state proxy;
5. enhancer-to-`Irx3` contact;
6. enhancer-to-`Fto/Rpgrip1l` contact;
7. developmental/POMC-lineage markers.

A compact causal discriminator is:

```text
allele
  x sex
  x cell identity
  x developmental/state markers
      -> TF occupancy
      -> chromatin contact
      -> Irx3
```

If the C allele changes `Irx3` without changing CUX1/ARID5B occupancy but changes contact/accessibility, H3/H4 rise.

If occupancy changes without contact change, H1/H2-like local grammar rises.

If only a developmental/POMC-like subpopulation responds, H9 rises sharply.

---

## Current conclusion

We have **not found the missing adult PH molecular mediator**.

But GAP-001 is now much narrower:

> The unresolved cause is likely not a generic FTO-locus mechanism. It is a context-specific regulatory state that distinguishes adult male posterior-hypothalamic cells from both generic neurons and developmental hypothalamic cells.

That is a substantially more falsifiable question than the original `rs1421085 -> ??? -> Irx3` gap.
