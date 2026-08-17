# CAUSAL-DNA

**CAUSAL-DNA** is an open research project for turning DNA variant associations into explicit, testable causal chains.

The core question is not only:

> Which variants correlate with a phenotype?

but:

> Through which molecular and cellular transitions does a specific DNA variant cause a measurable phenotype, and where does the evidence stop?

## Core model

```text
Variant
  -> regulatory / molecular effect
  -> target gene or transcript effect
  -> cellular state transition
  -> tissue / circuit effect
  -> organismal phenotype
  -> independent verification
```

Each arrow is a **causal edge**. A causal edge is never promoted to `established` only because two nodes are associated.

Every edge should carry:

- `claim`
- `from` / `to`
- `context` (species, tissue, cell type, sex, developmental state)
- `evidence_type`
- `intervention`
- `source`
- `confidence`
- `alternative_explanations`
- `falsification_test`
- `status`

See [`protocol/causal-proof-v0.1.md`](protocol/causal-proof-v0.1.md).

## First case: CDNA-001 — rs1421085

The obesity-associated non-coding variant **rs1421085 T>C** is a useful benchmark because one branch is unusually well supported experimentally, while another branch remains incomplete.

### Adipocyte branch — strong mechanistic evidence

```text
rs1421085 T>C
  -> disruption of an ARID5B repressor motif
  -> IRX3 / IRX5 derepression in adipocyte progenitors
  -> shift away from thermogenic beige adipocytes
  -> mitochondrial thermogenesis decreases
  -> lipid storage increases
```

The 2015 study used directed perturbations and endogenous CRISPR-Cas9 editing to test this pathway.

### CNS branch — open causal gap

A 2025 mouse study reported that the obesity-associated rs1421085 allele increases brain `Irx3`, including in the posterior hypothalamus (PH), and that increasing IRX3 in the PH is sufficient to increase food intake and body weight. Increased IRX3 also reduced excitability/activity of IRX3-positive PH neurons.

The unresolved edge we start with is:

```text
rs1421085
  -> ???
  -> IRX3 up in posterior hypothalamus
  -> PH neuronal activity down
  -> food intake up
  -> body weight up
```

The project will first test whether the known adipocyte regulatory mechanism generalizes to the relevant CNS cell types, or whether a different regulatory mechanism explains this branch.

See [`cases/CDNA-001-rs1421085.md`](cases/CDNA-001-rs1421085.md).

## Evidence policy

We distinguish:

- **association** — statistical co-occurrence only;
- **supportive** — mechanistically compatible observational evidence;
- **perturbed** — intervention changes the proposed mediator or outcome;
- **rescued** — reversing/restoring the proposed causal factor reverses the effect;
- **established** — convergent evidence with context and alternatives explicitly tested.

A causal claim must preserve negative, conflicting, and superseded evidence instead of deleting it.

## Initial sources

1. Claussnitzer M, et al. *FTO Obesity Variant Circuitry and Adipocyte Browning in Humans.* N Engl J Med. 2015. PMID: 26287746. https://pubmed.ncbi.nlm.nih.gov/26287746/
2. Laber S, et al. *Linking the FTO obesity rs1421085 variant circuitry to cellular, metabolic, and organismal phenotypes in vivo.* Sci Adv. 2021. PMID: 34290091. https://pubmed.ncbi.nlm.nih.gov/34290091/
3. Sullivan AI, et al. *Mice harboring the obesity-associated SNP rs1421085 exhibit increased body weight and reveal an IRX3 neuronal circuit regulating body weight.* Mol Metab. 2025. PMID: 40835181. https://pubmed.ncbi.nlm.nih.gov/40835181/

## Research discipline

CAUSAL-DNA is a computational and evidence-mapping research project. It does **not** treat exploratory hypotheses as medical advice, diagnostic conclusions, or validated therapeutic targets.

---

**Status:** bootstrap / protocol v0.1
