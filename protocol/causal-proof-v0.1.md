# Causal Proof Protocol v0.1

## Purpose

This protocol defines how CAUSAL-DNA represents and tests a proposed causal chain from a DNA variant to a phenotype.

The unit of analysis is a **causal edge**:

```text
A --[evidence / intervention / context]--> B
```

A chain is only as strong as its weakest load-bearing edge.

## 1. Edge record

Every proposed edge must record:

| Field | Meaning |
|---|---|
| `edge_id` | Stable identifier |
| `from` | Upstream state/entity |
| `to` | Downstream state/entity |
| `claim` | Human-readable causal statement |
| `context` | Species, tissue, cell type, sex, developmental/environmental state |
| `evidence_type` | Association, perturbation, rescue, temporal, functional, multi-omic, etc. |
| `intervention` | What was experimentally changed, if anything |
| `outcome_measure` | What was measured |
| `source` | DOI / PMID / dataset / accession |
| `confidence` | `low`, `moderate`, `high` |
| `status` | `candidate`, `supported`, `perturbed`, `rescued`, `established`, `rejected` |
| `alternative_explanations` | Competing mechanisms |
| `falsification_test` | Observation/intervention that would weaken or reject the edge |

Machine-readable representation: [`../schemas/causal-edge.schema.json`](../schemas/causal-edge.schema.json).

## 2. Evidence ladder

### L0 — Association

`A` and `B` covary.

This is **not** sufficient for a causal arrow.

Examples: GWAS association, eQTL association, expression correlation.

### L1 — Contextual support

The proposed mechanism is physically and biologically plausible in the relevant context.

Examples: chromatin accessibility, enhancer/promoter contact, transcription-factor motif, co-expression in the same cell population.

### L2 — Perturbation

Changing `A` changes `B` while relevant context is preserved.

Examples: CRISPR editing, knockdown, overexpression, pharmacological perturbation.

### L3 — Rescue / reversal

Restoring or reversing `A` reverses the effect on `B`.

This is stronger than one-directional perturbation because it tests whether the proposed mediator is load-bearing.

### L4 — Convergent causal support

Independent methods, datasets, models, or laboratories support the same edge and major alternatives have been explicitly challenged.

## 3. Minimum rule for an established edge

An edge SHOULD NOT be marked `established` unless:

1. the direction `A -> B` is explicit;
2. the biological context is recorded;
3. at least one intervention supports the direction;
4. temporal ordering is compatible with the claim;
5. a rescue/reversal or comparably strong causal test exists when feasible;
6. important alternative explanations are recorded;
7. source provenance is preserved;
8. contradictory evidence is retained.

## 4. Context is part of the claim

The following are different claims and must not be silently merged:

```text
variant -> IRX3 up in human adipocyte progenitor
variant -> Irx3 up in mouse posterior hypothalamus
variant -> IRX3 up in a specific hypothalamic neuronal subtype
```

A mechanism demonstrated in one tissue or species cannot automatically fill an edge in another.

## 5. Causal gap

A **causal gap** is a missing or weak edge between two supported nodes.

Example:

```text
rs1421085
  -> [unknown regulatory mechanism in PH]
  -> Irx3 up
```

A gap is scientifically useful when both boundary nodes are supported but the connecting mechanism is unresolved.

## 6. Hypothesis decomposition

For every causal gap, generate competing hypotheses rather than one preferred story.

For `rs1421085 -> Irx3 up in posterior hypothalamus`, initial hypothesis classes include:

- H1: ARID5B motif disruption operates in relevant PH cells.
- H2: another transcription factor binds the rs1421085-containing element in PH cells.
- H3: the variant changes enhancer activity or chromatin accessibility independent of ARID5B.
- H4: the variant changes long-range chromatin contact with `IRX3`.
- H5: the observed `Irx3` increase is indirect through another gene or signaling state.
- H6: the effect is cell-composition-dependent rather than cell-autonomous.
- H7: the effect depends on sex, diet, developmental state, or another interaction.

Each hypothesis must expose a falsification path.

## 7. Negative evidence

Negative results are first-class artifacts.

A failed replication, absent expression in the relevant cell type, lack of enhancer activity, or incompatible temporal ordering must be recorded rather than filtered out.

## 8. Chain score

For exploration only, a chain may be summarized by its weakest edge:

```text
chain_confidence = min(edge_confidence)
```

This prevents several strong downstream edges from hiding one speculative upstream transition.

The score is a triage aid, not a biological truth metric.

## 9. Verification loop

```text
claim
  -> evidence collection
  -> edge classification
  -> competing hypotheses
  -> falsification attempt
  -> independent evidence
  -> status update
```

The project favors **reproducible evidence trails over persuasive narratives**.
