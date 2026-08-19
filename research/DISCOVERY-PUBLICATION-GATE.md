# CAUSAL-DNA discovery → publication gate

This gate prevents a computationally interesting signal from being promoted into an "opening" or mechanistic discovery before the evidence supports it.

## Trigger

When a GAP-001 result appears novel relative to the primary literature, evaluate all gates below. If they pass at the appropriate claim level, create a full article package in `safal207/RESONANCE` and prepare an external submission plan.

## Gate A — novelty

- The claim is stated in one falsifiable sentence.
- A fresh primary-literature search finds no prior paper already establishing the same claim in the same biological context.
- Closely related published claims are listed explicitly, including negative/conflicting evidence.

## Gate B — provenance and exact context

- Exact variant/locus coordinates and genome assembly are verified.
- Dataset/sample/barcode provenance is reproducible.
- Cell/tissue/sex/developmental context is explicit.
- Analysis code and machine-readable outputs are committed and CI-verified.

## Gate C — replication / orthogonal support

At least one of:

1. the new signal reproduces in an independent dataset; or
2. the signal is supported by a complementary modality in matched cells (for example RNA + ATAC); or
3. an independent experimental source supports the same load-bearing edge.

A single thresholded peak, single motif score, single pooled expression result, or single biological sample is not enough.

## Gate D — discrimination

- The result distinguishes at least two live competing hypotheses, or materially narrows the mediator/cell-state search.
- Negative results and candidates that fail the gate are preserved.
- No ranking score is treated as biological probability.

## Gate E — claim strength

### Level 1: computational discovery candidate

May be claimed when A–D pass, but no causal language is used. Example shape:

> A reproducible exact-locus cell-state pattern not previously reported.

### Level 2: mechanistic computational evidence

Requires multiple modalities linking locus state to a candidate mediator and target-gene state, with alternative explanations explicitly tested. Still not equivalent to perturbational causality.

### Level 3: causal molecular mechanism

Requires exact-allele perturbation and mediator perturbation/rescue in the relevant biological context, or an equivalently strong causal identification design.

## Automatic RESONANCE action when the gate passes

Create a new RESONANCE issue package, separate from `001-age-of-agents`, with at minimum:

- issue README / editorial framing;
- full research article;
- sources/provenance file;
- claim/evidence table;
- limitations and falsification section;
- reproducibility links back to CAUSAL-DNA commits/workflows;
- external-submission abstract and cover-letter draft.

Do not create the article package merely because a candidate ranks first.

## Current GAP-001 status

- publication gate: **NOT YET PASSED**
- reason: exact adult-PH mediator/occupancy mechanism remains unresolved; current signals are motif, transcript-feasibility and raw-accessibility evidence.
- `GAP-001 = OPEN`
- `cause_found = false`
