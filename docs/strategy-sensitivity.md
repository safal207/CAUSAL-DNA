# Strategy Robustness and Planning-Prior Sensitivity

A strategy can look optimal only because a particular planning prior was chosen. CAUSAL-DNA therefore separates **strategy selection** from **strategy robustness**.

This layer asks:

> Across a broad, explicit region of alternative planning priors, how often does the same first experiment or conditional policy remain optimal?

It does **not** ask how probable a biological hypothesis is.

## Deterministic prior grid

For CDNA-001 we use four open planning hypotheses and enumerate positive integer compositions of 21 units with at least one unit per hypothesis.

```text
H1 + H2 + H3 + H4 = 21
Hi >= 1
```

This produces exactly 1,140 distinct prior states. The grid is deterministic, so CI replays the same sensitivity region on every run.

The smallest allowed share is `1/21 ≈ 4.8%`; no hypothesis is silently assigned zero plausibility in this first pass.

## What is measured

For every prior state the multi-step strategy planner is rerun under the same budget, depth, experiment costs, outcome planning probabilities, and likelihood-like multipliers.

We count:

- which **first experiment** wins;
- the exact **first-step policy signature**, including the contingent next experiment for each possible first outcome;
- support of the baseline strategy chosen under equal initial planning weights;
- the dominant first step and dominant policy over the whole grid.

Example policy signature:

```text
E3_ALLELE_SNATAC|
  accessibility_shift->E5_MATCHED_MULTIOMIC,
  no_accessibility_shift->E1_ARID5B_OCCUPANCY
```

## Operational robustness thresholds

`cases/CDNA-001.sensitivity-analysis.json` declares explicit thresholds for labeling planning stability:

- `first_step_robust_threshold`
- `policy_signature_robust_threshold`

These thresholds are engineering/planning criteria. They are **not** confidence intervals, p-values, posterior probabilities, or biological evidence.

A strategy that fails a threshold is not biologically false; it means the planning recommendation is sensitive to the assumed prior weights and should be presented as such.

## Scope boundary

This first sensitivity pass varies **planning priors only**. It holds outcome probabilities, likelihood multipliers, experiment costs, and assay definitions fixed.

A later robustness layer can independently vary those model assumptions. Keeping the dimensions separate makes it possible to identify whether instability comes from:

1. hypothesis priors;
2. outcome-model assumptions;
3. experimental costs;
4. the causal graph itself.

## Integrity rule

Sensitivity output is advisory only. It cannot set:

- `cause_found`;
- causal edge status;
- Bardo materialization;
- verification status.

Changing 1,140 planning priors does not create a single new biological observation.

## Implementation

- `causal_dna/sensitivity.py`
- `schemas/sensitivity-analysis.schema.json`
- `cases/CDNA-001.sensitivity-analysis.json`
- `tests/test_sensitivity.py`
- `analysis/validate_repository.py`
