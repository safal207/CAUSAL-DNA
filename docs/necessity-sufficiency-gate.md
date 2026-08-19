# Necessity & Sufficiency Gate (NSG) v0.2

The Causal Mediation Gate (CMG) asks whether a candidate mediator has crossed the basic intervention/rescue/verification boundary. The Necessity & Sufficiency Gate (NSG) asks a different question:

> What causal role does that mediator actually play in the tested genotype-to-outcome path?

This distinction matters because biological systems are frequently redundant, combinatorial and context-dependent. A mediator can be sufficient without being necessary, necessary without being sufficient, partially load-bearing, or merely correlated.

## Core counterfactual tests

For genotype `G`, mediator `M` and outcome `Y`:

### Necessity

Hold genotype fixed and block/normalize `M`:

```text
G(risk) -> M(blocked/normalized) -> Y ?
```

Interpretation:

- `abolished`: the genotype-associated outcome disappears under mediator blockade;
- `attenuated`: part of the effect remains, supporting partial mediation;
- `no_change`: `M` is not necessary under the tested context.

### Sufficiency

Hold the reference genotype and recreate the mediator state:

```text
G(reference) -> force M(risk-like) -> Y ?
```

Interpretation:

- `recapitulated`: the mediator reproduces the downstream outcome direction;
- `partial`: only part of the downstream outcome is reproduced;
- `no_change`: `M` is not sufficient under the tested context.

## Role matrix

| Necessity result | Sufficiency result | NSG role |
|---|---|---|
| not tested | not tested | `UNRESOLVED` |
| positive | not tested | `NECESSITY_SUPPORTED` |
| not tested | positive | `SUFFICIENCY_SUPPORTED` |
| attenuated | any | `PARTIAL_MEDIATOR` |
| abolished | no change | `NECESSARY_NOT_SUFFICIENT` |
| no change | positive | `SUFFICIENT_NOT_NECESSARY` |
| abolished | partial | `NECESSARY_PARTIALLY_SUFFICIENT` |
| abolished | recapitulated | `NECESSARY_AND_SUFFICIENT` |
| no change | no change | `NON_LOAD_BEARING` |

These labels are always contextual. `NECESSARY_AND_SUFFICIENT` means supported under the declared intervention, cell state, genotype, sex, age and environmental conditions; it is not a universal biological law.

## Residual genotype effect

NSG explicitly records whether the genotype still affects the outcome after mediator blockade.

```text
risk genotype effect after blocking M
    absent  -> compatible with full mediation
    present -> parallel/redundant path survives
```

This prevents a common overclaim: calling `M` *the* mechanism when the genotype still reaches `Y` through another route.

## Independent verification

A role inferred from one experimental program is not silently promoted to a verified causal role.

```text
NECESSARY_AND_SUFFICIENT
        |
        + independent verification
        v
VERIFIED_NECESSARY_AND_SUFFICIENT
```

The same rule applies to partial, necessity-only, sufficiency-only and non-load-bearing classifications.

## CDNA-001

The current `rs1421085 -> Irx3` GAP-001 candidate mediators are all `UNRESOLVED` at NSG v0.2 because the decisive exact-allele necessity/sufficiency interventions have not yet been performed in the relevant adult posterior-hypothalamic cell context.

Candidates currently represented:

- ARID5B occupancy/state;
- CUX1 occupancy/state;
- exact-element accessibility;
- exact-element-to-Irx3 contact;
- methylation-sensitive local regulatory state.

This is intentional. Existing public-data observations and planning scores do not count as necessity or sufficiency evidence.

## Files

- `causal_dna/necessity_sufficiency_gate.py`
- `schemas/necessity-sufficiency-gate.schema.json`
- `cases/CDNA-001.necessity-sufficiency-gate.json`
- `tests/test_necessity_sufficiency_gate.py`

## Scientific boundary

NSG evaluates the internal consistency and provenance of declared experiments. It does not infer biological truth from association, model scores, motif scans or planning priors. Negative tests are first-class outputs: a `NON_LOAD_BEARING` result can be scientifically valuable because it removes a plausible but incorrect causal path.
