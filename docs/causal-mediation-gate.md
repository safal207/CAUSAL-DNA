# Causal Mediation Gate v0.1

CAUSAL-DNA uses the Causal Mediation Gate (CMG) to prevent an associated molecular state from being promoted into a causal bridge merely because it sits between genotype and phenotype in a plausible graph.

The gate is epistemic and executable. It answers:

> What evidence is still missing before candidate mediator M may cross Bardo into material space?

## Core chain

For genotype `G`, candidate mediator `M`, and downstream molecular outcome `Y`:

```text
G -> M
G -> Y
intervene(M) -> Y changes
rescue(M | risk allele) -> Y moves toward reference state
independent verification -> same load-bearing relation survives
```

No individual line substitutes for another.

## State machine

```text
PROJECTIVE
    |
    | genotype->mediator + genotype->outcome
    v
BARDO_ASSOCIATED
    |
    | mediator intervention
    v
BARDO_INTERVENTION_SUPPORTED
    |
    | rescue
    v
BARDO_RESCUE_SUPPORTED
    |
    | independent verification
    v
MATERIAL
```

The final transition is intentionally stricter than a simple association or perturbation result.

## Why intervention is not enough

A mediator can change when the allele changes and still be a passenger. Perturbing that mediator may also change the outcome through a separate route. Rescue is therefore required to test whether normalizing the candidate mediator in the risk-allele background moves the outcome toward the reference state.

Even intervention plus rescue remains in Bardo until independent verification is attached. This preserves CAUSAL-DNA's observer separation: the system that proposes or executes a mechanism cannot silently certify the final causal bridge by itself.

## Invalid transitions

CMG rejects:

- manual declaration of `MATERIAL` without the complete evidence chain;
- rescue evidence without a mediator intervention;
- intervention evidence before genotype-to-mediator and genotype-to-outcome evidence;
- independent verification before rescue;
- `cause_found=true` while no material mediator certificate exists;
- any material mediator while the corresponding gap is still declared `OPEN`.

## CDNA-001 current state

`cases/CDNA-001.mediation-gate.json` tracks the candidate mediators for:

```text
rs1421085 T>C
    -> ???
    -> Irx3 transcription increase in adult male posterior hypothalamus
```

The current candidates are:

- ARID5B occupancy;
- CUX1 occupancy/state;
- exact-element chromatin accessibility;
- exact-element-to-Irx3 3D contact;
- methylation-sensitive local regulatory state / reader occupancy.

They remain `PROJECTIVE` in CMG because current public evidence does not provide the required exact-allele genotype-to-mediator measurement in the matched adult posterior-hypothalamic context. Existing feasibility, wild-type accessibility/contact, motif, transcript, or methylation observations are preserved elsewhere in CAUSAL-DNA but are not silently substituted for this missing contrast.

## Mediator certificate

`causal_dna.mediation_gate.MediationGate.certificate()` returns a machine-readable certificate containing:

- computed status;
- graph space;
- whether materialization is permitted;
- missing evidence requirements;
- provenance references already attached.

Example:

```bash
python -m causal_dna.mediation_gate \
  cases/CDNA-001.mediation-gate.json \
  --mediator M_ARID5B
```

A certificate is not a publication claim. It is a provenance-aware statement about which causal gates have or have not been crossed.

## Scientific invariant

The central invariant is:

```text
association != mediation
intervention != mediation
intervention + rescue != independently verified mediation
```

Only the complete chain may materialize the mediator edge.
