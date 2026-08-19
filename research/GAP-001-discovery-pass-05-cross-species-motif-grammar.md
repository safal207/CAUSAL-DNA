# GAP-001 discovery pass 05 — human→mouse motif grammar portability

## Question

Does the exact rs1421085 `T>C` substitution perturb the same TF-recognition grammar in the human locus and in the exact mouse ortholog used by the 2025 knock-in model?

This pass is a **sequence-model robustness gate**, not an occupancy or causal claim.

## Provenance

- Human: hg38 `chr16:53767042 T>C` (1-based).
- Mouse: mm10 `chr8:91374372 T>C` (1-based), established independently from the published WT CRISPR guide.
- Human 25-bp context: `CCTAAGGCATGATATTGATTAAGTG`.
- Mouse 25-bp context: `CCTAAGGCATCGTATTGATTAAGCG`.
- Motifs: verified Nature Codebook full atlas, 1,421 unique TFs / 1,685 published motif records.
- Workflow: `GAP-001 human-mouse motif stability`, successful run `32234286303`.

## Why the first cross-species criterion was insufficient

The initial comparator required the same motif record, the same allele-effect direction and `|Δ| >= 0.10` in both species. This produced an apparently strong result: many human and mouse `Δ` magnitudes were identical, including CXXC4, ARID5B and CUX1.

That is not independent replication. Under an additive PWM, a one-base edit at a fixed motif column changes the score according to the allele weights at that column; the surrounding bases can change while the `Δ` contributed by the edited column remains identical.

Therefore identical `Δ` magnitude was explicitly demoted from a replication criterion.

## Strict grammar contract

For a TF to pass the strict human→mouse grammar gate, this pass requires:

1. the most allele-sensitive **human** motif representative is selected first;
2. the exact same `motif_record_id` is evaluated in mouse;
3. the allele-effect direction is conserved;
4. `|Δ| >= 0.10` in both species;
5. the best variant-overlapping **strand and start offset are identical** in human and mouse.

Full-site strength is reported separately. No arbitrary biochemical site-strength cutoff is introduced.

## Global result

Across 1,685 motif records / 1,421 TFs:

- record-level direction conserved: `1591 / 1685 = 94.4%`;
- human-selected TF direction conserved: `1342 / 1421 = 94.4%`;
- human-selected TF placement conserved: `1199 / 1421 = 84.4%`;
- same-record direction + `|Δ| >= 0.10`: `863` TFs;
- strict same-record + same-placement grammar candidates: `760` TFs;
- strict conserved 0.75 threshold crossings: `0`.

The zero conserved threshold crossings reinforces the earlier result that the 0.75 crossing heuristic is source/context sensitive and should not be treated as a primary biological discriminator.

## Named candidates

### CXXC4 — demoted as literal same-grammar mouse mediator

Human:

- `gain_on_C`, `Δ = +0.2888`;
- best placement `strand -`, offset `-3`;
- max site strength `0.6504`.

Mouse:

- `gain_on_C`, `Δ = +0.2888`;
- best placement `strand +`, offset `-2`;
- max site strength `0.6293`.

The best placement changes strand and offset. Therefore:

`strict_grammar_portability = false`.

This **does not rule out CXXC4 biology**. It specifically demotes CXXC4 as a literal conserved motif-placement explanation for the mouse knock-in phenotype. It may remain human-specific or operate through a different alignment/state.

### ARID5B — strict grammar survives, but mouse site is weaker

Human:

- `loss_on_C`, `Δ = -0.2577`;
- `strand -`, offset `-5`;
- max site strength `0.6398`.

Mouse:

- `loss_on_C`, `Δ = -0.2577`;
- `strand -`, offset `-5`;
- max site strength `0.4082`.

`strict_grammar_portability = true`.

This re-elevates ARID5B **relative to CXXC4 for the mouse model**, but it does not remove the adult-hypothalamus evidence gap or the negative tension from the DEL82 bulk-hypothalamus result. Adult-PH occupancy remains unmeasured.

### CUX1 — placement survives, site becomes weak

Human site strength `0.4625`; mouse `0.2390`, with the same `strand - / offset -7` and `loss_on_C Δ = -0.1386`.

`strict_grammar_portability = true`, but the weak mouse site reduces enthusiasm until occupancy is measured.

## New strict candidates worth adult-PH feasibility filtering

Top strict human→mouse candidates include:

| TF | direction | min |Δ| | human site | mouse site | placement |
|---|---|---:|---:|---:|---|
| CXXC1 | gain_on_C | 0.3470 | 0.4371 | 0.5333 | `+ / -4` |
| TET1 | gain_on_C | 0.3054 | 0.6062 | 0.7090 | `+ / -4` |
| KDM2B | gain_on_C | 0.2925 | 0.5888 | 0.5932 | `- / -1` |
| DNMT1 | gain_on_C | 0.2911 | 0.4134 | 0.4918 | `+ / -5` |
| ARID5B | loss_on_C | 0.2577 | 0.6398 | 0.4082 | `- / -5` |
| MECP2 | loss_on_C | 0.2662 | 0.7091 | 0.7581 | `+ / -3` |

The appearance of `TET1`, `DNMT1`, `MECP2`, `CXXC1` and `KDM2B` makes an epigenetic/methylation-sensitive regulatory state a concrete follow-up direction, but the motif scan alone cannot establish such a mechanism.

## Updated discovery frontier

The useful next discriminator is no longer “which motif has the largest Δ?” It is:

`strict grammar portability × adult-PH expression × Irx3-cell compatibility × locus accessibility × methylation/occupancy evidence`.

The immediate next pass should therefore test the strict candidates in adult VPH / posterior-hypothalamic cell states, ideally including co-detection with `Irx3`, before investing in occupancy/contact experiments.

## Causal boundary

No causal edge is promoted by this pass.

`GAP-001 = OPEN`

`cause_found = false`
