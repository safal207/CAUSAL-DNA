# GAP-001 discovery pass 06 — adult Irx3+ candidate context

## Question

Do the strict human→mouse motif candidates survive a transcript-availability gate in adult ventral-posterior hypothalamus (VPH) `Irx3+` cells, and does simple sex-dependent transcript presence explain the male-biased phenotype reported for the exact rs1421085 T>C knock-in?

Candidate set fixed before this pass:

- `ARID5B`
- `CUX1`
- `TET1`
- `KDM2B`
- `DNMT1`
- `MECP2`
- `CXXC1`

## Data and provenance

- Public adult mouse VPH scRNA-seq: GSE146692 / Mickelsen et al. 2020 source matrices.
- Biological samples: `AJ18003` male, `AJ18004` female, `AJ19001` female, `AJ19002` male.
- Detection rule: >0 UMI in a filtered Cell Ranger cell.
- Workflow: GitHub Actions run `32234863794`, job `96012459756` — **SUCCESS**.
- Dedicated analysis tests: **5/5 PASS**.
- No expression cutoff or composite causal score was introduced.

## Result 1 — none of the seven candidates is eliminated by transcript availability

Every candidate co-detected with `Irx3` in **all 4/4 biological samples**.

Therefore this screen does **not** support a simple exclusion such as “candidate X is absent from adult VPH Irx3+ cells.”

Replicate-aware order by coverage (all tied at 4/4) then mean candidate detection among `Irx3+` cells:

| candidate | mean across replicates: candidate+ among Irx3+ | pooled candidate+ among Irx3+ | pooled candidate+ among Irx3- | pooled Irx3-context delta |
|---|---:|---:|---:|---:|
| MECP2 | 43.650% | 45.114% | 24.029% | +21.086 pp |
| CUX1 | 33.962% | 35.288% | 20.414% | +14.873 pp |
| ARID5B | 26.366% | 28.420% | 15.636% | +12.784 pp |
| DNMT1 | 22.567% | 23.506% | 11.046% | +12.460 pp |
| CXXC1 | 20.883% | 22.055% | 10.517% | +11.537 pp |
| KDM2B | 17.739% | 17.867% | 8.681% | +9.186 pp |
| TET1 | 14.600% | 15.299% | 8.029% | +7.270 pp |

All seven are more frequently detected in pooled `Irx3+` than pooled `Irx3-` cells in this dataset. This is a descriptive context signal only; it cannot establish TF binding or regulatory direction.

## Result 2 — simple sex-specific transcript availability is not a strong discriminator

Replicate-mean male minus female candidate detection among `Irx3+` cells:

- `ARID5B`: +4.098 pp
- `MECP2`: +3.815 pp
- `CUX1`: +0.688 pp
- `CXXC1`: -0.478 pp
- `TET1`: -0.670 pp
- `KDM2B`: -2.029 pp
- `DNMT1`: -3.893 pp

The observed differences are modest and mixed in direction. With only two biological samples per sex, they are prioritization signals, not genotype-by-sex inference.

Practical implication:

> the male-biased phenotype is unlikely to be explained by a binary “TF transcript present in males, absent in females” prerequisite within this candidate set.

A sex-dependent mechanism may instead live in occupancy, chromatin state, hormone/cofactor state, cell-state composition, developmental history, or downstream response.

## Convergence with strict human→mouse motif grammar

The preceding strict grammar pass required the same published motif record, allele-effect direction, strand and offset in human and mouse.

Within the present seven-candidate set:

- `ARID5B`, `CUX1`, `TET1`, `KDM2B`, `DNMT1`, `MECP2`, and `CXXC1` all survived that strict sequence-grammar gate.
- `CXXC4`, the earlier cross-source motif favorite, was demoted because its best human and mouse placements differed.

The combined result does **not** identify the mediator. It narrows the next useful evidence type: transcript availability is no longer the limiting uncertainty.

## Evidence-priority interpretation

Two mechanistic lanes now deserve explicit separation.

### Lane A — prior exact-site mechanism-linked candidates

`ARID5B` and `CUX1` have prior exact-site biological evidence in other contexts and now also pass adult-VPH transcript feasibility plus strict human→mouse motif grammar.

- `ARID5B`: strongest prior causal positive-control at this locus in adipocyte progenitors, but adult-PH reuse remains unproven.
- `CUX1`: prior allele-dependent exact-site binding evidence, but the relevant target-gene relationship in adult PH remains unresolved.

### Lane B — epigenetic/chromatin-state candidates

`MECP2`, `DNMT1`, `CXXC1`, `KDM2B`, and `TET1` survive both sequence-grammar and adult-Irx3+ transcript-feasibility gates, but currently lack locus-specific adult-PH occupancy evidence.

`MECP2` has the strongest adult `Irx3+` co-detection/context delta in this group, while `TET1`, `KDM2B`, `DNMT1`, and `CXXC1` preserve the broader methylation/chromatin-state hypothesis.

## Next decisive gate

Do **not** spend another pass on generic transcript presence.

Highest-information next measurement is locus/cell-state specific and should compare the two lanes in the same biological context:

`exact T/C allele × adult PH Irx3+ cell state × sex`

with, ideally in the same matched preparation:

1. locus accessibility;
2. allele-specific occupancy for `ARID5B`, `CUX1`, and a prioritized epigenetic candidate such as `MECP2`;
3. local DNA methylation state;
4. enhancer↔`Irx3` contact;
5. `Irx3` RNA;
6. perturb/rescue of any mediator that survives the occupancy gate.

For open-data computational work, the next useful gate is to push from pooled VPH into the most specific available posterior-hypothalamic `Irx3` cell types and overlay locus accessibility, rather than adding more TF-presence screens.

## Causal status

`GAP-001 = OPEN`

`cause_found = false`

No causal edge is promoted by this pass.
