# GAP-001 strict broad cell-context screen

Status: **INSUFFICIENT_RESOLUTION — raw expression programs cannot resolve native cell identity**

GitHub Actions run: `31988615816`
Artifact: `gap001-strict-cell-context`
Artifact digest: `sha256:26054a8171475f61a297fc5f99979c517de30eeaa84b308ab05f0eb1e18571d3`
Head SHA: `3e2f3df661792408843c64f97756ae217a3b9bd0`

## Method

Each broad expression program requires at least **2 concordant detected markers**:

- neuronal: `Snap25`, `Syt1`, `Rbfox3`
- astrocyte: `Aqp4`, `Agt`, `Slc1a2`, `Gja1`, `Aldoc`
- oligodendrocyte: `Mbp`, `Mog`, `Plp1`, `Cnp`
- OPC: `Pdgfra`, `Cspg4`, `Olig1`, `Olig2`
- microglia: `Tmem119`, `P2ry12`, `Cx3cr1`, `Aif1`
- endothelial: `Pecam1`, `Cldn5`, `Kdr`, `Emcn`

An exclusive broad label is assigned only when exactly one program passes. Multiple passing programs are `ambiguous`.

## Key result

| candidate subset | n | neuronal program | astrocyte program | exclusive neuronal | exclusive astrocyte | ambiguous |
|---|---:|---:|---:|---:|---:|---:|
| `Irx3+` | 1,791 | 48.30% | 82.19% | 2.23% | 7.82% | **83.70%** |
| `Irx3+/Arid5b+` | 509 | 76.23% | 83.10% | 2.16% | 1.57% | **93.91%** |
| `Irx3+/Cux1+` | 632 | 69.15% | 87.03% | 2.06% | 2.22% | **93.83%** |
| `Irx3+/Arid5b+/Cux1+` | 277 | 82.31% | 87.73% | 1.44% | 0.36% | **97.83%** |

## Interpretation

The earlier high-sensitivity OR-marker screen suggested strong neuronal compatibility of the candidate state. This stricter analysis shows why that observation **cannot be upgraded to a cell-identity claim**: multiple broad programs co-occur in the raw Cell Ranger matrices at implausibly high rates.

The triple-positive candidate state is especially unresolved: 82.3% pass the neuronal program and 87.7% pass the astrocyte program, with 97.8% classified as ambiguous once competing programs are preserved.

Possible contributors include ambient RNA, doublets, broad regional state programs, and the limitations of assigning cell identity from a handful of raw-count markers. Whatever the source, the evidence is insufficient for an exclusive neuron-vs-astrocyte conclusion.

### Causal consequence

We therefore **do not know yet which native PH cell type carries the rs1421085 → Irx3 molecular effect**.

This matters because Sullivan et al. showed that manipulating IRX3 in PH IRX3+ neurons is sufficient to affect neuronal excitability, food intake and body weight, but sufficiency of a downstream neuronal intervention does not by itself prove that the SNP's upstream regulatory action originates in the same neuronal population.

## Methodological decision

`raw-count marker inference` is now marked **insufficient** for GAP-001 cell-of-action resolution.

The next valid evidence tier is one of:

1. published barcode-to-cluster annotations from the Mickelsen/JAX VPH atlas;
2. a curated AnnData object with original cluster labels;
3. spatial transcriptomic localization;
4. exact-SNP single-cell / single-nucleus expression and chromatin data.

Until one of those resolves the population, both neuronal and glial upstream routes remain open.

## Guardrails

- This negative result does not say the cells are biologically mixed identities.
- It says the **screening method cannot separate them reliably**.
- No cell-of-action claim should be promoted from these raw marker screens.
