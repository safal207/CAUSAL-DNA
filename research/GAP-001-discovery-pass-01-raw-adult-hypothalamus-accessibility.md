# GAP-001 discovery pass 01 — raw adult hypothalamus accessibility at the exact rs1421085 ortholog

## Question

Does the exact mouse ortholog of human `rs1421085 T>C` show adult hypothalamic chromatin accessibility, or is the element effectively closed outside the previously established adipocyte/developmental contexts?

This pass addresses only the accessibility boundary. It does **not** test allele dependence, TF occupancy, enhancer target, or causal mediation.

## 1. Coordinate gate

The published mouse CRISPR guide from Laber et al. was used as the sequence anchor rather than assuming that a published fosmid interval represented the exact SNP coordinate.

A full `mm10 chr8` scan found one exact match to the reverse complement of the published guide:

```text
mm10 chr8:91374372 (1-based)
mm10 chr8:91374371 (0-based)
T>C on the assembly strand
```

The separately lifted published fosmid starts ~13.27 kb downstream of this exact site. Therefore the fosmid interval is useful as historical/regulatory-region provenance but cannot be substituted for the exact variant coordinate.

Evidence: `discovery/CDNA-001-rs1421085-mouse-ortholog.json`.

## 2. Thresholded candidate-enhancer catalog

The exact coordinate was scanned against the compact `GSE246791` adult whole-brain enhancer/cCRE archive.

After auditing and correcting the archive parser for the documented `chr:start-end` candidate-enhancer format:

```text
chr8 candidate-enhancer intervals parsed: 14,232
chr8 positive-PDC endpoints parsed:        65,700
exact-locus overlaps:                           0
nearest candidate enhancer:               +4,363 bp
```

This is a negative result only for the thresholded catalog layer.

It demotes:

> the exact rs1421085 ortholog is a strong constitutive adult-brain candidate enhancer under the atlas calling thresholds.

It does **not** imply zero raw accessibility.

Evidence: `discovery/CDNA-001-GSE246791-candidate-enhancer-locus.json`.

## 3. Raw 500-bp tile and Tn5 insertion matrices

The atlas authors' per-sample SnapATAC2 H5AD files retain a 500-bp tile matrix and raw insertion matrix. Four hypothalamus dissections were identified from the authors' metadata:

```text
6H = HY-1
7J = HY-2
8K = HY-3
9L = HY-4
```

Two biological replicates were analyzed for each dissection: eight adult male hypothalamus samples total.

The exact variant lies in:

```text
chr8:91374000-91374500
```

Every sample was queried with the same predeclared measurements:

1. non-zero nuclei in the exact 500-bp tile;
2. nuclei with Tn5 insertion signal within +/-250 bp of the exact base;
3. local 21-bin neighborhood for context.

### Replicate result

```text
target 500-bp tile positive:  8 / 8 samples
+/-250 bp insertion positive: 8 / 8 samples
```

Per-region target-tile fractions:

| dissection | replicate range |
|---|---:|
| 6H | 0.160% - 0.212% |
| 7J | 0.158% - 0.186% |
| 8K | 0.126% - 0.201% |
| 9L | 0.204% - 0.230% |

Per-region +/-250 bp insertion fractions:

| dissection | replicate range |
|---|---:|
| 6H | 0.141% - 0.233% |
| 7J | 0.138% - 0.167% |
| 8K | 0.135% - 0.155% |
| 9L | 0.184% - 0.212% |

The exact bin is generally **not** an obvious local accessibility maximum. Target/local-neighborhood median ratios range from approximately `0.68` to `1.54` across samples.

Evidence: `discovery/CDNA-001-GSE246791-hypothalamus-all8.json`; GitHub Actions run `32109378429`, aggregate artifact `9314321285`.

## 4. Discovery statement

### OBSERVED

The exact mouse orthologous rs1421085 element carries low-frequency but reproducible raw chromatin-accessibility signal in every tested adult male hypothalamus dissection and both biological replicates in `GSE246791`, despite not being nominated in the atlas union candidate-enhancer catalog.

### CONTEXT

```text
species: Mus musculus
assembly: mm10
sex: male
age: 8 weeks
tissue: hypothalamus dissections 6H / 7J / 8K / 9L
modality: snATAC
variant state: wild-type atlas; no rs1421085 genotype contrast
```

### NOVELTY / PRIOR BOUNDARY

Prior evidence established the obesity interval, developmental/hypothalamic regulatory relevance, and the exact knock-in phenotype, but did not resolve the adult posterior-hypothalamic molecular regulatory state linking the exact substitution to `Irx3` up-regulation.

This pass contributes a narrower descriptive fact: the orthologous locus is not simply inaccessible in adult male hypothalamus, while also failing the stronger criterion of a prominent constitutive atlas cCRE/peak.

A literature novelty claim remains provisional until a dedicated primary-literature search confirms that this exact raw adult-hypothalamus accessibility observation has not already been reported.

### REPLICATION

Replication is currently **within one dataset only**:

```text
4 anatomical dissections x 2 biological replicates = 8 / 8 positive
```

This strengthens D1 but does not satisfy independent-dataset D2 replication.

### ALTERNATIVES STILL COMPATIBLE

- rare-cell accessibility;
- broadly distributed weak/background accessibility;
- cell-state-specific accessibility diluted at sample level;
- accessibility unrelated to `Irx3` targeting;
- accessibility that is permissive but not causal;
- a T>C-induced change that cannot be observed in a wild-type atlas;
- ARID5B, CUX1, or another interpreter acting only in a subset of accessible nuclei.

### FALSIFICATION / NEXT TEST

The next high-information test is to recover the identities of the signal-carrying nuclei and join them to the atlas integrated subclass/cluster annotations.

Key questions:

1. Is the signal concentrated in a reproducible neuronal or glial subclass?
2. Does that subclass align with an `Irx3`-positive posterior-hypothalamic state from `GSE146692` / HypoMap?
3. Does an independent hypothalamus Multiome dataset reproduce locus accessibility in the corresponding RNA-defined state?

If the signal distributes uniformly across unrelated cell types and does not reproduce in an independent dataset, the mechanistic value of H3 should be strongly demoted.

## 5. Effect on H3_ACCESSIBILITY

### Demoted

```text
H3a: exact locus is completely closed in wild-type adult hypothalamus
H3b: exact locus is a strong constitutive adult-brain cCRE/peak
```

### Still compatible

```text
H3c: rare-cell / subthreshold accessibility
H3d: cell-state-specific accessibility
H3e: T>C changes accessibility relative to wild type
H3f: accessibility is permissive while the causal discriminator is TF occupancy or 3D contact
```

## Status

```text
discovery level: D1_REPLICATED_WITHIN_DATASET_RAW_ACCESSIBILITY
GAP-001: OPEN
cause_found: false
```

This result changes the shape of the search space; it does not close the causal gap.
