# Hypoxia persistence and reversion (E14 vs E15)

E14S = never-hypoxic (normoxic). E15S = hypoxic exposure. HIF fate-mapping GFP marks
**history** of hypoxia; the transcriptome can still be hypoxic (**persistent**) or look
like E14 (**reverted**).

This is an iteration branch — no PR until the labels stabilize.

## Literature (paperclip)

Searches: `paperclip search -s pmc` on hypoxia fate-mapping GFP / reoxygenation / scRNA-seq
signatures. Result set `s_0775bf51`, methods `s_34ceb1c5`.

| paper | takeaway for us |
|---|---|
| [PMC8449249](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8449249/) Rocha/Godet — HRE–ODD–Cre switches DsRed→**GFP permanently**. GFP+ = ever-hypoxic. **Currently** hypoxic = GFP+ **and** Hypoxyprobe+. Post-hypoxic = GFP+ Hypoxyprobe−. |
| [PMC7406318](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7406318/) HIF-MARCer — ODD–eGFP is **short-lived** (current HIF); tdTomato after Cre is **permanent** (history). Same two-clock idea. |
| [PMC4688563](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4688563/) Danhier — HRE-EGFP (t½ ~15 h) vs HRE-ODD-luc (t½ ~30 min). EGFP+/luc− = **reoxygenated**. |
| [PMC11872601](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11872601/) 70 hypoxia signatures; Tardon/Buffa/Ragnum work in scRNA-seq. |
| [PMC10334835](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10334835/) Zhang — ssGSEA on hypoxia sets, then **GMM** to call high-confidence hypoxic cells. |
| [PMC9153107](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9153107/) trained hypoxia **classifier** vs simple gene-set scores. |
| [PMC8916770](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8916770/) bulk ssGSEA is biased in cancer scRNA-seq; prefer single-cell scoring. |

Mapping onto this experiment:

| reporter / score | persistent | reverted | never hypoxic |
|---|---|---|---|
| GFP (history) | + | + | − |
| hypoxia transcriptome | high | low (E14-like) | low |
| sample | E15 | E15 | E14 |

No GFP feature is in the sailfish GEX/ADT matrices (no `eGFP` gene, no GFP CITE channel).
Until a GFP column is passed, **all E15 cells are treated as exposed**.

## Cell groups

Four groups on `E14SE15S_gex_adt_placed.h5ad` (RNA only; CITE protein UMIs are essentially empty).

```bash
python analysis/hypoxia_persistence/annotate_cell_groups.py --write-h5ad
```

Writes `cell_groups.csv` and `obs['cell_group']`.

| group | how it is called | n |
|---|---|---|
| Tumor | Ptprc-low clusters, typically high UMI | 3085 |
| TAM | C1q / Csf1r / Cd68 / Adgre1 clusters (not S100a8/Cxcr2) | 4389 |
| neutrophil | S100a8/S100a9/Cxcr2 clusters | 1476 |
| others | monocytes (Ly6c2/Ccr2/Chil3), T/NK, DC, low-UMI | 3118 |

Cycling vs non-cycling tumor are both `Tumor`. Classical monocytes are `others`, not TAM.

## UMAP site

https://cozy-koan-vdzb.here.now/

UMAP of the 12,068 placed cells, colored by `cell_group` (toggle sample). Built from `site_umap/` — **new slug**, does not overwrite other here.now sites.

## Genes expected to move when cells actively revert

Ranked list (mouse symbols for MC38): `reversion_key_genes.tsv`.

Reversion is mostly **HIF-1α protein destroyed by PHD/VHL within minutes of O2**, then HIF target mRNAs fall (hours; Lai peak ~8–12 h). A ROS/NRF2 pulse can **raise** Hmox1/Nqo1. Godet **Muc1** staying high is memory, not reversion.

## Tumor kinetics (v2)

```bash
python analysis/hypoxia_persistence/kinetics.py --write-h5ad
```

Tumor only (`cell_group == Tumor`), spliced/unspliced 1-D model. E14S is the never-hypoxic anchor. Outputs `hypoxia_kinetics.csv`, `hypoxia_kinetics_gene_qc.csv`, and obs columns on the placed h5ad.

| step | what |
|---|---|
| QC | Tumor spliced UMI ≥ 5000 (E14 libraries are ~20× smaller than E15) |
| θ | d-weighted HIF-down genes (log1p size-normalized spliced). Cycle is **not** subtracted from these genes (glycolysis collinear with growth). |
| size factor | spliced and unspliced share the **spliced** UMI total (E14 tumor U/S ≈ 2× E15; separate layer CPM invented fake velocity) |
| κ | `mean(u)/mean(s)` on **E15 persistent** only (lowest 20% θ). E14 is not mixed in. |
| v_reox | kNN-smooth u,s on HIF spliced PCA, then −(u − κ s) scaled by persist MAD (not E14 SD). Center persist median at 0. |
| gate | \|v\| cut = E15 persist 95%. Reverting also needs θ ≤ 0.7; inducing needs θ ≥ 0.3. |
| memory | Muc1 is undetectable; Sod2 z vs E14 |
| states | Phenotype-constrained velocity, then E14 → `never_hypoxic`. E15: `memory` / `reverted` / `persistent` / `partial`. Low-UMI tumor → `tumor_low_umi`. |

v1 of this script z-scored each gene’s residual to **E14** (shallow, different unspliced capture) and gated on the E14 95% \|v\|. That called a handful of high-unspliced outliers `reverting` and hid any real E15 transients under E14 noise. Control table: `hypoxia_kinetics_control.csv` (constrained + unconstrained tails). Without GFP, E15 + high θ is **reverted or never-exposed**; E14 is the only clean never-hypoxic class.

Corrected E15 QC tumor (n = 2327): reverted 717, persistent 704, partial 505, inducing 277, memory 122, reverting **2**. E14 QC (n = 152): never_hypoxic 131, inducing 21, reverting 0. E15 inducing (12%) is **not** above the E14 floor (14%); those cells are high-θ with u > κ_persist s, the expected low-expression artifact of a persist-only κ. Reverting is 2 cells (0.09%). Active HIF transitions are still absent after the estimator fix.

## Neutrophil 1-D model (same Tumor θ genes)

Neutrophil depth is matched E14/E15 (~3.5–4k spliced); min UMI 1500 keeps 1387/1476. Native HIF Cohen’s d vs E14 is ~0 (only Slc2a1 ≥ 0.2), so a neutrophil-refit θ set does not exist. The Tumor module is transferred (`--theta-from hypoxia_kinetics_gene_qc.csv`).

```bash
python analysis/hypoxia_persistence/kinetics.py --cell-group neutrophil \
  --theta-from analysis/hypoxia_persistence/hypoxia_kinetics_gene_qc.csv \
  --out-cells analysis/hypoxia_persistence/neutrophil_kinetics.csv \
  --out-qc analysis/hypoxia_persistence/neutrophil_kinetics_gene_qc.csv \
  --out-control analysis/hypoxia_persistence/neutrophil_kinetics_control.csv
```

E15 QC neutrophils (n = 1024): **reverted 639 (62%)**, persistent 190 (19%), memory 73, partial 93, inducing 20, reverting 9. E14 (n = 363): never_hypoxic 336, inducing 25, reverting 2. Velocity is only P4ha1+Ero1a and is cycle-entangled (b_S = −1.56); treat states as **θ**, not v. Glycolytic Tumor-weighted genes (Eno1 d = −0.47) are lower in E15 neutrophils, so most score E14-like (reverted), unlike Tumor (~30% persistent / 31% reverted).

## Spliced / unspliced quality (is missing velocity a data limit?)

`python analysis/hypoxia_persistence/us_qc.py` → `us_qc_library.csv`, `us_qc_hif_genes.csv`, `us_qc_velocity_panel.csv`.

**The 10x velocyto run is fine. The HIF panel is not.**

| slice | median spliced | median unspliced | U/S | mapping |
|---|---|---|---|---|
| E14 all cells | 6.9k | 1.4k | 0.20 | 0.89 |
| E15 all cells | 9.2k | 2.1k | 0.19 | 0.89 |
| E14 TAM | 7.6k | 1.5k | 0.19 | 0.89 |
| E15 TAM | 7.8k | 1.8k | 0.22 | 0.88 |
| E14 Tumor (all / QC ≥5k) | 2.8k / **14.2k** | 1.1k / 2.3k | 0.28 / 0.20 | 0.80 / 0.90 |
| E15 Tumor (all / QC) | **59.6k / 67.6k** | 7.9k / 8.5k | 0.13 / 0.12 | 0.91 / 0.92 |

Genome-wide this is ordinary 10x (~15% unspliced of total, ~12% ambiguous). E15 Tumor QC has **2355 genes with mean unspliced ≥ 1** — enough for a general RNA-velocity analysis. E14 vs E15 depth is matched for TAM/neutrophil/others; only Tumor is ~20× deeper in E15 (large cells), and 282/434 E14 tumor fail the 5k spliced gate (those failures map at 0.73, i.e. junk, not just shallow).

The kinetic genes are the bottleneck:

| | E14 Tumor QC (n=152) | E15 Tumor QC (n=2327) |
|---|---|---|
| 9 velocity genes, median spliced UMIs | 36.5 | 176 |
| same panel, median **unspliced** UMIs | **1** | **7** |
| cells with **zero** unspliced on the whole panel | **40%** | 4% |

Per-gene E15 means: P4ha1 and Ero1a are the only HIF genes with mean u ≳ 1. Ldha is 82 spliced vs **0.52 unspliced** (κ ≈ 0.006); **ambiguous 7.6** — most intron-overlapping Ldha reads were not assigned to unspliced. Cited2 is 27 spliced and **u = 0 in every cell** (10x 3′ never sees that intron). Car9 / Aldoa / Angptl4 / Ankrd37 were dropped from v for the same reason.

So θ (spliced HIF) is well measured on E15 and adequate on the 152 E14 QC cells. **v is counting ~7 intron UMIs per cell, then kNN-smoothing them.** That can rule out a huge transcription burst; it cannot call a few-percent `reverting` class. The “no active reversion” result is therefore only a weak negative: a real reox wave could hide in Poisson + ambiguous assignment. It is **not** evidence that velocyto failed.

## scVelo vs 1-D model (cycle on/off)

`python analysis/hypoxia_persistence/experiments_scvelo.py`

scVelo stochastic on all 12,068 cells (HVGs ∪ HIF ∪ cycle genes). Cycle correction residualizes spliced+unspliced on S/G2M before moments (all cells, or E14-fit only).

| run | v_cycle vs S | v_HIF vs S | \|v_cycle\|/\|v_HIF\| | % cells cycle dominates | vs our v_reox (tumor) |
|---|---|---|---|---|---|
| scVelo, no CC | **0.50** | 0.14 | 0.52 | **25%** | r = 0.17 |
| scVelo, CC (all cells) | 0.25 | 0.05 | 0.31 | 12% | r = 0.15 |
| scVelo, CC (E14 fit) | 0.39 | 0.15 | 0.33 | 14% | r = 0.17 |
| our 1-D, CC on v | — | v_reox vs S **0.006** | — | — | (self) |
| our 1-D, CC off | — | v_reox vs S −0.016 | — | — | vs CC-on **r = 0.986** |

scVelo on all cells is a **cell-cycle velocity** (v_cycle tracks S-phase at r = 0.50). Residualizing cycle cuts that in half but does not make scVelo match the HIF 1-D model (r ≈ 0.15–0.17). Our v_reox is already nearly orthogonal to cycle; turning cycle residualization off barely changes it. E14 `inducing` is the noise floor for “moving toward hypoxia.”

## Run

```bash
python analysis/hypoxia_persistence/hypoxia_states.py \
  --h5ad /ix1/ylee/kor11/MC38/E14SE15S/E14SE15S_gex_adt_placed.h5ad \
  --out analysis/hypoxia_persistence/hypoxia_states.csv
```

If GFP (or ODD-GFP ADT) exists as an obs column:

```bash
python analysis/hypoxia_persistence/hypoxia_states.py \
  --h5ad /ix1/ylee/kor11/MC38/E14SE15S/E14SE15S_gex_adt_placed.h5ad \
  --gfp-obs GFP \
  --out analysis/hypoxia_persistence/hypoxia_states.csv
```

## What v1 does

1. Score Hallmark / Buffa / Winter hypoxia and mouse HIF-1α targets (mean log1p of matched genes).
2. Logistic regression: E14 vs E15 cells in the **top quartile** of Hallmark hypoxia → `P_hypoxic`.
3. 3-component GMM on the Hallmark axis within the data; E15 low / mid / high → `reverted` / `partial` / `persistent`. E14 → `never_hypoxic`.
4. Optional GFP gate: GFP− E15 cells are `unlabeled_*` instead of reverted/persistent.

Outputs: `hypoxia_states.csv`. Only input is `E14SE15S_gex_adt_placed.h5ad`.

## Did the E15 Image-IT / GFP sort scoop TAM via autofluorescence?

Image-iT LIVE Green is carboxy-H2DCFDA; the oxidized product is a fluorescein (FITC/GFP band, ~495/529 nm). Literature on **unstained** FITC autofluorescence, not DCF biology:

| FITC AF (high → low) | evidence | our E15 vs E14 |
|---|---|---|
| **Eosinophils** | FAD/FMN granules emit ~520 nm; used to FACS-sort unlabeled eos ([Weil 1981](https://pubmed.ncbi.nlm.nih.gov/7460387/); [PMC11617454](https://pmc.ncbi.nlm.nih.gov/articles/PMC11617454/)) | top 5% eos-score **depleted** (OR 0.64) |
| **Tissue macrophages / TAM** | highest myeloid AF after eos; flavin/NADH in FITC; alveolar/Kupffer AF is a named unmixing component ([Abcam](https://www.abcam.com/en-us/technical-resources/applications/flow-cytometry/flow-cytometry-buffers-reagents-equipment/autofluorescence-in-flow-cytometry); [Mitchell 2010 JLB](https://doi.org/10.1189/jlb.0310184); [PMC8965042](https://pmc.ncbi.nlm.nih.gov/articles/PMC8965042/)) | TAM **depleted** (OR 0.66) |
| **Neutrophils** | moderate FITC AF, weaker than eos; AF-based neutrophil sorts exist ([Dorward 2013 JLB](https://doi.org/10.1189/jlb.0113040); [Yakimov 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6701549/)) | **enriched** (OR 1.94) |
| **Tumor / epithelium** | elevated, variable metabolic AF ([Smith 2006 Cytometry](https://doi.org/10.1002/cyto.b.20090)) | **enriched** (OR 5.88) |
| **Monocytes** | > lymphocytes, << tissue macs | mono-like others **depleted** (OR 0.29) |
| **Lymphocytes (T/NK)** | lowest FITC AF; PBMC unmixing treats them as the dim signature | tnk-like others **depleted** (OR 0.39) |

`python analysis/hypoxia_persistence/ros_sort_tam.py` and `ros_sort_fitc_literature.csv`.

If the GFP gate were dominated by **autofluorescence**, the literature ranking predicts eosinophils and TAM in, lymphocytes out. Lymphocytes are out, but the two highest-AF populations (eos, TAM) are **not** in — they are depleted. What is in is Tumor and neutrophils, which is the ranking for **true DCF/ROS** (hypoxic tumor ROS + neutrophil burst), not the ranking for FITC AF.

Caveat: TAM are still 32% of E15, so AF can contribute at the margin. It does not explain the sort. There are almost no SiglecF/Epx-high cells in either sample, so eosinophil AF is not a hidden E15 contaminant.

## Public validation of θ

`analysis/hypoxia_persistence/validation/` scores human orthologs of the Tumor HIF module
on GEO accessions (priority: GSE200207, GSE296547, GSE227508 skipped ~3 GB + BAM,
GSE240212 Visium Moran I, GSE292771 HIF KO, GSE30019 bulk reox). Velocity and GFP
history are not testable on these files. Write-up and numbers: `validation/README.md`.
