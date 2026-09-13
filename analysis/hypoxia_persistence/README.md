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
| κ | `mean(u)/mean(s)` on E14 ∪ most-hypoxic E15 (zeros make median u/s = 0) |
| v_reox | −(u − κ s), then residualized on S/G2M fit in **E14 only**. Same |v| gate on E14: `inducing` = toward hypoxia (control). |
| memory | Muc1 is undetectable; Sod2 z vs E14 |
| states | Velocity first (`reverting` / `inducing`). Else E14 → `never_hypoxic`. E15: `memory` / `reverted` / `persistent` / `partial`. Low-UMI tumor → `tumor_low_umi`. |

Control table: `hypoxia_kinetics_control.csv`. Without GFP, E15 + high θ is **reverted or never-exposed**; E14 is the only clean never-hypoxic class.

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
