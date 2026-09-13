# Public validation of the Tumor HIF θ module

The E14/E15 1-D model’s **θ** is a Cohen’s-d-weighted HIF-down spliced module
(mouse: Car9, Vegfa, Ldha, Pdk1, Egln1, Pgk1, Eno1, Aldoa, Angptl4, P4ha1,
Ero1a, Ankrd37, Serpine1, Cited2). These scripts score the **human orthologs**
on labeled public matrices. They do **not** re-estimate unspliced velocity
(GSE200207 BAMs are withheld; GSE227508 velocyto was not rebuilt).

GEO archives stay in `data/` (gitignored). Commit only code and `results/*.csv`.

```bash
cd analysis/hypoxia_persistence/validation
python3 robustness.py
python3 -m pytest test_score.py -q
```

Seurat RDS → gene subset (once, needs Seurat 5):

```bash
Rscript --vanilla extract_seurat.R data/GSE200207_seurat.rds data/extracted/GSE200207 hif_genes.txt
Rscript --vanilla extract_seurat.R data/GSE296547_seurat.rds data/extracted/GSE296547 hif_genes.txt
```

`run_public.py` parses GSE292771 TCM, GSE30019 series matrix, and GSE240212 MTX.

## What each accession can test

| priority | accession | θ (Hx vs Nx) | reox / time | velocity | GFP history | used |
|---|---|---|---|---|---|---|
| 1 | GSE200207 | normal kidney 0.5% vs 21% O2; ccRCC at 21% as constitutive HIF | no | no (BAMs withheld) | no | yes |
| 2 | GSE296547 | Parse AT2 hypoxia days vs NOR | HYPO30 then NOR 6/15/30 d | no (no velocyto) | no | yes |
| 3 | GSE227508 | 10x HBE 1% O2 6 h vs 5 d | duration, not reox | Seurat RDS ~3 GB + SRA BAMs; not downloaded | no | **not run** |
| 4 | GSE240212 | Visium HIF spatial structure | no | no | GFP/RFP images **not** in the MTX tar | Moran I only |
| 5 | GSE292771 | MCF-7 48 h Hx vs Nx; HIF-1/2 KO | no | no | no | yes |
| extra | GSE30019 | bulk MCF-7 reox 0–24 h | yes | no | no | yes |

Robustness on every labeled split: AUROC, 400-bootstrap 95% CI, 400 label-permutation p,
Cohen’s d, MWU, module vs median gene, leave-one-gene-out. Visium uses **hex** neighbors
(rook 4-neighbors on the Visium grid find almost no pairs).

## Results

| test | metric | 95% bootstrap CI | perm p |
|---|---|---|---|
| GSE200207 normal kidney 0.5% vs 21% | AUROC **0.768** (d = 0.97) | 0.746–0.789 | 0.0025 |
| GSE200207 vs median gene | median gene 0.579; best PGK1 0.836; worst ENO1 **0.183** | — | — |
| GSE200207 LOO min (drop PGK1) | AUROC 0.655 | — | — |
| GSE200207 drop ENO1 | AUROC **0.835** (ENO1 anti-correlates in kidney) | — | — |
| GSE200207 ccRCC 21% vs normal 21% | AUROC **0.851** (d = 1.27) | 0.836–0.867 | 0.0025 |
| GSE296547 all cells HYPO vs NOR | AUROC 0.636 | 0.621–0.650 | 0.0025 |
| GSE296547 AT2 by timepoint HYPO vs NOR | AUROC 0.592 | 0.573–0.609 | 0.0025 |
| GSE296547 AT2 HYPO30 vs reox | AUROC **0.690** | 0.666–0.714 | 0.0025 |
| GSE296547 all-cell HYPO30 vs reox | AUROC 0.579 (basal mix) | 0.563–0.594 | 0.0025 |
| GSE296547 AT2 median vs hypoxia day | Spearman **1.0** (4 timepoints) | — | 0.087 |
| GSE292771 WT Hx vs Nx | AUROC **0.821** vs median gene 0.533 | 0.812–0.829 | 0.0025 |
| GSE292771 HIF1KO / HIF2KO | AUROC 0.622 / 0.657 | 0.610–0.633 / 0.646–0.667 | 0.0025 |
| GSE292771 WT LOO min (drop LDHA) | AUROC 0.786 | — | — |
| GSE30019 0 h vs 24 h (n = 6) | AUROC 1.0; Spearman vs hours **−0.969** | 1–1 | AUROC 0.075; Spearman 0.0025 |
| GSE240212 Tumor-897 / 899 Moran I | **0.51 / 0.66** | — | 0.005 |

θ transfers as a hypoxia / constitutive-HIF score and falls after reox in bulk MCF-7 and in AT2. It is not a 1-gene score: ENO1 is inverted in kidney; PGK1/LDHA carry most of the signal. HIF-1 KO knocks the MCF-7 AUROC from 0.82 to 0.62. All-cell Parse reox looks weak because timepoints change cell-type mix.

## What this does *not* validate

- The unspliced **v** estimator and `reverting`/`inducing` gates (needs matched
  velocyto + known reox time or Hypoxyprobe).
- GFP history vs transcriptome (Godet images not in GEO MTX).
- Mouse MC38-specific weights on every human tissue equally (ENO1 is inverted in
  GSE200207 normal kidney; the module can still work if other genes compensate).
