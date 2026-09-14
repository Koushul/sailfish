# A223 tumor 1-D HIF persistence (E27 / E29)

Applies the E14/E15 Tumor HIF-down 1-D axis (`θ_normoxic`, logistic calibration on E14 Tumor) to Palak-labeled Tumor cells.

- Persistent: `θ ≤ 0.3` (high HIF)
- Partial: `0.3 < θ < 0.7`
- Reverted: `θ ≥ 0.7` (low HIF)
- QC: Palak `Tumor`, spliced UMI ≥ 5000
- θ-only (no RNA velocity). Do not use neutrophil-run `hif_state` (neu-scale).

```bash
/ix1/ylee/kor11/tools/af_tutorial/conda_env/bin/python analysis/a223_tumor_kinetics/run.py
```

Full outputs (including per-cell `tumor_qc.csv` and plots): `/ix1/ylee/kor11/A223/tumor_kinetics/`

Summary tables in `results/` (this folder).

## Result

Pooled A223 QC Tumor **n=9691**: **28.2% persistent**, 15.1% partial, **56.6% reverted** (median θ 0.82).

| Cohort | n | persistent | partial | reverted |
|--------|--:|----------:|--------:|---------:|
| E15S MC38 (same θ-only gates) | 2327 | 26.7% | 16.8% | 56.5% |
| A223 E27+E29 | 9691 | 28.2% | 15.1% | 56.6% |
| E27 | 5066 | 30.5% | 16.3% | 53.2% |
| E29 | 4625 | 25.8% | 13.8% | 60.4% |

OCM gates (DN and lactate+ Tumor n are tiny):

| Gate | n | persistent | reverted |
|------|--:|----------:|---------:|
| DP | 5371 | 31.0% | 53.2% |
| DCF+ (`hypoxia_plus`) | 4283 | 25.0% | 60.7% |
| DN | 25 | 16.0% | 76.0% |
| lactate+ | 12 | 0.0% | 91.7% |

Palak `Tumor (hypoxic)` is the HIF-high subset: **59.3% persistent**, median θ **0.15** (n=1640). Proliferating / epithelial / ribo-high tumors are mostly reverted.

## Why cell cycle, and with vs without it

Glycolytic HIF targets (`Ldha`, `Eno1`, `Pgk1`) are also growth genes. On E14 Tumor, `Ldha` spliced vs Tirosh S/G2M has R² ≈ 0.16. scVelo on E14/E15 tracks S-phase (v_cycle vs S r = 0.50), which is why cycle is scored as a **covariate** and residualized only for velocity — not subtracted from θ. Residualizing those genes from θ would throw out the hypoxia program (it is collinear with growth).

In A223 that collinearity is spatial/niche, not “cycling = HIF-high”: Palak proliferating tumors have the **highest** S score (median 0.15) and are **mostly reverted** (19% persistent). Hypoxic tumors have low S (0.04) and 59% persistent. Growing cells are the oxygenated pool.

```bash
/ix1/ylee/kor11/tools/af_tutorial/conda_env/bin/python analysis/a223_tumor_kinetics/cycle_ablation.py
```

| cohort | model | persistent | reverted | Spearman θ vs S |
|--------|--------|----------:|---------:|----------------:|
| E15S | no cycle (default) | 26.7% | 56.5% | 0.13 |
| E15S | residualize S+G2M on HIF genes | 26.2% | 57.2% | 0.23 |
| E15S | drop genes with cycle R²≥0.1 | 24.5% | 58.5% | 0.23 |
| A223 | no cycle (default) | 28.2% | 56.6% | 0.23 |
| A223 | residualize S+G2M on HIF genes | 33.4% | 48.4% | 0.19 |
| A223 | drop genes with cycle R²≥0.1 | 45.4% | 40.3% | 0.23 |

Dropped genes: `Ldha`, `Pgk1`, `Eno1`, `Angptl4` (the glycolytic HIF core, not a pure cycle list). That jump to 45% persistent is deleting the axis, not “correcting cycle.” Residualizing S+G2M on A223 only moves persist 28% → 33%; every default-persistent cell stays persistent. Default remains **no cycle correction on θ**.

## Combining vs separating neutrophils

Same 14 HIF-down genes. **Separate** = Tumor θ calibrated on E14 Tumor (UMI ≥ 5000), neutrophil θ on E14 neutrophils (UMI ≥ 1500). **Combined** = one μ/σ and logistic anchors from E14 Tumor ∪ neutrophils.

Neutrophils are only **3.8%** of A223 QC Tumor+neu cells (388 / 10,079), so dumping them into the tumor report barely moves the headline percentage. What *does* move the numbers is putting both lineages on one axis.

| cells | calibration | n | persistent | reverted |
|--------|-------------|--:|----------:|---------:|
| A223 Tumor | separate (default) | 9691 | **28.2%** | 56.6% |
| A223 neutrophil | separate (neu E14) | 388 | **49.0%** | 42.3% |
| A223 Tumor | combined Tumor+neu E14 | 9691 | 35.4% | 43.1% |
| A223 neutrophil | combined Tumor+neu E14 | 388 | 30.4% | 57.5% |
| A223 Tumor+neu | pooled on Tumor axis | 10079 | 28.2% | 57.0% |
| A223 Tumor+neu | pooled on neu axis | 10079 | 66.3% | 21.6% |
| A223 Tumor+neu | pooled on combined axis | 10079 | 35.2% | 43.7% |
| A223 Tumor+neu | separate θ, then concat | 10079 | 29.0% | 56.1% |

E15S is the same pattern: Tumor 26.7% persist separate vs 39.6% if you recalibrate on Tumor+neu; neutrophils 20.3% separate vs 9.8% on the combined axis. E14 combined anchors are neutrophil-heavy (363 neu vs 152 Tumor), so tumors look more HIF-high.

Keep lineages **separate**. Do not score tumors on the neutrophil axis.

```bash
/ix1/ylee/kor11/tools/af_tutorial/conda_env/bin/python analysis/a223_tumor_kinetics/lineage_combine.py
```
