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
