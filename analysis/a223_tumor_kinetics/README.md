# A223 tumor 1-D HIF persistence (E27 / E29)

Applies the E14/E15 Tumor HIF-down 1-D axis (`θ_normoxic`) to Palak-labeled Tumor cells in the A223 technical replicates.

- Persistent: `θ ≤ 0.3` (high HIF / not reverted)
- Reverted: `θ ≥ 0.7`
- QC: spliced UMI ≥ 5000
- θ-only (no RNA velocity)

Run:

```bash
conda run -p /ix1/ylee/kor11/tools/af_tutorial/conda_env \
  python analysis/a223_tumor_kinetics/run.py
```

Outputs: `/ix1/ylee/kor11/A223/tumor_kinetics/`
