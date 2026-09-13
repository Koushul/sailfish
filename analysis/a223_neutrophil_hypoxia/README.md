# A223 neutrophil HIF θ

Script: `python analysis/a223_neutrophil_hypoxia/run.py`

Finds neutrophils in sailfish E27/E29 OCM GEX (Palak labels + S100a8/Cxcr2 markers) and scores the E14/E15 Tumor HIF-down gene set. **Persistent vs reverted gates use E14 neutrophil μ/σ** (same 14 genes / weights as Tumor θ). Tumor-scale θ is stored for a between-lineage check.

Results: `/ix1/ylee/kor11/A223/neutrophil_hypoxia/`
