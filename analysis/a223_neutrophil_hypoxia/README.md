# A223 neutrophil HIF θ

Script: `python analysis/a223_neutrophil_hypoxia/run.py`

Finds neutrophils in sailfish E27/E29 OCM GEX (Palak labels + S100a8/Cxcr2 markers) and scores the E14/E15 Tumor HIF-down gene set. **A223 `hypoxia+` is Image-iT LIVE Green ROS (DCF), same FITC-band probe as E15**, used as a hypoxia surrogate; lactate+ is a second stain; DP = both. Persistent vs reverted gates use E14 neutrophil μ/σ (same 14 genes / weights as Tumor θ).

Results: `/ix1/ylee/kor11/A223/neutrophil_hypoxia/`

## Neutrophils in DCF+/hypoxia-sorted gates

There is **no GFP** in A223 or E14/E15 GEX. The FITC-band sort is Image-iT LIVE Green ROS (DCF). vs DN (not DCF+, not lactate+):

Palak-matched cells, Palak `Neutrophil`:

| gate | n | % neutrophil | OR vs DN | p |
|------|--:|-------------:|---------:|--:|
| DN | 902 | **1.22%** | 1 | — |
| DCF+ (`hypoxia_plus`) | 14459 | **5.28%** | 4.52 | 4.9e-10 |
| DP | 12794 | 2.45% | 2.03 | 0.017 |
| DCF+ or DP | 27253 | 3.95% | 3.33 | 2.5e-6 |
| all Palak-matched | 28410 | 3.85% | — | — |

DN matches unsorted A223 PMN-MDSC ~1.22% (Strait et al.). DCF+ over-represents neutrophils ~4.5-fold vs DN. E15S vs E14S placed cells: neutrophils **14.9% vs 8.3%** (OR 1.94).

Table: `neutrophil_gate_enrichment.csv`.
