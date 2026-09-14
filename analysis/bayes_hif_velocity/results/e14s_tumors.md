# E14S tumors (single cohort)

Sample `E14S`, lineage `Tumor`, spliced UMI ≥ 5000. No second library. No forced phenotype. \(h\) is this cohort's own median/MAD. Persist means ≥1.5 MAD above a typical E14S tumor, not an external never-hypoxic arm.

## QC

- n_cells: 152
- median spliced UMI: 14166
- detected panel genes: 20 / 32
- velocity genes: 20

## Phenotype (θ-gate)

- persistent: 0.132 (n=20)
- partial: 0.118
- reverted: 0.750
- GMM persist (argmax): 0.158
- mean p_persist (GMM): 0.186
- mean phenotype entropy (nats): 0.732

Relative ranks: about 1.5 MAD above the cohort median is the persist cut. In a unimodal sample that cut is a tail, not an independent hypoxia class.

## Bayesian lag uncertainty

- σ_v: 1.649 (clip is 1.649)
- median Laplace sd(ξ): 0.376
- mean |ξ|: 0.678
- mean |ξ|/sd: 1.940
- fraction whose 95% ξ interval includes 0: 0.684
- mean p_none / p_away / p_toward: 0.545 / 0.371 / 0.084
- fraction p_away > 0.5: 0.289
- mean flux entropy (nats): 0.684
- mix π none/away/toward: 0.568 / 0.337 / 0.094

Laplace sd(ξ) is the local posterior width of the lag, combined with a spike-slab. If almost every cell's ξ interval covers 0 and p_none is high, the data do not support directed flux.

## Bootstrap uncertainty on h

- n_boot: 80 (cells resampled; ALS factor only)
- median bootstrap sd(h): 0.180
- persist fraction across boots: mean 0.123, 2.5–97.5% 0.085–0.165
- mean per-cell bootstrap P(h≥1.5): 0.127

## Cell cycle

- Spearman (h, S): 0.219
- Spearman (h, G2M): 0.346
- Spearman (h, log L): 0.352
- Spearman (ξ, S): 0.032
- Spearman (ξ, G2M): 0.032
- Spearman (h_cycle, h_nocycle): 0.999
- Spearman (h_nocycle, S): 0.222
- θ-gate label flip with vs without cycle: 0.013
- persist with cycle: 0.132; without: 0.125

Cycle on spliced z is a gene-specific loading (glycolytic HIF targets track S). Cycle on unspliced is a gene-specific intercept for lag. ξ is residualized on S/G2M. Ablation zeros both S and G2M scores.

### Persist by Tirosh S quintile (with cycle in the model)

| s_quintile       |   n |     mean_h |   frac_persist |   mean_p_away |   mean_cycle_s |
|:-----------------|----:|-----------:|---------------:|--------------:|---------------:|
| (-0.305, -0.241] |  31 | -0.187318  |      0.193548  |      0.315954 |     -0.269903  |
| (-0.241, -0.151] |  30 |  0.187761  |      0.266667  |      0.334542 |     -0.1977    |
| (-0.151, 0.0714] |  30 |  0.464469  |      0.133333  |      0.449632 |     -0.0476728 |
| (0.0714, 0.255]  |  30 |  0.248762  |      0.0333333 |      0.412955 |      0.169372  |
| (0.255, 0.554]   |  31 | -0.0447963 |      0.0322581 |      0.342523 |      0.343451  |

## Top factor genes (β)

Slc2a1=0.76, Aldoa=0.74, Eno1=0.73, Tpi1=0.73, Pgk1=0.72, Bnip3=0.71, Egln3=0.69, Gapdh=0.67

## How to read this

- This is one library of tumor cells. Persist/reverted are tails of that library's HIF-target factor.
- Lag posteriors are Bayesian (Laplace + spike-slab). h bootstrap is frequentist resampling of the spliced factor.
- If cycle ablation barely moves h ranks but persist % moves, the 1.5 MAD cut is sensitive, not the factor.

