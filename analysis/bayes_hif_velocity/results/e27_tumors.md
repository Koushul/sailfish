# E27 tumors (single cohort)

A223 E27 (3′ OCM GEX), lineage `Tumor`, spliced UMI ≥ 5000. Tumor barcodes come from the tumor QC table, not from every droplet. No E14S, E15S, or E29 cells. No neutrophil pool. Image-iT gates (`hypoxia_plus` is the same GFP-channel probe as E15S; `DN` is Image-iT−) are **not** used as control vs exposed. \(h\) is this cohort's own median/MAD among E27 tumors. Persist means ≥1.5 MAD above a typical E27 tumor in this fit. Do not write into the h5ad.

## QC

- n_cells: 5066
- median spliced UMI: 10490
- detected panel genes: 22 / 32
- velocity genes: 22

## Phenotype (θ-gate)

- persistent: 0.085 (n=431)
- partial: 0.204
- reverted: 0.711
- GMM persist (argmax): 0.272
- mean p_persist (GMM): 0.321
- mean phenotype entropy (nats): 0.790

Relative ranks: about 1.5 MAD above the cohort median is the persist cut. In a unimodal sample that cut is a tail, not an independent hypoxia class.

## Bayesian lag uncertainty

- σ_v: 1.649 (clip is 1.649)
- median Laplace sd(ξ): 0.127
- mean |ξ|: 0.248
- mean |ξ|/sd: 1.776
- fraction whose 95% ξ interval includes 0: 0.721
- mean p_none / p_away / p_toward: 0.961 / 0.023 / 0.016
- fraction p_away > 0.5: 0.001
- mean flux entropy (nats): 0.175
- mix π none/away/toward: 0.921 / 0.042 / 0.036

Laplace sd(ξ) is the local posterior width of the lag, combined with a spike-slab. If almost every cell's ξ interval covers 0 and p_none is high, the data do not support directed flux.

## Bootstrap uncertainty on h

- n_boot: 80 (cells resampled; ALS factor only)
- median bootstrap sd(h): 0.027
- persist fraction across boots: mean 0.085, 2.5–97.5% 0.076–0.095
- mean per-cell bootstrap P(h≥1.5): 0.085

## Cell cycle

- Spearman (h, S): -0.236
- Spearman (h, G2M): -0.194
- Spearman (h, log L): -0.082
- Spearman (ξ, S): -0.009
- Spearman (ξ, G2M): -0.002
- Spearman (h_cycle, h_nocycle): 0.999
- Spearman (h_nocycle, S): -0.267
- θ-gate label flip with vs without cycle: 0.012
- persist with cycle: 0.085; without: 0.085

Cycle on spliced z is a gene-specific loading (glycolytic HIF targets track S). Cycle on unspliced is a gene-specific intercept for lag. ξ is residualized on S/G2M. Ablation zeros both S and G2M scores.

### Persist by Tirosh S quintile (with cycle in the model)

| s_quintile         |    n |     mean_h |   frac_persist |   mean_p_away |   mean_cycle_s |
|:-------------------|-----:|-----------:|---------------:|--------------:|---------------:|
| (-0.188, -0.128]   | 1014 |  0.46817   |      0.176529  |     0.0239782 |     -0.151889  |
| (-0.128, -0.0851]  | 1013 |  0.289905  |      0.128332  |     0.0225578 |     -0.10765   |
| (-0.0851, 0.00513] | 1013 | -0.0155417 |      0.053307  |     0.0221075 |     -0.0446534 |
| (0.00513, 0.12]    | 1013 | -0.152458  |      0.0325765 |     0.0224239 |      0.0586525 |
| (0.12, 0.574]      | 1013 | -0.154302  |      0.0345508 |     0.026133  |      0.24569   |

### Image-iT / OCM gates (not used in the fit)

Gates are recorded after fitting. They do not define control vs exposed and do not enter \(h\) or \(\xi\).

| ocm_gate     |    n |    mean_h |   frac_persist |   mean_p_away |   median_spliced_umi |
|:-------------|-----:|----------:|---------------:|--------------:|---------------------:|
| hypoxia_plus | 2535 | -0.017044 |      0.0686391 |     0.0281546 |               8425   |
| DN           |   13 | -0.634241 |      0         |     0.013274  |              20014   |
| DP           | 2518 |  0.195932 |      0.102065  |     0.0187465 |              13598.5 |

## Top factor genes (β)

Bnip3=0.77, Pgk1=0.73, Ndrg1=0.70, Tpi1=0.70, Adm=0.69, Slc2a1=0.68, Ldha=0.65, Gapdh=0.62

## How to read this

- This is one library of tumor cells. Persist/reverted are tails of that library's HIF-target factor.
- Lag posteriors are Bayesian (Laplace + spike-slab). h bootstrap is frequentist resampling of the spliced factor.
- If cycle ablation barely moves h ranks but persist % moves, the 1.5 MAD cut is sensitive, not the factor.

