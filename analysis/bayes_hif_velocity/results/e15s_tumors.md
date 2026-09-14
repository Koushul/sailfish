# E15S tumors (single cohort)

Sample `E15S`, lineage `Tumor`, spliced UMI ≥ 5000. No second library. No forced phenotype. Parameters are estimated from these cells only. \(h\) is this cohort's own median/MAD. Persist means ≥1.5 MAD above a typical E15S tumor, not a comparison to another sample.

## QC

- n_cells: 2327
- median spliced UMI: 67639
- detected panel genes: 24 / 32
- velocity genes: 23

## Phenotype (θ-gate)

- persistent: 0.139 (n=324)
- partial: 0.188
- reverted: 0.673
- GMM persist (argmax): 0.231
- mean p_persist (GMM): 0.363
- mean phenotype entropy (nats): 0.844

Relative ranks: about 1.5 MAD above the cohort median is the persist cut. In a unimodal sample that cut is a tail, not an independent hypoxia class.

## Bayesian lag uncertainty

- σ_v: 1.649 (clip is 1.649)
- median Laplace sd(ξ): 0.152
- mean |ξ|: 0.308
- mean |ξ|/sd: 1.832
- fraction whose 95% ξ interval includes 0: 0.666
- mean p_none / p_away / p_toward: 0.895 / 0.078 / 0.027
- fraction p_away > 0.5: 0.034
- mean flux entropy (nats): 0.296
- mix π none/away/toward: 0.865 / 0.089 / 0.046

Laplace sd(ξ) is the local posterior width of the lag, combined with a spike-slab. If almost every cell's ξ interval covers 0 and p_none is high, the data do not support directed flux.

## Bootstrap uncertainty on h

- n_boot: 80 (cells resampled; ALS factor only)
- median bootstrap sd(h): 0.039
- persist fraction across boots: mean 0.138, 2.5–97.5% 0.120–0.154
- mean per-cell bootstrap P(h≥1.5): 0.138

## Cell cycle

- Spearman (h, S): -0.139
- Spearman (h, G2M): -0.181
- Spearman (h, log L): 0.117
- Spearman (ξ, S): -0.024
- Spearman (ξ, G2M): -0.038
- Spearman (h_cycle, h_nocycle): 0.998
- Spearman (h_nocycle, S): -0.147
- θ-gate label flip with vs without cycle: 0.014
- persist with cycle: 0.139; without: 0.142

Cycle on spliced z is a gene-specific loading (glycolytic HIF targets track S). Cycle on unspliced is a gene-specific intercept for lag. ξ is residualized on S/G2M. Ablation zeros both S and G2M scores.

### Persist by Tirosh S quintile (with cycle in the model)

| s_quintile        |   n |    mean_h |   frac_persist |   mean_p_away |   mean_cycle_s |
|:------------------|----:|----------:|---------------:|--------------:|---------------:|
| (-0.448, -0.146]  | 466 | 0.722681  |      0.33691   |     0.108442  |     -0.266387  |
| (-0.146, -0.0306] | 465 | 0.213688  |      0.139785  |     0.0828072 |     -0.0869576 |
| (-0.0306, 0.0591] | 465 | 0.0726866 |      0.0903226 |     0.0721547 |      0.0143497 |
| (0.0591, 0.154]   | 465 | 0.0698877 |      0.0602151 |     0.0558523 |      0.103861  |
| (0.154, 0.559]    | 466 | 0.0761598 |      0.0686695 |     0.0707846 |      0.235201  |

## Top factor genes (β)

Bnip3=0.76, Tpi1=0.72, Pgk1=0.71, Bnip3l=0.71, Pkm=0.71, Ldha=0.65, Aldoa=0.63, Slc2a1=0.62

## How to read this

- This is one library of tumor cells. Persist/reverted are tails of that library's HIF-target factor.
- Lag posteriors are Bayesian (Laplace + spike-slab). h bootstrap is frequentist resampling of the spliced factor.
- If cycle ablation barely moves h ranks but persist % moves, the 1.5 MAD cut is sensitive, not the factor.

