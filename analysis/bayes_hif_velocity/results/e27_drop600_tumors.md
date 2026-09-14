# E27 drop-top-600 tumors (single cohort)

Control: A223 E27 tumors after deleting the 600 highest-\(h\) cells from the previous single-cohort fit (`e27_tumor_cells.tsv`). Removed \(h\) range 1.29–3.38; 431 of those 600 were θ-gate persistent in the first fit (first-fit persist n=431). Remaining n=4466. The model is fit from scratch on the remaining cells only: new gene detection, new ALS \(h\), new MAD, new lag. It is not told that any cells were removed. Image-iT gates still do not enter the fit. Do not write into the h5ad.

## QC

- n_cells: 4466
- median spliced UMI: 10739
- detected panel genes: 22 / 32
- velocity genes: 22

## Phenotype (θ-gate)

- persistent: 0.032 (n=144)
- partial: 0.247
- reverted: 0.721
- GMM persist (argmax): 0.299
- mean p_persist (GMM): 0.276
- mean phenotype entropy (nats): 0.607

Relative ranks: about 1.5 MAD above the cohort median is the persist cut. In a unimodal sample that cut is a tail, not an independent hypoxia class.

## Bayesian lag uncertainty

- σ_v: 1.649 (clip is 1.649)
- median Laplace sd(ξ): 0.135
- mean |ξ|: 0.256
- mean |ξ|/sd: 1.743
- fraction whose 95% ξ interval includes 0: 0.734
- mean p_none / p_away / p_toward: 0.957 / 0.027 / 0.016
- fraction p_away > 0.5: 0.001
- mean flux entropy (nats): 0.185
- mix π none/away/toward: 0.918 / 0.045 / 0.036

Laplace sd(ξ) is the local posterior width of the lag, combined with a spike-slab. If almost every cell's ξ interval covers 0 and p_none is high, the data do not support directed flux.

## Bootstrap uncertainty on h

- n_boot: 80 (cells resampled; ALS factor only)
- median bootstrap sd(h): 0.028
- persist fraction across boots: mean 0.030, 2.5–97.5% 0.016–0.046
- mean per-cell bootstrap P(h≥1.5): 0.030

## Cell cycle

- Spearman (h, S): -0.153
- Spearman (h, G2M): -0.131
- Spearman (h, log L): -0.028
- Spearman (ξ, S): -0.004
- Spearman (ξ, G2M): -0.004
- Spearman (h_cycle, h_nocycle): 0.999
- Spearman (h_nocycle, S): -0.194
- θ-gate label flip with vs without cycle: 0.026
- persist with cycle: 0.032; without: 0.026

Cycle on spliced z is a gene-specific loading (glycolytic HIF targets track S). Cycle on unspliced is a gene-specific intercept for lag. ξ is residualized on S/G2M. Ablation zeros both S and G2M scores.

### Persist by Tirosh S quintile (with cycle in the model)

| s_quintile         |   n |     mean_h |   frac_persist |   mean_p_away |   mean_cycle_s |
|:-------------------|----:|-----------:|---------------:|--------------:|---------------:|
| (-0.197, -0.131]   | 894 |  0.224661  |      0.0425056 |     0.0298408 |     -0.15761   |
| (-0.131, -0.0823]  | 893 |  0.108952  |      0.0347144 |     0.024673  |     -0.109481  |
| (-0.0823, 0.00907] | 893 | -0.0437571 |      0.0347144 |     0.0254557 |     -0.0390118 |
| (0.00907, 0.122]   | 893 | -0.136185  |      0.0302352 |     0.0255697 |      0.0607435 |
| (0.122, 0.565]     | 893 | -0.130819  |      0.019037  |     0.0281427 |      0.245537  |

### Image-iT / OCM gates (not used in the fit)

Gates are recorded after fitting. They do not define control vs exposed and do not enter \(h\) or \(\xi\).

| ocm_gate     |    n |     mean_h |   frac_persist |   mean_p_away |   median_spliced_umi |
|:-------------|-----:|-----------:|---------------:|--------------:|---------------------:|
| hypoxia_plus | 2292 | -0.0595132 |      0.0309773 |     0.0318708 |               8600.5 |
| DN           |   13 | -0.568279  |      0         |     0.0146353 |              20014   |
| DP           | 2161 |  0.0760869 |      0.0337807 |     0.021365  |              13943   |

## Top factor genes (β)

Bnip3=0.72, Tpi1=0.70, Pgk1=0.67, Ldha=0.61, Ndrg1=0.60, Slc2a1=0.57, Aldoa=0.55, Adm=0.53

## How to read this

- This is one library of tumor cells. Persist/reverted are tails of that library's HIF-target factor.
- Lag posteriors are Bayesian (Laplace + spike-slab). h bootstrap is frequentist resampling of the spliced factor.
- If cycle ablation barely moves h ranks but persist % moves, the 1.5 MAD cut is sensitive, not the factor.

## Versus the untrimmed E27 fit

- Removed 600 cells with highest first-fit \(h\) (min removed \(h\)=1.288). Overlap of remaining barcodes with the removed list: 0.
- First-fit persist fraction among all 5066: 0.085. This refit's persist fraction: 0.032.
- If persist is only the upper tail of a unimodal factor, trimming the old tail and re-standardizing should grow a new tail of similar size.
