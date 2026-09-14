# Critique: circularity, consistency, red flags

This is an audit of the HIF-α lag **method** and the E14/E15 application. Dataset names below are only as **design labels** (never-hypoxic vs hypoxia-exposed). They are not part of the likelihood.

## What is consistent

- **Do not pool lineages.** Tumor vs neutrophil panels, UMI floors, and *Car9* loading differ; a shared \(\kappa,\lambda\) would mix scales.
- **Cell-wide \(L\) + CPM + within-sample \(\ell_n\).** E15 tumors are \(\sim 5\times\) deeper than E14. After that correction, \(\mathrm{corr}(h,\log L)\approx 0.10\) in tumors (was \(\sim 0.34\)). CA9/VEGFA/BNIP3 CPM are **not** induced in E15; glycolysis is a modest \(\sim 1.2{-}1.3\times\). Persist after the fix is a small tail (\(7.6\%\) tumors, \(6.1\%\) neutrophils), not a 30% block. That matches the biology of the panel more than the historical depth-mixed \(\theta\).
- **Lag is weak on real data.** Tumor \(p^{\mathrm{away}}\) is \(0.16\) (E14) vs \(0.14\) (E15). There is no E15-specific reoxygenation wave. Synthetic tests already said hard `transitioning_out` is under-called; trust probabilities. \(\sigma_v\) sitting on the \(e^{0.5}\) clip is a warning that lag scale is not identified, not a feature.
- **Cycle on unspliced only.** Spearman \((\xi, S)\approx -0.03\) (tumors). Do not residualize \(\theta\) with one global S/G2M regression.
- **Do not treat DCF+ as HIF.** Image-iT LIVE Green ROS is not GFP and not HIF-α protein. A223 `hypoxia_plus` is a different assay.

## Circular or tautological pieces

### 1. Control persist = 0% is forced

`phenotype_calls(..., exposed)` writes every control cell to **reverted**. E14 \(0\%\) persist / \(100\%\) reverted is therefore **not a result**. It is the design constraint “never-hypoxic cells are not persistent.”

Empirical θ-gate (same \(\theta\le 0.3 / \ge 0.7\), **no** force) on the current E14 fit:

| lineage | empirical persist | empirical partial | empirical reverted |
|---------|-------------------|-------------------|--------------------|
| Tumor (n=152) | 11.2% | 10.5% | 78.3% |
| neutrophil (n=319) | 5.3% | 19.4% | 75.2% |

That is the number that can move: about one in nine E14 tumors still look “persistent” on the HIF panel after CPM. Some of that is noise and residual depth; some is that \(h=1.5\) MAD is not a huge separation when the program barely shifts. **Always report empirical control rates next to forced labels.**

The synthetic check `ctrl_persist_theta < 0.08` used the **forced** gate, so it could not fail. Empirical persist on synthetic never-hypoxic cells is \(\sim 14\%\) — the same order as E14 tumors (\(11\%\)). That is the θ=0.3 / 1.5-MAD gate overlapping the control tail, not a library-specific bug. The test now checks forced persist = 0 and empirical persist \(< 25\%\).

### 2. Two phenotypes that are allowed to disagree

Application hard labels = θ-gate + force-reverted control. Soft GMM is fit on **unscaled ALS \(h\)** with a control boost on the low-\(h\) component. They are not the same estimator.

- Tumors, E15: θ-gate persist \(7.6\%\) vs GMM persist \(1.8\%\). Agreement on exposed \(\approx 80\%\).
- Neutrophils: GMM is **100% reverted in E14 and E15**. The mixture collapsed to one component. Neutrophil “6% persist” exists only on the θ-gate. Do not cite GMM and θ-gate as independent confirmation.

### 3. \(\kappa\) from the cells you later score

\(\hat\kappa_g\) is the mean \(u/s\) in the lowest 20% \(\theta\) of **exposed** cells. Those cells are the persistent-like tail. Lag \(\xi\) is then a residual from that QSS. Mild circularity: persist-like cells help set the “already on” ratio, so they are biased toward \(\xi\approx 0\). That is intentional QSS, but it is not an independent flux measurement in the persist tail.

### 4. Control mean \(\xi = 0\) by construction

Identifiability subtracts the control mean of \(\xi^{\mathrm{raw}}\). “The control has no net HIF flux” is an **assumption**, not a test of net flux. What *can* be tested is spread: tumor control \(p^{\mathrm{away}}\approx 0.16\) is not small. If that is the noise floor, E15 \(0.14\) is not “less reoxygenation.” Do not interpret a drop in mean \(p^{\mathrm{away}}\) as biology when it is within the control floor.

Joint state then folds **reverted × away → reverted**. Away flux in cells already called HIF-off is absorbed. That hides the same control \(p^{\mathrm{away}}\) in the hard table (`frac_hard_flux_control = 0`).

### 5. Balanced WLS vs sample size

E14 tumors n=152 vs E15 n=2327, but WLS gives **equal total weight** to control and exposed. The factor location is half-determined by 152 cells. Empirical 11% E14 persist is partly that: a small, shallower control sets the MAD. Neutrophils are better balanced (319 vs 938).

### 6. Wagner is not an independent experiment

`WagnerCollab/mc38_velocity.h5ad` is the **same** E14S/E15S libraries with different cell typing and Ensembl gene keys. Existing `hypoxia_state` there calls 91% of E15 Tumor cells persistent. Using that column as truth would be circular (agreement with the new θ-gate is 23%).

The fit is a **pipeline replicate**. Empirical control persist matches the placed object (11.5% vs 11.2%). Exposed persist does not (20% vs 7.6%), and neither does depth (UMI ratio 1.01 vs 4.8). Persist after hypoxia is not a stable number across the two objects. “Hypoxic Tumor” on E14 (n=299 before UMI filter) is a cluster name, not exposure.


### 7. A223 DN cannot identify never-hypoxic persist (fit anyway)

DN (DCF−) tumors: **13 (E27) + 12 (E29)**. DN neutrophils: **4 + 3**. The model still ran per chemistry. It does not make persist identifiable. E27 tumors look \(35{-}43\%\) persist on DCF+/DP; E29 tumors look \(6{-}8\%\). Neutrophils on E29 are \(\sim 93\%\) persist because three DN cells set the MAD. DCF is ROS, not HIF. Do not transfer the E14 control location. Do not pool 3′ with 5′. Details: `results/a223.md`.

### 8. Historical object columns

`theta_normoxic` / `hypoxia_kinetics_state` mixed library size into the program and were tumor-only. Neutrophil-run scores must not label tumors. Agreement with the new tumor θ-gate is \(58\%\); old persist was \(30\%\) vs new \(7.6\%\). Treat historical labels as a **confounded baseline**, not as supervision.

## Other flaws (not circular, still real)

- **Glycolysis vs HIF-α.** After depth correction the E15 shift is glycolytic more than CA9/VEGFA/BNIP3. Persist may be a metabolic/HIF-adjacent program, not CA9-high hypoxia.
- **Neutrophil \(h\) vs \(\log L\).** Spearman \(\approx -0.30\) after the same residualization. Neutrophil persist is more depth-entangled than tumor persist.
- **Panel membership is a modeling choice.** HALLMARK_HYPOXIA is not ChIP; CA9 is kept anyway. Changing the panel changes \(h\).
- **Hard flux rarity + \(\sigma_v\) clip** together mean the lag half of the model is barely identified on these libraries. Phenotype (spliced factor) is doing almost all the work.

## What to report going forward

1. Forced phenotype (design) **and** empirical control persist/partial.
2. GMM vs θ-gate; if GMM collapses (neutrophils), say so and do not average them.
3. \(p^{\mathrm{away}}\) in control vs exposed; do not claim reoxygenation if exposed ≤ control.
4. \(\mathrm{corr}(h,\log L)\) and UMI ratio every dataset.
5. Skip the dataset if \(n_{\mathrm{control}} < 40\) (or say the MAD is unidentified).
