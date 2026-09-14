# wagner application

Same E14S (never-hypoxic) vs E15S (hypoxia-exposed) MC38 experiment as the placed E14SE15S object, in the Wagner velocity AnnData. Gene keys are Ensembl `_index` plus `var.gene_name`. Tumor lineage is the union of Tumor / Hypoxic Tumor / Tumor Proliferating cluster names (cluster names are not hypoxia labels). Existing `hypoxia_state` is a previous analysis, not ground truth.

Control sample `E14S`; exposed `E15S`. Forced control labels are design constraints. Empirical θ-gates on control are the diagnostic for circular persist calls.

## QC

| lineage    |   min_spliced_umi |   n_in_lineage |   n_dropped_umi |   n_fit |   n_control |   n_exposed |   n_genes |   n_detected_genes |   n_velocity_genes |   median_spliced_umi | status   |
|:-----------|------------------:|---------------:|----------------:|--------:|------------:|------------:|----------:|-------------------:|-------------------:|---------------------:|:---------|
| Tumor      |              5000 |           2400 |             359 |    2041 |         192 |        1849 |        32 |                 24 |                 23 |              8728.98 | fit      |
| neutrophil |              2000 |           1725 |               0 |    1725 |         495 |        1230 |        32 |                 13 |                 13 |              8268.36 | fit      |

## Sample-level fractions

| lineage    | sample   |   n_cells |   frac_pheno_persistent |   frac_pheno_partial |   frac_pheno_reverted |   frac_pheno_empirical_persistent |   mean_p_away |   mean_p_toward |   frac_hard_flux |
|:-----------|:---------|----------:|------------------------:|---------------------:|----------------------:|----------------------------------:|--------------:|----------------:|-----------------:|
| Tumor      | E14S     |       192 |               0         |             0        |              1        |                         0.114583  |     0.0826006 |       0.0271201 |       0          |
| Tumor      | E15S     |      1849 |               0.203894  |             0.345592 |              0.450514 |                         0.203894  |     0.0201051 |       0.0212993 |       0          |
| neutrophil | E14S     |       495 |               0         |             0        |              1        |                         0.0464646 |     0.0606576 |       0.0288338 |       0          |
| neutrophil | E15S     |      1230 |               0.0341463 |             0.290244 |              0.67561  |                         0.0341463 |     0.0422872 |       0.0342235 |       0.00325203 |

## Skipped

None.

## Diagnostics

### Tumor

- `n_cells`: 2041
- `n_control`: 192
- `n_exposed`: 1849
- `median_spliced_umi_control`: 8665.0901
- `median_spliced_umi_exposed`: 8736.4164
- `umi_ratio_exposed_over_control`: 1.0082
- `spearman_xi_cycle_s`: -0.0270
- `spearman_h_cycle_s`: -0.1406
- `spearman_h_logL`: 0.0336
- `mean_abs_xi_control`: 0.4522
- `mean_p_away_control`: 0.0826
- `mean_p_toward_control`: 0.0271
- `frac_hard_flux_control`: 0.0000
- `control_frac_persistent_forced`: 0.0000
- `control_frac_persistent_empirical`: 0.1146
- `control_frac_partial_empirical`: 0.1354
- `control_frac_reverted_empirical`: 0.7500
- `control_gmm_frac_persistent`: 0.0521
- `exposed_gmm_frac_persistent`: 0.0973
- `mean_p_away_exposed`: 0.0201
- `mean_p_toward_exposed`: 0.0213
- `frac_hard_flux_exposed`: 0.0000
- `sigma_v`: 1.6487
- `pi_none`: 0.9143
- `pi_away`: 0.0446
- `pi_toward`: 0.0411
- `mean_lambda_used`: 1.8038
- `exposed_frac_persistent`: 0.2039
- `exposed_n_persistent`: 377
- `exposed_mean_p_away_persistent`: 0.0190
- `exposed_mean_p_toward_persistent`: 0.0181
- `exposed_frac_partial`: 0.3456
- `exposed_n_partial`: 639
- `exposed_mean_p_away_partial`: 0.0198
- `exposed_mean_p_toward_partial`: 0.0211
- `exposed_frac_reverted`: 0.4505
- `exposed_n_reverted`: 833
- `exposed_mean_p_away_reverted`: 0.0208
- `exposed_mean_p_toward_reverted`: 0.0229
- `hist_n_compared`: 1849
- `hist_agree_pheno`: 0.2347
- `hist_frac_persistent`: 0.9064
- `new_frac_persistent_on_hist_cells`: 0.2039
- `hist_frac_partial`: 0.0011
- `new_frac_partial_on_hist_cells`: 0.3456
- `hist_frac_reverted`: 0.0925
- `new_frac_reverted_on_hist_cells`: 0.4505
- `n_velocity_genes`: 23
- `loss_final`: 16716.0996

### neutrophil

- `n_cells`: 1725
- `n_control`: 495
- `n_exposed`: 1230
- `median_spliced_umi_control`: 8388.4100
- `median_spliced_umi_exposed`: 8206.3500
- `umi_ratio_exposed_over_control`: 0.9783
- `spearman_xi_cycle_s`: -0.0136
- `spearman_h_cycle_s`: -0.1339
- `spearman_h_logL`: 0.1336
- `mean_abs_xi_control`: 0.3966
- `mean_p_away_control`: 0.0607
- `mean_p_toward_control`: 0.0288
- `frac_hard_flux_control`: 0.0000
- `control_frac_persistent_forced`: 0.0000
- `control_frac_persistent_empirical`: 0.0465
- `control_frac_partial_empirical`: 0.2162
- `control_frac_reverted_empirical`: 0.7374
- `control_gmm_frac_persistent`: 0.0000
- `exposed_gmm_frac_persistent`: 0.0000
- `mean_p_away_exposed`: 0.0423
- `mean_p_toward_exposed`: 0.0342
- `frac_hard_flux_exposed`: 0.0033
- `sigma_v`: 1.6487
- `pi_none`: 0.8868
- `pi_away`: 0.0629
- `pi_toward`: 0.0503
- `mean_lambda_used`: 3.2826
- `exposed_frac_persistent`: 0.0341
- `exposed_n_persistent`: 42
- `exposed_mean_p_away_persistent`: 0.0400
- `exposed_mean_p_toward_persistent`: 0.0450
- `exposed_frac_partial`: 0.2902
- `exposed_n_partial`: 357
- `exposed_mean_p_away_partial`: 0.0492
- `exposed_mean_p_toward_partial`: 0.0403
- `exposed_frac_reverted`: 0.6756
- `exposed_n_reverted`: 831
- `exposed_mean_p_away_reverted`: 0.0394
- `exposed_mean_p_toward_reverted`: 0.0311
- `hist_n_compared`: 1230
- `hist_agree_pheno`: 0.4187
- `hist_frac_persistent`: 0.2146
- `new_frac_persistent_on_hist_cells`: 0.0341
- `hist_frac_partial`: 0.1561
- `new_frac_partial_on_hist_cells`: 0.2902
- `hist_frac_reverted`: 0.6293
- `new_frac_reverted_on_hist_cells`: 0.6756
- `n_velocity_genes`: 13
- `loss_final`: 19415.7923

## Versus the placed E14SE15S object

This is the same exposure, not a new experiment. *Lox* is missing from `gene_name`. Tumor cells are the union of Tumor / Hypoxic Tumor / Tumor Proliferating (cluster names, not hypoxia labels).

| check | placed E14SE15S | Wagner |
|-------|-----------------|--------|
| tumor UMI ratio E15/E14 | 4.77 | 1.01 |
| tumor empirical E14 persist | 11.2% | 11.5% |
| tumor E15 persist (forced gates on exposed) | 7.6% | 20.4% |
| tumor \(p^{\mathrm{away}}\) E14 vs E15 | 0.16 vs 0.14 | 0.083 vs 0.020 |
| neutrophil E15 persist | 6.1% | 3.4% |
| neutrophil GMM persist | 0% both samples | 0% both samples |

Empirical never-hypoxic persist matching at \(\sim 11\%\) is the stable diagnostic: the \(\theta=0.3\) gate always clips a control tail. E15 persist does **not** match (7.6% vs 20%), and library-size ratio does not match. Wagner spliced medians are \(\sim 8.7\mathrm{k}\) in both samples; the placed object is \(\sim 14\mathrm{k}\) vs \(\sim 68\mathrm{k}\). Different quantification / filtering, not two independent hypoxia experiments.

Wagner `hypoxia_state` called 91% of E15 tumor cells persistent. Agreement with the new θ-gate is 23%. Do not treat that column as truth.

Lag still does not show an E15 reoxygenation wave (exposed \(p^{\mathrm{away}}\) ≤ control). \(\sigma_v\) is again at the clip.

