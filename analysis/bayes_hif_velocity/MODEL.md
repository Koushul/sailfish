# Bayesian 1-D hypoxia velocity

A lineage-restricted model for whether transcription is moving **toward** or **away from** hypoxia. It is not genome-wide RNA velocity. Phenotype \(\theta\) is spliced expression of **direct HIF-α targets**. Direction is residual unspliced lag of those same genes after library size, unspliced capture, and cell-cycle covariates.

Fit separately per lineage (tumor, neutrophil, …). Do not pool lineages on one kinetic scale.

---

## 1. HIF-α target panel

HIF-1α / HIF-2α protein is stabilized in low oxygen and destroyed in high oxygen. The RNA the model uses is **downstream of that switch**: genes with hypoxia-response elements that are induced when HIF-α is high and fall when it is degraded (reoxygenation or PHD-dependent turnover). *HIF1A* / *EPAS1* mRNA and VHL/PHD machinery are **excluded**; protein abundance, not those transcripts, is the oxygen sensor.

The default panel is `hif_targets.tsv` (32 genes). It is the intersection of well-studied direct targets (glycolysis, PDK1, VEGFA, CA9, BNIP3/BNIP3L, NDRG1, DDIT4, P4HA1/2, LOX, EGLN3, CXCR4, SERPINE1, …) with hypoxia-up evidence. Membership in MSigDB **HALLMARK_HYPOXIA** (Liberzon et al. 2015; genes up in low oxygen, not a TF ChIP set) is recorded but is not required: CA9, BNIP3, PKM, EGLN3, and SLC16A3 (MCT4) are kept as core HIF-α targets even when they are missing from that hallmark list.

Sources: Mole et al. 2009 (*J. Biol. Chem.*) HIF-1α/HIF-2α ChIP; Benita et al. 2009 (*NAR*) core HIF-1 response across cell types; Semenza glycolytic HIF-1 targets; HALLMARK_HYPOXIA https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/HALLMARK_HYPOXIA (CC BY 4.0).

Why this panel for velocity:

- persistent hypoxia: high spliced HIF-target program, unspliced at quasi-steady \(\kappa s\) \(\Rightarrow \xi\approx 0\);
- reversion (high \(\mathrm{O}_2\)): program turning off, \(u < \kappa s\) \(\Rightarrow \xi>0\);
- entry into hypoxia: program turning on, \(u > \kappa s\) \(\Rightarrow \xi<0\);
- never-hypoxic control: low program, \(\xi\) pinned to mean 0.

Glycolytic targets are retained even though they correlate with S-phase; cycle is a covariate on unspliced only.

---

## 2. Inputs

For one lineage after QC (minimum spliced UMI):

- integer counts \(U_{ng}, S_{ng}\) for cells \(n=1{\ldots}N\) and HIF-α target genes \(g=1{\ldots}G\);
- sample indicator \(r_n\in\{0,1\}\) (0 = never-hypoxic **control**, 1 = hypoxia-**exposed**);
- observed cell-cycle covariates \(\tilde S_n, \tilde G_n\) (Tirosh S and G2M, centered on the control);
- optional extra covariates in \(x_n\) (the same linear slot).

Spliced library \(L_n\) is the **cell-wide** spliced UMI count (not the sum over the HIF panel). Size-normalizing by the panel sum cancels a coordinated HIF program.

---

## 3. Phenotype \(h,\theta\) (spliced factor)

Size-normalize spliced to the median cell-wide \(L_n\) and control-standardize:

\[
z_{ng}=\frac{\log(1+s^{\mathrm{norm}}_{ng})-\mu_{g,0}}{\sigma_{g,0}}.
\]

Fit a non-negative hypoxia factor (MAP least squares), **not** a global residualization of \(\theta\) on cycle:

\[
z_{ng}=\beta_g h_n+\gamma_g\tilde S_n+\delta_g\tilde G_n+\alpha_g+\varepsilon_{ng},
\qquad \beta_g\ge 0.
\]

\(h_n\) is the shared residual HIF program after gene-specific cycle loadings. Genes with low cycle correlation (CA9, VEGFA, ADM, …) are up-weighted when initializing \(h\). Glycolysis still informs \(h\) through leftover coordinated variation. Sign: high \(h\) = HIF program on. Control cells pin the low-\(h\) end.

\[
\theta_n=\frac{1}{1+\exp\bigl((h_n-h_0)/\tau\bigr)}.
\]

Soft phenotype: a 3-component Gaussian mixture on \(h\) (control boosted on the low-\(h\) component) gives

\[
\bigl(p^{\mathrm{persist}}_n,\,p^{\mathrm{partial}}_n,\,p^{\mathrm{reverted}}_n\bigr).
\]

Hard phenotype is \(\arg\max\) of those probabilities. Cycle is **not** subtracted from \(h\) by a single regression of the score on S/G2M.

\(h\), \(\theta\), and the GMM are frozen during unspliced inference.

---

## 4. Quasi-steady lag (direction)

Full splicing ODEs \(\dot u=\alpha-\beta u\), \(\dot s=\beta u-\gamma s\) are not identifiable on a short, dropout-heavy HIF panel. Use the same lag the point estimator uses.

Let \(\kappa_g\approx\mathbb{E}[u]/\mathbb{E}[s]\) in cells whose spliced phenotype is already hypoxic (low \(\theta\)). If transcription \(\alpha\) changes slowly compared with splicing,

\[
\log \mu^u_{ng}
=\log L_n+\log\rho_{g,r_n}+\log\kappa_g+\log\tilde s_{ng}
+\lambda_g\,\xi_n+c_g^\top x_n,
\]

\[
\tilde s_{ng}=\frac{S_{ng}+\varepsilon}{L_n},\qquad
x_n=(\tilde S_n,\,\tilde G_n,\,r_n).
\]

- \(\rho_{g,1}=1\) (exposed is the capture reference); \(\rho_{g,0}\) is extra unspliced capture in the control library.
- \(\lambda_g\ge 0\) on HIF-α targets that pass unspliced QC; \(\lambda_g=0\) if unspliced is too sparse or control/exposed \(\kappa\) ratios look like capture artifacts.
- \(\xi_n\) is **hypoxia-directed residual lag** after cycle and sample.

Sign: \(v=\mathrm{d}\theta/\mathrm{d}t\). Reversion (\(\theta\) up, HIF transcription down) makes \(u<\kappa s\) after confounders. With \(\lambda_g>0\),

- \(\xi_n>0\): **away from hypoxia** (reverting);
- \(\xi_n<0\): **toward hypoxia** (inducing);
- \(\xi_n\approx 0\): no detectable transcriptional flux (persistent or reverted is then \(\theta\) only).

---

## 5. Hierarchical prior on \(\xi\)

Cycle, sample, and capture are **not** allowed inside \(\xi\). They enter only as gene-specific terms \(c_g^\top x_n\) and \(\rho_{g,r}\). Putting S/G2M into a shared lag would make glycolytic HIF targets look like hypoxia flux. After the hierarchical draw, \(\xi\) is residualized on \((\tilde S,\tilde G)\) so direction is orthogonal to cycle.

\[
\xi^{\mathrm{raw}}_n=\sigma_v\hat\xi_n,
\qquad
\hat\xi_n\sim\mathcal{N}(0,1),
\qquad
\xi^{\mathrm{raw}}\leftarrow \xi^{\mathrm{raw}}-X(X^\top X)^{-1}X^\top\xi^{\mathrm{raw}},
\]

where \(X=[\tilde S,\,\tilde G]\) (columns centered).

Identifiability: subtract the control mean,

\[
\xi_n=\xi^{\mathrm{raw}}_n-\frac{1}{N_0}\sum_{n:r_n=0}\xi^{\mathrm{raw}}_n.
\]

The control has no net HIF flux. Sample-level unspliced shifts go into \(\rho_g\) and \(c_{g,r}\).

Priors (weakly informative):

| parameter | prior |
|-----------|--------|
| \(\log\kappa_g\) | \(\mathcal{N}(\log\hat\kappa_g,\,0.3^2)\), \(\hat\kappa_g\) from low-\(\theta\) exposed cells |
| \(\lambda_g=\mathrm{softplus}(\ell_g)\) | \(\ell_g\sim\mathcal{N}(0,1)\); \(\lambda_g=0\) if gene fails unspliced QC |
| \(c_g\) | \(\mathcal{N}(0,0.3^2 I)\) |
| \(\log\rho_{g,0}\) | \(\mathcal{N}(0,0.5^2)\) |
| \(\log\sigma_v\) | \(\mathcal{N}(-0.2,0.7^2)\) |
| \(\log\phi_g\) (NB concentration) | \(\mathcal{N}(2,1^2)\) |
| \(\mathrm{logit}\,\omega_g\) (zero inflation) | \(\mathcal{N}(-1.2,1^2)\) |
| \(\hat\xi_n\) | \(\mathcal{N}(0,1)\) |

---

## 6. Likelihood

Unspliced is **zero-inflated** negative binomial (dropout is not reversion):

\[
U_{ng}=0 \text{ with extra mass }\omega_g,\quad
\text{else }U_{ng}\sim\mathrm{NegBin}(\mu^u_{ng},\phi_g).
\]

\[
P(U=0)=\omega_g+(1-\omega_g)\,p_0(\mu,\phi),
\qquad
P(U=k>0)=(1-\omega_g)\,\mathrm{NB}(k\mid\mu,\phi).
\]

Spliced counts are conditioned on. Velocity is only in \(U\mid S,x,h\).

---

## 7. Inference

MAP of the ZINB posterior with Adam and analytic gradients. \(\xi\) is control-centered and cycle-orthogonal at every step. A spike-and-slab is **not** used as a MAP prior (it pinned \(\xi\) at 0). After MAP, a Laplace variance \(\mathrm{Var}(\xi_n)\) is combined with a three-component slab

\[
\xi\sim \pi_0\mathcal{N}(0,\tau_0^2)+\pi_+\mathcal{N}(+\mu,\tau_1^2)+\pi_-\mathcal{N}(-\mu,\tau_1^2)
\]

with \(\mu=0.7\), \(\tau_0=0.28\), \(\tau_1=0.55\). Mixing weights \(\pi\) are empirical-Bayes from the Laplace posterior. Component posterior:

\[
P(k\mid U)\propto \pi_k\,\mathcal{N}\bigl(\xi^{\mathrm{MAP}}_n; m_k,\,\mathrm{Var}(\xi_n)+\tau_k^2\bigr).
\]

\(p^{\mathrm{away}}=P(+)\), \(p^{\mathrm{toward}}=P(-)\), \(p^{\mathrm{none}}=P(0)\).

Joint state probability is the product of GMM phenotype and flux components (reverted \(\times\) away is folded into reverted). Hard call is \(\arg\max\). Use probabilities, not only hard labels.

---

## 8. Phenotype, flux, and joint states

| state | phenotype component | flux component |
|-------|---------------------|----------------|
| persistent | persist | none |
| persistent_exiting | persist | away |
| persistent_deepening | persist | toward |
| partial | partial | none |
| transitioning_out | partial | away |
| transitioning_in | partial | toward |
| reverted | reverted | none or away |
| reverted_entering | reverted | toward |

Partial vs transitioning is the **flux posterior**, not a \(\theta\in(0.3,0.7)\) cut.

---

## 9. Fitting recipe

1. QC cells in one lineage by spliced UMI.
2. Score Tirosh S/G2M; center on control.
3. Restrict to the HIF-α panel; fit the spliced hypoxia factor and GMM phenotype.
4. Estimate \(\hat\kappa_g\) on exposed cells with low \(\theta\).
5. MAP ZINB lag (Adam); Laplace + spike-slab for flux probabilities. Clip control-standardized spliced \(z\) to \([-6,6]\) and \(\sigma_v\le e^{0.5}\) so a few cells cannot explode \(h\) or \(\xi\).
6. Pin control mean \(\xi\) to 0; orthogonalize \(\xi\) to cycle.
7. Write posterior phenotype, flux, joint state, and \(\mathbb{E}[\xi]\).
8. Write per-cell \(\theta\), \(\mathbb{E}[\xi]\), \(p^{\mathrm{away}}\), \(p^{\mathrm{toward}}\), joint state, and gene parameters \(\kappa_g,\lambda_g,c_g,\rho_g\).

---

## 10. Checks (synthetic and real)

- Control: fraction of confident toward/away calls near the nominal 5% tail.
- Spearman \((\mathbb{E}[\xi],\tilde S)\) should be weak if cycle was absorbed in \(c_g\).
- On synthetic HIF-α panels: GMM should recover partial vs persist/reverted; flux AUROC should separate transitioning_out from partial; control false-flux low; scrambled \(U\) should destroy lag AUROC.
- Weak lag (true \(\xi\) scaled \(\sim 0.4\)) is below detection — report probabilities, do not force hard transition calls.
- Strong cycle should not delete the GMM partial component if CA9/VEGFA-like anchors are in the panel.

Sweep (`benchmark_synthetic.py`, 3 seeds, overlapping states, discrete cycle, decoys, dropout, mixed lag): realistic Spearman \(\hat\xi\sim 0.24{-}0.32\); transitioning_out vs partial AUROC \(\sim 0.93\); GMM partial recall \(\sim 0.9\); inducing AUROC \(\sim 0.69{-}0.77\); control false flux \(\sim 1\%\). Hard `transitioning_out` recall stays low (\(\sim 0.1\)); use \(p^{\mathrm{away}}\). Strong cycle: partial phenotype holds (\(\sim 0.8{-}0.93\)), lag Spearman collapses. Weak lag: not recovered. Scrambled \(U\): lag AUROC \(\sim 0.5\). Table: `results/synthetic_benchmark.tsv`.

---

## 11. Application: E14 vs E15 tumors and neutrophils

E14S = never-hypoxic control, E15S = hypoxia-exposed. Fits are lineage-specific (`fit_e14e15.py`; tables in `results/e14e15_*.tsv`). Tirosh S/G2M is recomputed from spliced counts (the stored `cycle_s` column is tumor-only). Cell-wide spliced UMI floors: tumors \(\ge 5000\) (606 / 3085 dropped), neutrophils \(\ge 2000\) (219 / 1476 dropped). All 32 panel genes are present as mouse symbols; unspliced-sparse genes keep \(\lambda_g=0\) (23 velocity genes in tumors, 14 in neutrophils). *Car9* loads the tumor factor and lag; it is silent in neutrophils.

On this object the 3-component GMM on \(h\) collapses the partial component. Reported phenotype uses the same \(\theta\le 0.3\) / \(\ge 0.7\) gate as the earlier tumor-only analysis. E15 tumors (n=2327): persist / partial / reverted \(\approx 27.5\% / 16.4\% / 56.2\%\), vs historical \(\theta\) gates \(30.3\% / 23.0\% / 46.7\%\) (71% cell-level agreement). E14 tumors (n=152) are not a clean never-hypoxic pile: \(\approx 21\%\) still gate as persistent, so a glycolytic HIF-like program is already on in some controls.

Lag does **not** add a clear E15-specific reoxygenation wave. Tumor \(\mathbb{E}[\xi]\) is nearly uncorrelated with cycle (Spearman \(-0.03\)). Mean \(p^{\mathrm{away}}\) is similar in E14 (\(0.15\)) and E15 (\(0.13\)); hard flux is \(2.6\%\) of E14 and \(1.1\%\) of E15. Reverted E15 tumors have modestly higher \(p^{\mathrm{away}}\) (\(0.17\)) than persistent (\(0.09\)). \(\sigma_v\) sits at the upper clip (\(\approx 1.65\)).

Neutrophils have a weaker panel (no *Car9* unspliced; glycolysis + *Vegfa* / *Ndrg1* / *P4ha1* / *Egln3*). E15 θ-gates \(\approx 26\% / 21\% / 54\%\) persist / partial / reverted, but that is **not** tumor kinetics: \(h\) is weakly shifted, control \(p^{\mathrm{away}}\) is low (\(0.02\)), and hard flux is \(0\). Do not use neutrophil states as tumor labels.

Figures: `results/e14e15_tumor_overview.png`, `results/e14e15_neutrophil_overview.png`. Narrative: `results/e14e15.md`.

---

## 12. What this does not do

- Infer a shared latent time across the transcriptome.
- Invent flux for genes with no unspliced counts (those posteriors stay wide).
- Treat HIF-1α protein destruction as RNA velocity; if mRNA is already down, \(\xi\approx 0\) and class is \(\theta\) only.
