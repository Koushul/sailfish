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

## 3. Phenotype \(\theta\) (spliced only)

Size-normalize spliced to the median \(L_n\). For each gene,

\[
z_{ng}=\frac{\log(1+s^{\mathrm{norm}}_{ng})-\mu_{g,0}}{\sigma_{g,0}},
\]

where \(\mu_{g,0},\sigma_{g,0}\) are the control mean and sd. With fixed non-negative weights \(w_g\) (normalized Cohen’s \(d\) of exposed vs control among HIF-α targets, clipped at 0),

\[
h_n=\sum_g w_g z_{ng},\qquad
\theta_n=\frac{1}{1+\exp\bigl((h_n-h_0)/\tau\bigr)}.
\]

Anchors: \(h_0\) is the midpoint of the control median of \(h\) and the exposed upper-quartile median of \(h\); \(\tau=|h_{\mathrm{hyp}}-h_{\mathrm{ctrl}}|/6\). High \(\theta\) = low HIF (reverted / never hypoxic). Low \(\theta\) = high HIF (persistent).

**Do not residualize cycle out of \(\theta\).** Glycolytic HIF targets are collinear with growth; subtracting S/G2M deletes the hypoxia axis. Cycle enters only the unspliced lag (below).

\(\theta\) and \(w_g\) are **frozen** during velocity inference.

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
| \(\log\sigma_v\) | \(\mathcal{N}(-2,0.5^2)\) |
| \(\log\phi_g\) (NB concentration) | \(\mathcal{N}(2,1^2)\) |
| \(\hat\xi_n\) | \(\mathcal{N}(0,1)\) |

---

## 6. Likelihood

\[
U_{ng}\sim\mathrm{NegBin}\bigl(\mu^u_{ng},\,\phi_g\bigr),
\qquad
\mathrm{Var}(U)=\mu+\mu^2/\phi.
\]

Spliced counts are conditioned on (they already defined \(\theta\)). The velocity information is only in \(U\mid S,x,\theta\).

---

## 7. Inference

Let \(z\) collect all parameters (\(\kappa,\lambda,c,\rho,a,\sigma_v,\phi,\hat\xi\)). Maximize the joint log posterior

\[
\mathcal{L}(z)=\sum_{n,g}\log\mathrm{NB}\bigl(U_{ng}\mid\mu^u_{ng}(z),\phi_g\bigr)+\log p(z)
\]

with Adam (analytic gradients; no finite differences). For \(\mathrm{NB}(k\mid\mu,\phi)\) with \(\mathrm{Var}=\mu+\mu^2/\phi\),

\[
\frac{\partial\log p}{\partial\log\mu}
=(k-\mu)\,\frac{\phi}{\phi+\mu}.
\]

The remaining chain rule through \(\log\mu^u\) is in the implementation (`model.py`). Identifiability \(\xi\leftarrow\xi-\mathrm{mean}_{\mathrm{control}}(\xi)\) is applied at every Adam step, and the control-mean subtraction is included in the \(\xi\) gradient.

After the MAP \(\hat z\), cell-level uncertainty is a Laplace approximation treating other parameters as fixed:

\[
\mathrm{Var}(\xi_n)
\approx
\sigma_v^2\Big/\Big(1+\sum_g \lambda_g^2\,\frac{\phi_g\,\mu_{ng}}{\phi_g+\mu_{ng}}\Big),
\]
\[
p^{\mathrm{away}}_n=\Pr(\xi_n>0\mid U)
\approx 1-\Phi\bigl(0;\,\xi_n^{\mathrm{MAP}},\,\sqrt{\mathrm{Var}(\xi_n)}\bigr),
\qquad
p^{\mathrm{toward}}_n=1-p^{\mathrm{away}}_n.
\]

If all \(\lambda_g=0\) (no usable unspliced), set \(\xi_n=0\) and return \(\theta\)-only states.

---

## 8. Phenotype, flux, and joint states

Two axes are scored separately, then combined. \(\theta\) is the spliced HIF-α program (where the cell sits). \(\xi\) is residual unspliced lag (whether it is moving).

Phenotype from \(\theta\):

| phenotype | \(\theta\) | meaning |
|-----------|------------|---------|
| persistent | \(\le 0.3\) | HIF program on |
| partial | \((0.3,0.7)\) | intermediate / partially reverted program |
| reverted | \(\ge 0.7\) | HIF program off |

Flux from control-calibrated tails. Let \(q_{95}^{0}(p^{\mathrm{away}})\) be the 95th percentile of \(p^{\mathrm{away}}\) among control cells (empirical null). Analogous for toward.

| flux | rule |
|------|------|
| away | \(p^{\mathrm{away}}\ge \max\bigl(q_{95}^{0}(p^{\mathrm{away}}),\,0.8\bigr)\) |
| toward | \(p^{\mathrm{toward}}\ge \max\bigl(q_{95}^{0}(p^{\mathrm{toward}}),\,0.8\bigr)\) |
| none | neither |

Joint call (flux only changes the label when it is biologically possible at that \(\theta\)):

| state | phenotype | flux | meaning |
|-------|-----------|------|---------|
| persistent | persistent | none | stable hypoxia |
| persistent_exiting | persistent | away | still HIF-high spliced, transcription already shutting off |
| persistent_deepening | persistent | toward | already HIF-high, still inducing |
| partial | partial | none | **partially reverted**, no detectable flux |
| transitioning_out | partial | away | leaving hypoxia through the intermediate |
| transitioning_in | partial | toward | entering hypoxia through the intermediate |
| reverted | reverted | none or away | stable reversion (away lag ignored once spliced is already off) |
| reverted_entering | reverted | toward | spliced already reverted, transcription turning the program back on |

Partial vs transitioning is the lag test: the same mid-\(\theta\) bin is **stuck intermediate** if \(\xi\approx 0\), and **in transit** if the control-calibrated tail fires.

---

## 9. Fitting recipe

1. QC cells in one lineage by spliced UMI.
2. Score Tirosh S/G2M; center on control.
3. Restrict to the HIF-α target panel; build \(\theta\) from those spliced counts (frozen weights).
4. Estimate \(\hat\kappa_g\) on exposed cells with \(\theta\) in the lowest quintile. Drop genes with unspliced detection below a floor, or with control/exposed \(\kappa\) ratio above a cap (capture, not kinetics).
5. Run MAP (Adam, analytic gradients) on \(U\mid S,x,\theta\); Laplace for \(p^{\mathrm{away}}\).
6. Pin control mean \(\xi\) to 0 at every optimization step.
7. Calibrate flux tails on the control posterior; write phenotype, flux, and joint state.
8. Write per-cell \(\theta\), \(\mathbb{E}[\xi]\), \(p^{\mathrm{away}}\), \(p^{\mathrm{toward}}\), joint state, and gene parameters \(\kappa_g,\lambda_g,c_g,\rho_g\).

---

## 10. Checks (synthetic and real)

- Control: fraction of confident toward/away calls near the nominal 5% tail.
- Spearman \((\mathbb{E}[\xi],\tilde S)\) should be weak if cycle was absorbed in \(c_g\).
- On synthetic HIF-α panels with known states: recover rank correlation of \(\mathbb{E}[\xi]\) vs truth; AUROC of \(p^{\mathrm{away}}\) for cells leaving hypoxia and of \(p^{\mathrm{toward}}\) for cells entering; spliced \(\theta\) should recover persistent / **partial** / reverted when \(\xi\approx 0\); mid-\(\theta\) cells with true lag should be called transitioning, not partial.
- If unspliced is permuted (independent of \(\theta\) by construction), lag recovery should collapse while partial vs persist/reverted from \(\theta\) remains.
- If unspliced is permuted (independent of \(\theta\) by construction), lag recovery should collapse.
- Cycle-only decoy genes added to the panel should not make \(\xi\) track S-phase.

Sweep (`benchmark_synthetic.py`, 3 seeds): default panel Spearman \(\hat\xi\) vs truth \(\approx 0.72\); AUROC of \(p^{\mathrm{away}}\) for transitioning_out vs partial \(\approx 1\). Spliced \(\theta\) recovers persist/reverted near 1 and partial \(\approx 0.4{-}0.6\). Joint hard calls for partial vs transitioning_out are noisier than the AUROC (mid-\(\theta\) is a thin band). Weak lag drops transition recall; scrambled unspliced drops lag AUROC to \(\sim 0.5\) while \(\theta\) partial remains. Table: `results/synthetic_benchmark.tsv`.

---

## 11. What this does not do

- Infer a shared latent time across the transcriptome.
- Invent flux for genes with no unspliced counts (those posteriors stay wide).
- Treat HIF-1α protein destruction as RNA velocity; if mRNA is already down, \(\xi\approx 0\) and class is \(\theta\) only.
