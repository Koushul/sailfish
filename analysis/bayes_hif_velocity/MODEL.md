# Bayesian 1-D hypoxia velocity

A lineage-restricted model for whether transcription is moving **toward** or **away from** hypoxia. It is not genome-wide RNA velocity. Phenotype \(\theta\) is a spliced HIF-down axis. Direction is residual unspliced lag of those same genes after library size, unspliced capture, and cell-cycle covariates.

Fit separately per lineage (tumor, neutrophil, …). Do not pool lineages on one kinetic scale.

---

## 1. Inputs

For one lineage after QC (minimum spliced UMI):

- integer counts \(U_{ng}, S_{ng}\) for cells \(n=1{\ldots}N\) and HIF-down genes \(g=1{\ldots}G\);
- sample indicator \(r_n\in\{0,1\}\) (0 = never-hypoxic **control**, 1 = hypoxia-**exposed**);
- observed cell-cycle covariates \(\tilde S_n, \tilde G_n\) (Tirosh S and G2M, centered on the control);
- optional extra covariates in \(x_n\) (the same linear slot).

Spliced library \(L_n=\sum_g S_{ng}\) is treated as observed. Unspliced is **not** given its own library size (that confuses capture with kinetics).

---

## 2. Phenotype \(\theta\) (spliced only)

Size-normalize spliced to the median \(L_n\). For each gene,

\[
z_{ng}=\frac{\log(1+s^{\mathrm{norm}}_{ng})-\mu_{g,0}}{\sigma_{g,0}},
\]

where \(\mu_{g,0},\sigma_{g,0}\) are the control mean and sd. With fixed non-negative weights \(w_g\) (normalized Cohen’s \(d\) of exposed vs control, clipped at 0),

\[
h_n=\sum_g w_g z_{ng},\qquad
\theta_n=\frac{1}{1+\exp\bigl((h_n-h_0)/\tau\bigr)}.
\]

Anchors: \(h_0\) is the midpoint of the control median of \(h\) and the exposed upper-quartile median of \(h\); \(\tau=|h_{\mathrm{hyp}}-h_{\mathrm{ctrl}}|/6\). High \(\theta\) = low HIF (reverted / never hypoxic). Low \(\theta\) = high HIF (persistent).

**Do not residualize cycle out of \(\theta\).** Glycolytic HIF targets are collinear with growth; subtracting S/G2M deletes the hypoxia axis. Cycle enters only the unspliced lag (below).

\(\theta\) and \(w_g\) are **frozen** during velocity inference.

---

## 3. Quasi-steady lag (direction)

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
- \(\lambda_g\ge 0\) on HIF-down genes that pass unspliced QC; \(\lambda_g=0\) if unspliced is too sparse or control/exposed \(\kappa\) ratios look like capture artifacts.
- \(\xi_n\) is **hypoxia-directed residual lag** after cycle and sample.

Sign: \(v=\mathrm{d}\theta/\mathrm{d}t\). Reversion (\(\theta\) up, HIF transcription down) makes \(u<\kappa s\) after confounders. With \(\lambda_g>0\),

- \(\xi_n>0\): **away from hypoxia** (reverting);
- \(\xi_n<0\): **toward hypoxia** (inducing);
- \(\xi_n\approx 0\): no detectable transcriptional flux (persistent or reverted is then \(\theta\) only).

---

## 4. Hierarchical prior on \(\xi\)

\[
\xi^{\mathrm{raw}}_n
= a_\theta(\theta_n-\bar\theta)
+ a_S\tilde S_n
+ a_{G2M}\tilde G_n
+ a_r r_n
+ \sigma_v\hat\xi_n,
\qquad
\hat\xi_n\sim\mathcal{N}(0,1).
\]

Identifiability: subtract the control mean,

\[
\xi_n=\xi^{\mathrm{raw}}_n-\frac{1}{N_0}\sum_{n:r_n=0}\xi^{\mathrm{raw}}_n.
\]

The control has no net HIF flux. Sample-level shifts go into \(a_r\) and \(\rho_g\), not into every cell’s direction.

Priors (weakly informative):

| parameter | prior |
|-----------|--------|
| \(\log\kappa_g\) | \(\mathcal{N}(\log\hat\kappa_g,\,0.3^2)\), \(\hat\kappa_g\) from low-\(\theta\) exposed cells |
| \(\lambda_g=\mathrm{softplus}(\ell_g)\) | \(\ell_g\sim\mathcal{N}(0,1)\); \(\lambda_g=0\) if gene fails unspliced QC |
| \(c_g\) | \(\mathcal{N}(0,0.3^2 I)\) |
| \(\log\rho_{g,0}\) | \(\mathcal{N}(0,0.5^2)\) |
| \(a_\theta,a_S,a_{G2M},a_r\) | \(\mathcal{N}(0,0.5^2)\) |
| \(\log\sigma_v\) | \(\mathcal{N}(-2,0.5^2)\) |
| \(\log\phi_g\) (NB concentration) | \(\mathcal{N}(2,1^2)\) |
| \(\hat\xi_n\) | \(\mathcal{N}(0,1)\) |

---

## 5. Likelihood

\[
U_{ng}\sim\mathrm{NegBin}\bigl(\mu^u_{ng},\,\phi_g\bigr),
\qquad
\mathrm{Var}(U)=\mu+\mu^2/\phi.
\]

Spliced counts are conditioned on (they already defined \(\theta\)). The velocity information is only in \(U\mid S,x,\theta\).

---

## 6. Inference

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

## 7. Direction calls

Let \(q_{95}^{0}(p^{\mathrm{away}})\) be the 95th percentile of \(p^{\mathrm{away}}\) among **control** cells (empirical null). Analogous for toward.

| call | \(\theta\) | lag |
|------|------------|-----|
| persistent | \(\mathbb{E}[\theta]\le 0.3\) | neither tail exceeds the control 95th percentile |
| reverted | \(\mathbb{E}[\theta]\ge 0.7\) | same |
| partial | \(0.3<\mathbb{E}[\theta]<0.7\) | same |
| reverting (away) | \(\mathbb{E}[\theta]\le 0.7\) | \(p^{\mathrm{away}}\ge q_{95}^{0}(p^{\mathrm{away}})\) |
| inducing (toward) | \(\mathbb{E}[\theta]\ge 0.3\) | \(p^{\mathrm{toward}}\ge q_{95}^{0}(p^{\mathrm{toward}})\) |

Phenotype still blocks the wrong direction of travel (high \(\theta\) cannot be called reverting).

---

## 8. Fitting recipe

1. QC cells in one lineage by spliced UMI.
2. Score Tirosh S/G2M; center on control.
3. Build \(\theta\) from HIF-down spliced (frozen weights).
4. Estimate \(\hat\kappa_g\) on exposed cells with \(\theta\) in the lowest quintile. Drop genes with unspliced detection below a floor, or with control/exposed \(\kappa\) ratio above a cap (capture, not kinetics).
5. Run MAP (Adam, analytic gradients) on \(U\mid S,x,\theta\); Laplace for \(p^{\mathrm{away}}\).
6. Pin control mean \(\xi\) to 0 at every optimization step.
7. Calibrate direction calls on the control posterior tail.
8. Write per-cell \(\theta\), \(\mathbb{E}[\xi]\), \(p^{\mathrm{away}}\), \(p^{\mathrm{toward}}\), and gene parameters \(\kappa_g,\lambda_g,c_g,\rho_g\).

---

## 9. Checks (synthetic and real)

- Control: fraction of confident toward/away calls near the nominal 5% tail.
- Spearman \((\mathbb{E}[\xi],\tilde S)\) should be weak if cycle was absorbed in \(c_g\).
- On synthetic data with known \(\xi^{\mathrm{true}}\): recover rank correlation of \(\mathbb{E}[\xi]\) vs truth; higher \(p^{\mathrm{away}}\) in true reverting cells than in control.
- If unspliced is independent of \(\theta\) by construction, posterior \(\xi\) should collapse toward 0.

---

## 10. What this does not do

- Infer a shared latent time across the transcriptome.
- Invent flux for genes with no unspliced counts (those posteriors stay wide).
- Treat HIF-1α protein destruction as RNA velocity; if mRNA is already down, \(\xi\approx 0\) and class is \(\theta\) only.
