# Bayesian 1-D hypoxia velocity

See `MODEL.md` for the generative model, priors, MAP + Laplace fit, and direction calls.

```bash
/ix1/ylee/kor11/tools/af_tutorial/conda_env/bin/python analysis/bayes_hif_velocity/test_synthetic.py
```

Synthetic data: control vs hypoxia-exposed cells, Tirosh-like cycle covariates, gene-specific unspliced capture, two genes with \(\lambda=0\), and a minority of exposed cells with true toward/away lag.
