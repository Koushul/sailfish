# Bayesian 1-D hypoxia velocity

See `MODEL.md` for the HIF-α target panel, generative model, MAP + Laplace fit, and direction calls. The gene list is `hif_targets.tsv`.

```bash
/ix1/ylee/kor11/tools/af_tutorial/conda_env/bin/python analysis/bayes_hif_velocity/test_synthetic.py
/ix1/ylee/kor11/tools/af_tutorial/conda_env/bin/python analysis/bayes_hif_velocity/benchmark_synthetic.py
```

Synthetic data are control vs hypoxia-exposed cells whose spliced/unspliced counts are generated from the curated HIF-α panel (glycolysis, VEGFA, CA9, BNIP3, …), with gene-specific capture, dropout, and cycle coupling. The benchmark varies panel size, dropout, cycle strength, decoy genes, lag scale, and scrambled unspliced counts.
