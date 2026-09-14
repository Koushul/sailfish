# Bayesian 1-D hypoxia velocity

See `MODEL.md` for the HIF-α panel, spliced factor + GMM phenotype, ZINB lag, and Laplace spike-slab flux posteriors. Gene list: `hif_targets.tsv`.

```bash
/ix1/ylee/kor11/tools/af_tutorial/conda_env/bin/python analysis/bayes_hif_velocity/test_synthetic.py
/ix1/ylee/kor11/tools/af_tutorial/conda_env/bin/python analysis/bayes_hif_velocity/benchmark_synthetic.py
/ix1/ylee/kor11/tools/af_tutorial/conda_env/bin/python analysis/bayes_hif_velocity/fit_e14e15.py
/ix1/ylee/kor11/tools/af_tutorial/conda_env/bin/python analysis/bayes_hif_velocity/fit_dataset.py --config analysis/bayes_hif_velocity/datasets/wagner.json
/ix1/ylee/kor11/tools/af_tutorial/conda_env/bin/python analysis/bayes_hif_velocity/report_a223.py
```

Circularity audit: `results/critique.md`. A223 is not fitted (DN n too small).

Synthetic cells have overlapping HIF programs, discrete G1/S/G2M, cycle-only decoys, silent targets, mixed weak/strong lag, capture shift, and unspliced dropout. Trust posterior probabilities more than hard `transitioning_out` labels.

`fit_e14e15.py` reads the placed E14SE15S h5ad **read-only** and writes tumor/neutrophil tables under `results/` (`e14e15.md`). Do not write back into the h5ad. Fit the two lineages separately.

Synthetic data are control vs hypoxia-exposed cells whose spliced/unspliced counts are generated from the curated HIF-α panel (glycolysis, VEGFA, CA9, BNIP3, …), with gene-specific capture, dropout, and cycle coupling. True states include persistent, partially reverted (mid \(\theta\), no flux), reverted, transitioning in/out, persistent-exiting, and reverted-entering. The benchmark varies panel size, dropout, cycle strength, decoy genes, lag scale, and scrambled unspliced counts.
