# Dataset-agnostic simpleaf 10x pipeline

JSON configs drive two modes:

- `quant` — pooled GEX (± ADT) like the E14S/E28S `run_fast_sbatch.sh` jobs
- `ocm` — same quantification plus Cell Ranger-compatible OCM demux (2 bp overhang at barcode positions 7–8)

```bash
export PATH="/ix1/ylee/kor11/tools/af_tutorial/conda_env/bin:$PATH"
export LD_LIBRARY_PATH="/ix1/ylee/kor11/tools/af_tutorial/conda_env/lib:${LD_LIBRARY_PATH:-}"
export PYTHONNOUSERSITE=1

# E28S OCM (reuses existing simpleaf quants; compares to cellranger multi)
python /ix1/ylee/kor11/tools/af_tutorial/pipeline/run.py run \
  --config /ix1/ylee/kor11/tools/af_tutorial/pipeline/configs/e28s_ocm.json

# E28S pooled (no demux)
python /ix1/ylee/kor11/tools/af_tutorial/pipeline/run.py run \
  --config /ix1/ylee/kor11/tools/af_tutorial/pipeline/configs/e28s_quant.json

# Fresh quantification: set skip_quant to false (or omit it) in the config
python /ix1/ylee/kor11/tools/af_tutorial/pipeline/run.py run --config /path/to/dataset.json
```

Copy `configs/e28s_ocm.json` and point `gex` / `adt` FASTQs, piscem `index`, `ocm.samples`, and optional `cellranger.per_sample_outs` at the new dataset. For a non-OCM 10x run, use `mode: "quant"` and omit `ocm`.

Required environment: simpleaf/piscem/alevin-fry from `af_tutorial/conda_env`, plus `ALEVIN_FRY_HOME` from the config (E28S home already has `e14s-adt-10xv4`).

## Tests that were run

| Test | Result |
|------|--------|
| E28S `quant` | 5003 cells (GEX UMI≥500), GEX+ADT merged |
| E28S `ocm` vs cellranger multi | 100% recall of CR filtered barcodes in every OB sample |
| PBMC 1k toy (2e6 reads), OCM union vs `cellranger count` | 100% recall of 1129 CR cells; cell UMI Pearson 0.999 |

Toy FASTQs: `pipeline/runs/pbmc1k_toy/fastq/`. Cell Ranger outs: `pipeline/runs/pbmc1k_toy/cellranger/pbmc1k_toy_cr/outs/`.
