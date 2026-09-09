# Sailfish

simpleaf / alevin-fry 10x GEX (± ADT, ± VDJ TCR/BCR), optional Cell Ranger–compatible OCM demux.

Needs `simpleaf`, `piscem`, `alevin-fry` on `PATH` (or `export SAILFISH_CONDA=/path/to/conda_env`). Python: `anndata`, `numpy`, `pandas`, `scipy`. Splici indexing also needs `roers` (`cargo install roers`). Sibling `conda_env` / `cargo_tools` are auto-detected. `ALEVIN_FRY_HOME` defaults to `{output}/af_home`.

```bash
python run.py my_sample.json
python run.py my_sample.json -o /tmp/out --skip-quant
sbatch run_sbatch.sh configs/e28s_ocm.json
```

`mode: "quant"` is pooled; `"ocm"` splits by the 2 bp overhang at barcode positions 7–8 (GEM-X: **GT=OB1, CA=OB2, TC=OB3, AG=OB4**).

## Indices

Config `*.index` must be the folder `simpleaf index --output DIR` writes as `DIR/index` (`piscem_idx*`, `t2g_3col.tsv`). Use a local `--work-dir`.

GEX splici (genome+GTF → `spliced` / `unspliced` / `ambiguous`). `--rlen` ≈ cDNA length (GEM-X v4 ~90; 10xv3 often 91). Default `--ref-type` is `spliced+intronic`.

```bash
export ALEVIN_FRY_HOME="${ALEVIN_FRY_HOME:-./af_home}"
simpleaf set-paths
simpleaf index \
  --output mouse-2024-A_splici \
  --fasta /path/to/refdata-gex-GRCm39-2024-A/fasta/genome.fa \
  --gtf /path/to/refdata-gex-GRCm39-2024-A/genes/genes.gtf \
  --ref-type spliced+intronic --rlen 91 --threads 16 \
  --work-dir /tmp/${USER}_simpleaf_index
# gex.index = mouse-2024-A_splici/index
```

ADT from a 10x feature CSV (`id`, `name`, `sequence`). Repo copy: `refs/new_feature_ref_quant.csv` (15 bp spatial + antibody). k=7, minimizer 5.

```bash
simpleaf index \
  --output adt_feature_index \
  --feature-csv refs/new_feature_ref_quant.csv \
  --kmer-length 7 --minimizer-length 5 --overwrite --threads 16
# adt.index = adt_feature_index/index
```

## VDJ (TCR / BCR)

10x Immune Profiling (typically **5′**) adds VDJ-T and VDJ-B FASTQs that use the same barcode+UMI chemistry as GEX. Cell Ranger `vdj` / `multi` maps those reads to a **V(D)J reference** (`refdata-cellranger-vdj-…/fasta/regions.fa`) and then **assembles** contigs, CDR3s, and clonotypes.

Sailfish uses the same 10x `regions.fa` but maps with **simpleaf / piscem / alevin-fry** (k=21, minimizer 11, `--keep-duplicates`). It does **not** assemble CDR3s or clonotypes. Per barcode it reports the dominant V / J / C gene (UMI) at each locus and a coarse paired flag (TRA V+J and TRB V+J for TCR; IGH V+J plus IGK or IGL V+J for BCR).

Index is built under `{output}/vdj/vdj_index/` from `vdj.reference` (the 10x ref directory or a `regions.fa`). D-segments and sequences shorter than *k* are dropped. `--t2g-map` collapses alleles to 10x gene names.

If TCR/BCR FASTQs are missing, GEX (±ADT, ±OCM) still runs and VDJ is skipped with a warning.

```json
"vdj": {
  "reference": "/path/to/refdata-cellranger-vdj-GRCm38-alts-ensembl-7.0.0",
  "chemistry": "10xv3-5p",
  "tcr": { "fastqs": "/path/to/E23TB_TCR" },
  "bcr": { "fastqs": "/path/to/E23TB_BCR" }
}
```

`tcr` / `bcr` may instead set `reads1` / `reads2` (Illumina `*_R1_001.fastq.gz` / `*_R2_001.fastq.gz`). Optional `vdj.index` points at a prebuilt `…/index` folder. Example: `configs/e23tb_ocm_vdj.json`.

**Outputs:** `{output}/vdj/{tcr,bcr}_quant/` (simpleaf gene-level counts), `{tcr,bcr}_annotations.csv`, `vdj_run.json`. Columns `tcr_*` / `bcr_*` are attached to the GEX/OCM h5ad `obs`.

## Config

Omit `adt` for GEX-only; omit `ocm` (or `"mode": "quant"`) to skip demux. Default OCM map is OB1–OB4 if `ocm.samples` is missing.

```json
{
  "sample": "E28S",
  "mode": "ocm",
  "output": "runs/E28S",
  "threads": 16,
  "min_gex_umi": 500,
  "gex": {
    "reads1": ["fastqs/E28S_GEX_S3_R1_001.fastq.gz"],
    "reads2": ["fastqs/E28S_GEX_S3_R2_001.fastq.gz"],
    "index": "mouse-2024-A_splici/index",
    "chemistry": "10xv4-3p"
  },
  "adt": {
    "reads1": ["fastqs/E28S_ADT_S4_R1_001.fastq.gz"],
    "reads2": ["fastqs/E28S_ADT_S4_R2_001.fastq.gz"],
    "index": "adt_feature_index/index",
    "chemistry": "e14s-adt-10xv4",
    "geometry": "1{b[16]u[12]x:}2{r[15]x:}",
    "feature_ref": "refs/new_feature_ref_quant.csv"
  },
  "ocm": {
    "samples": {
      "GT": "edge_gfp_plus",
      "CA": "edge_gfp_minus",
      "TC": "core_gfp_plus",
      "AG": "core_gfp_minus"
    }
  }
}
```

| Field | Default | Meaning |
|-------|---------|---------|
| `sample` | `output` basename | Name in output files |
| `mode` | `quant` | `quant` pooled; `ocm` demux by overhang |
| `output` | required | Run dir |
| `threads` | `16` | simpleaf threads (GEX; ADT ~¼ unless `adt.threads`) |
| `skip_quant` | `false` | Skip mapping; reuse `gex.h5ad` or `{output}/gex_quant` |
| `min_gex_umi` | `500` | Cell filter (`filter.min_gex_umi` still works) |
| `alevin_fry_home` | `{output}/af_home` | Chemistry registry + permit lists |
| `tools.conda` | `$SAILFISH_CONDA` or sibling `conda_env` | Env with simpleaf/piscem/alevin-fry |
| `tools.cargo` | `$SAILFISH_CARGO` or sibling `cargo_tools` | Extra `bin/` (e.g. `roers`) |
| `gex.reads1` / `reads2` | | R1/R2 FASTQs (string or list, paired order) |
| `gex.index` | | `…/index` from `simpleaf index` |
| `gex.chemistry` | `10xv3` | e.g. `10xv3`, `10xv4-3p` |
| `gex.min_reads` | `10` | `--min-reads` with unfiltered permit list |
| `gex.resolution` | `cr-like` | UMI resolution (`cr-like`, `cr-like-em`, …) |
| `gex.h5ad` | | Existing GEX `quants.h5ad` (skip GEX map) |
| `gex.geometry` | | Register `gex.chemistry` with this geometry |
| `vdj.reference` | | 10x `refdata-cellranger-vdj-*` dir or `fasta/regions.fa` |
| `vdj.chemistry` | GEX chemistry | Usually `10xv3-5p` (same barcode/UMI as 5′ GEX) |
| `vdj.tcr` / `vdj.bcr` | omitted | `fastqs` dir or `reads1`/`reads2`; skipped if files are missing |
| `vdj.index` | `{output}/vdj/vdj_index/index` | Prebuilt simpleaf VDJ index |
| `vdj.kmer_length` | `21` | piscem k for V/J/C segments |
| `vdj.minimizer_length` | `11` | piscem minimizer |
| `adt` | omitted | Same keys as `gex` |
| `adt.feature_ref` | | Feature CSV to rename ADT (`id`/`name`) |
| `adt.geometry` | | Custom chemistry (needed for `e14s-adt-10xv4`) |
| `ocm.samples` | GEM-X OB1–OB4 | Overhang → id `{"GT": "edge_gfp_plus"}`, or `{"OB1": "…"}`, or `[{sample_id, overhang}]` |
| `ocm.overhang_start` | `7` | 0-based overhang start in the 16 bp barcode |
| `ocm.overhang_len` | `2` | Overhang length |
| `ocm.include_unassigned` | `false` | Keep other overhangs as `unassigned` |
| `cellranger.per_sample_outs` | | `cellranger multi` `per_sample_outs` (OCM compare) |
| `cellranger.filtered_mtx` | | `cellranger count` `filtered_feature_bc_matrix` |
| `existing_quants.gex_h5ad` / `adt_h5ad` | | Same as `gex.h5ad` / `adt.h5ad` |

**Outputs:** quant → `{output}/quant/{sample}_gex[_adt].h5ad` + MTX. OCM → `{output}/ocm/per_sample_outs/<id>/` and `{sample}_gex_adt_ocm.h5ad`. VDJ → `{output}/vdj/` and `tcr_*`/`bcr_*` columns on those h5ads. `cellranger.*` adds `compare_to_cellranger.json`.

## Microwell placement

`placement/cell_placement.py` is the entropy localizer used by
https://lucid-crystal-kmqy.here.now/ (MAP well + discrete/spatial entropy).
Layout: `placement/layouts/chip_layout.xlsx` (source oligos) and
`placement/layouts/data.js` (published 48×48 barcode map). See `placement/README.md`.

