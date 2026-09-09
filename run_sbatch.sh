#!/bin/bash
#SBATCH -J sf_pipeline
#SBATCH -M gpu
#SBATCH -p a100_nvlink
#SBATCH -q gpu-a100_nvlink-s
#SBATCH -t 4:00:00
#SBATCH -c 64
#SBATCH --mem=200G
#SBATCH -o /ix1/ylee/kor11/tools/af_tutorial/pipeline/runs/slurm_%j.out
#SBATCH -e /ix1/ylee/kor11/tools/af_tutorial/pipeline/runs/slurm_%j.err

set -euo pipefail
PIPE=/ix1/ylee/kor11/tools/af_tutorial/pipeline
CONFIG="${CONFIG:?set CONFIG=/path/to/config.json}"
export PATH="/ix1/ylee/kor11/tools/af_tutorial/conda_env/bin:/ix1/ylee/kor11/tools/af_tutorial/cargo_tools/bin:$PATH"
export LD_LIBRARY_PATH="/ix1/ylee/kor11/tools/af_tutorial/conda_env/lib:${LD_LIBRARY_PATH:-}"
export PYTHONNOUSERSITE=1
python "$PIPE/run.py" run --config "$CONFIG" "$@"
