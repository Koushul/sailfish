#!/usr/bin/env bash
set -euo pipefail
PIPE=/ix1/ylee/kor11/tools/af_tutorial/pipeline
export PATH="/ix1/ylee/kor11/tools/af_tutorial/conda_env/bin:/ix1/ylee/kor11/tools/af_tutorial/cargo_tools/bin:$PATH"
export LD_LIBRARY_PATH="/ix1/ylee/kor11/tools/af_tutorial/conda_env/lib:${LD_LIBRARY_PATH:-}"
export PYTHONNOUSERSITE=1

python "$PIPE/run.py" run --config "$PIPE/configs/e28s_quant.json"
python "$PIPE/run.py" run --config "$PIPE/configs/e28s_ocm.json"
echo "E28S quant + OCM done"
