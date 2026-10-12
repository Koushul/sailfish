#!/bin/bash
#SBATCH -J vdj-e20-e24
#SBATCH -p rtx6k
#SBATCH -A ylee
#SBATCH --qos=normal
#SBATCH --gres=gpu:1
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH -t 1-00:00:00
#SBATCH -o /ix1/ylee/kor11/TCR_BCR/logs/vdj-e20-e24-%j.out
#SBATCH -e /ix1/ylee/kor11/TCR_BCR/logs/vdj-e20-e24-%j.err

set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
PY="${PY:-/ix1/ylee/kor11/tools/af_tutorial/conda_env/bin/python}"
mkdir -p /ix1/ylee/kor11/TCR_BCR/logs /ix1/ylee/kor11/TCR_BCR/runs
cd "$ROOT"

run_one() {
  local cfg="$1"
  echo "START $(date -Is) $cfg"
  "$PY" "$ROOT/run.py" "$cfg"
  echo "DONE $(date -Is) $cfg"
}

run_one "$ROOT/configs/e20_e24/e20t.json" &
run_one "$ROOT/configs/e20_e24/e21t.json" &
run_one "$ROOT/configs/e20_e24/e22b.json" &
run_one "$ROOT/configs/e20_e24/e24t.json" &
wait
run_one "$ROOT/configs/e20_e24/e23tb.json"
echo "ALL DONE $(date -Is)"
