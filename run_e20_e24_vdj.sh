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
LOGS=/ix1/ylee/kor11/TCR_BCR/logs
mkdir -p "$LOGS" /ix1/ylee/kor11/TCR_BCR/runs
cd "$ROOT"
EXTRA=("$@")

run_one() {
  local cfg="$1" receptor="$2"
  local name
  name="$(basename "$cfg" .json)_${receptor}"
  echo "START $(date -Is) $name"
  "$PY" "$ROOT/run.py" "$cfg" --vdj-receptor "$receptor" "${EXTRA[@]}" > "$LOGS/${name}.log" 2>&1
  echo "DONE $(date -Is) $name"
}

C="$ROOT/configs/e20_e24"
run_one "$C/e20t.json" tcr &
run_one "$C/e21t.json" tcr &
run_one "$C/e22b.json" bcr &
run_one "$C/e23tb.json" tcr &
run_one "$C/e23tb.json" bcr &
run_one "$C/e24t.json" tcr &
wait
echo "VDJ DONE $(date -Is)"
"$PY" "$ROOT/run.py" "$C/ocm_split.json" --split-vdj-ocm > "$LOGS/ocm_split.log" 2>&1
echo "ALL DONE $(date -Is)"
