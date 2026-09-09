#!/bin/bash
#SBATCH -J sailfish
#SBATCH -t 4:00:00
#SBATCH -c 64
#SBATCH --mem=200G
#SBATCH -o sailfish_%j.out
#SBATCH -e /dev/null

set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
CONFIG="${1:?usage: sbatch run_sbatch.sh config.json}"
shift || true
python "$ROOT/run.py" "$CONFIG" "$@"
