#!/usr/bin/env bash
set -euo pipefail
PIPE="$(cd "$(dirname "$0")" && pwd)"
python "$PIPE/run.py" "$PIPE/configs/e28s_quant.json"
python "$PIPE/run.py" "$PIPE/configs/e28s_ocm.json"
echo "E28S quant + OCM done"
