#!/usr/bin/env bash
set -euo pipefail
if [[ ${1:-} == --help || $# -ne 2 ]]; then
    echo 'Usage: bash remote/raptgen/run_remote.sh BUNDLE OUTPUT_DIRECTORY'
    exit 0
fi
BUNDLE=$(realpath "$1")
OUTPUT=$(mkdir -p "$2" && realpath "$2")
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
UPSTREAM="$OUTPUT/upstream"
if [[ ! -d "$UPSTREAM/.git" ]]; then
    git clone --filter=blob:none https://github.com/hmdlab/raptgen.git "$UPSTREAM"
fi
git -C "$UPSTREAM" checkout --detach c4986ca9fa439b9389916c05829da4ff9c30d6f3
export PYTHONPATH="$UPSTREAM${PYTHONPATH:+:$PYTHONPATH}"
export OMP_NUM_THREADS=1 MPLBACKEND=Agg
python "$SCRIPT_DIR/adapter.py" --bundle "$BUNDLE" --output "$OUTPUT"
