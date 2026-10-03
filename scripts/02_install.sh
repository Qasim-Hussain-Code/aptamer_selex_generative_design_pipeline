#!/usr/bin/env bash
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/lib.sh"
BOOTSTRAP_PYTHON=${BOOTSTRAP_PYTHON:-$PYTHON}
if [[ ${1:-} == --help ]]; then
    echo 'Usage: bash scripts/02_install.sh [--record-only]'
    exit 0
fi
if [[ ${1:-} != --record-only ]]; then
    "$BOOTSTRAP_PYTHON" scripts/install_preflight.py
    if [[ ! -d .venv ]]; then
        "$BOOTSTRAP_PYTHON" -m venv .venv
    fi
    source scripts/lib.sh
    export PIP_CACHE_DIR="$ROOT/.cache/pip" TMPDIR="$ROOT/.cache/tmp" TEMP="$ROOT/.cache/tmp" TMP="$ROOT/.cache/tmp"
    mkdir -p "$TMPDIR"
    "$BOOTSTRAP_PYTHON" scripts/resource_wrapper.py --stage install -- "$PYTHON" -m pip install --no-cache-dir -r config/requirements.txt
fi
stage software scripts/software_versions.py
