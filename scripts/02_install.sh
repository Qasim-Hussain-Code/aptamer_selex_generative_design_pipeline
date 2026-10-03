#!/usr/bin/env bash
set -euo pipefail
# shellcheck source=scripts/lib.sh
source "$(dirname -- "${BASH_SOURCE[0]}")/lib.sh"
BOOTSTRAP_PYTHON=${BOOTSTRAP_PYTHON:-$PYTHON}
if [[ ${1:-} == --help ]]; then
    echo 'Usage: bash scripts/02_install.sh [--record-only]'
    exit 0
fi
if [[ ${1:-} != --record-only ]]; then
    "$BOOTSTRAP_PYTHON" scripts/resource_wrapper.py --stage package_preflight -- "$BOOTSTRAP_PYTHON" scripts/acquire_reference_inputs.py --part packages
    "$BOOTSTRAP_PYTHON" scripts/install_preflight.py
    if [[ ! -d .venv ]]; then
        "$BOOTSTRAP_PYTHON" -m venv .venv
    fi
    source scripts/lib.sh
    export PIP_CACHE_DIR="$ROOT/.cache/pip" TMPDIR="$ROOT/.cache/tmp" TEMP="$ROOT/.cache/tmp" TMP="$ROOT/.cache/tmp"
    mkdir -p "$TMPDIR"
    requirements=config/requirements.txt
    if [[ -f config/requirements-lock.txt ]]; then requirements=config/requirements-lock.txt; fi
    "$BOOTSTRAP_PYTHON" scripts/resource_wrapper.py --stage install -- "$PYTHON" -m pip install --no-cache-dir -r "$requirements"
fi
stage software scripts/software_versions.py
