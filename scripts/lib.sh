#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
cd "$ROOT"
# An explicit interpreter avoids old molecular-graphics Python installations on PATH.
if [[ -f "$ROOT/.venv/Scripts/python.exe" ]]; then
    PYTHON="$ROOT/.venv/Scripts/python.exe"
elif [[ -f "$ROOT/.venv/bin/python" ]]; then
    PYTHON="$ROOT/.venv/bin/python"
else
    PYTHON=${PYTHON:-python3}
fi
if [[ -f project.conf ]]; then
    # shellcheck disable=SC1091
    source project.conf
fi
export PYTHONPATH="$ROOT/scripts${PYTHONPATH:+:$PYTHONPATH}"
export OMP_NUM_THREADS=${THREADS:-1} OPENBLAS_NUM_THREADS=${THREADS:-1} MKL_NUM_THREADS=${THREADS:-1}
export PIPELINE_BASH_VERSION="$BASH_VERSION"
if command -v cygpath >/dev/null 2>&1; then
    export PIPELINE_BASH_EXECUTABLE
    PIPELINE_BASH_EXECUTABLE=$(cygpath -am "$BASH")
else
    export PIPELINE_BASH_EXECUTABLE="$BASH"
fi
export PATH="$ROOT/.venv/Scripts:$ROOT/.venv/bin:$PATH"
stage() {
    local label=$1
    shift
    "$PYTHON" scripts/resource_wrapper.py --stage "$label" -- "$PYTHON" "$@"
}
