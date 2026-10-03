#!/usr/bin/env bash
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/lib.sh"
"$PYTHON" scripts/inspect_upstream.py "$@"
