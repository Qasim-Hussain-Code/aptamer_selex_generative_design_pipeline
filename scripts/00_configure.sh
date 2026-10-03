#!/usr/bin/env bash
set -euo pipefail
# shellcheck source=scripts/lib.sh
source "$(dirname -- "${BASH_SOURCE[0]}")/lib.sh"
"$PYTHON" scripts/resource_wrapper.py --stage configure -- "$PYTHON" scripts/configure.py "$@"
