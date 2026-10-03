#!/usr/bin/env bash
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/lib.sh"
stage references scripts/acquire_reference_inputs.py "$@"
