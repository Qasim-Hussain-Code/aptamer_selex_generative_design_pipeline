#!/usr/bin/env bash
set -euo pipefail
# shellcheck source=scripts/lib.sh
source "$(dirname -- "${BASH_SOURCE[0]}")/lib.sh"
if [[ ${1:-} == --help ]]; then
    echo 'Usage: bash scripts/18_verify.sh (runs clean-clone smoke, tests, shellcheck, core restart and full refusal audit)'
    exit 0
fi
stage tests -m pytest -q --junitxml=logs/pytest.xml
mapfile -t SHELL_SCRIPTS < <(find scripts remote -type f -name '*.sh' -print)
shellcheck -x run_all.sh "${SHELL_SCRIPTS[@]}" > logs/shellcheck.txt 2>&1
# These commands prove restart and explicit refusal; stdout is retained for the audit.
bash run_all.sh --mode core --yes > logs/core_verification.txt 2>&1
if bash run_all.sh --mode full --gpu-mode hosted --from 11 --yes > logs/full_verification.txt 2>&1; then
    echo 'Unexpected full-mode success without all benchmark arms' >&2
    exit 1
fi
stage clean_clone scripts/clean_clone.py
stage verification scripts/verify_repository.py --clean-clone-status passed
