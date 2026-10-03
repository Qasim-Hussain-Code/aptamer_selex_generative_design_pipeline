#!/usr/bin/env bash
set -euo pipefail
if [[ ${1:-} == --help ]]; then
    echo 'AptaDiff execution is blocked. Read remote/aptadiff/BLOCKER.md.'
    exit 0
fi
echo 'Refused: inspected AptaDiff commit has no licence file; paper MIT declaration conflicts with repository evidence.' >&2
echo 'No AptaDiff implementation was installed, copied or approximated. See remote/aptadiff/BLOCKER.md.' >&2
exit 78
