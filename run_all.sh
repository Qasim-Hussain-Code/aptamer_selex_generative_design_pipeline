#!/usr/bin/env bash
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/scripts/lib.sh"
MODE=core FROM=0 TARGET=all SEED=1729
CONFIG_ARGS=()
while (($#)); do
    case "$1" in
        --mode) MODE=$2; shift 2 ;;
        --from) FROM=$2; shift 2 ;;
        --target) TARGET=$2; shift 2 ;;
        --seed) SEED=$2; shift 2 ;;
        --threads|--ram|--disk|--gpu-mode) CONFIG_ARGS+=("$1" "$2"); shift 2 ;;
        --yes) CONFIG_ARGS+=(--yes); shift ;;
        --help) echo 'Usage: bash run_all.sh --mode core|full|smoke --from 00..18 --target all|tg2|integrin --threads N --ram GB --disk GB --gpu-mode none|hosted|local --seed INTEGER --yes'; exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done
[[ $MODE == core || $MODE == full || $MODE == smoke ]] || { echo 'Invalid mode' >&2; exit 2; }
[[ $TARGET == all || $TARGET == tg2 || $TARGET == integrin ]] || { echo 'Invalid target' >&2; exit 2; }
[[ $FROM =~ ^[0-9]+$ && $SEED =~ ^[0-9]+$ ]] || { echo 'Stage and seed must be integers' >&2; exit 2; }
FROM=$((10#$FROM))
((FROM <= 18)) || { echo 'Stage must be between 00 and 18' >&2; exit 2; }
if ((FROM == 18)); then
    bash scripts/18_verify.sh
    exit 0
fi
# Serial stages give stable measurements and avoid simultaneous large temporary files.
if [[ $MODE == smoke ]]; then
    stage smoke scripts/smoke.py
    exit 0
fi
bash scripts/00_configure.sh "${CONFIG_ARGS[@]}"
source scripts/lib.sh
export MASTER_SEED="$SEED"
stage resolve scripts/resolve_config.py
if ((FROM <= 1)); then
    stage bootstrap_sources scripts/bootstrap_sources.py
    for part in packages metadata upstream supplement; do
        bash scripts/acquire_reference_inputs.sh --part "$part"
    done
    # Initial repository records and paper XML must be present; this bootstrap is restartable.
    stage sources scripts/01_verify_sources.py --offline
fi
if ((FROM <= 3)); then stage metadata scripts/03_fetch_metadata.py; fi
if ((FROM <= 4)); then bash scripts/04_fetch_selex.sh --target "$TARGET"; fi
if ((FROM <= 6)); then stage ground_truth scripts/06_build_ground_truth.py; fi
# Target selection limits retrieval. The publication comparison always requires both
# targets so a partial fetch cannot silently replace the other target's analysis.
if ((FROM <= 7)); then stage splits scripts/07_build_splits.py; fi
if ((FROM <= 8)); then stage baselines scripts/08_baselines.py; fi
if ((FROM <= 9)); then stage secondary scripts/09_secondary_structure.py; fi
if ((FROM <= 10)); then stage gpu_jobs scripts/10_prepare_gpu_jobs.py; fi
if [[ $MODE == full ]]; then stage full_preflight scripts/full_preflight.py; fi
if ((FROM <= 13)); then stage ablation scripts/13_structure_ablation.py; fi
if ((FROM <= 12 || FROM == 13)); then stage candidates scripts/12_score_candidates.py; fi
if ((FROM <= 14)); then stage tertiary scripts/14_optional_tertiary.py; fi
if ((FROM <= 15)); then stage evaluation scripts/15_evaluate.py; fi
if ((FROM <= 15)); then stage tests -m pytest -q; fi
if ((FROM <= 16)); then stage figures scripts/16_figures.py; fi
if ((FROM <= 17)); then stage report scripts/write_report.py; fi
echo "Completed $MODE workflow. Run bash scripts/18_verify.sh for the independent audit."
