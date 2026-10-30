#!/usr/bin/env bash
# ANALYSIS_PLAN §3 tuning pipeline (split == 'tune' only), in dependency order, with sentinels.
#
#   bash scripts/tuning/run_tuning.sh [--root R] [--only "step1 step2 ..."] [--from STEP] [--force] [--dry-run]
#
# Steps (sentinel ROOT/logs/tuning/<step>.done, log ROOT/logs/tuning/<step>.log; a step with a sentinel is skipped
# unless --force; --only restricts to the named steps; --from starts at STEP):
#   validate_profile_alpha  profile_alpha_validation.py      (validation only; no yaml keys)
#   tune_alpha              tune_alpha.py                    -> results/tuning/alpha.json
#   yaml_alpha              write_tuning_yaml.py             (alpha into tuning.yaml; read by the next steps)
#   floors_tuned            floors_tuned.py                  -> floors_tuned.json (eta, S0 at tuned alpha)
#   select_bootsub          select_bootsub.py                -> data/interim/pairs_bootsub/bootsub_pairs.tsv
#   bootsub_pairs           scripts/pairs/run_pairs.py       (boot_sub pairs + column dumps; $PAIRS_THREADS threads)
#   se_calibration          se_calibration.py                -> se_calibration.json (kappa, rho, lddt SE model)
#   yaml_se                 write_tuning_yaml.py             (eta/S0/kappa/rho/SE model into tuning.yaml)
#   apply_tuning            apply_tuning.py                  (pair tables recomputed IN PLACE at tuned alpha, lddt_D(_se))
#   tau_misfit              tau_misfit.py                    -> tau_misfit.json
#   b_e4_timing             b_e4_timing.py                   -> b_e4_timing.json   (uses tau_misfit)
#   e5_power                e5_power.py                      -> e5_power.json      (uses tau_misfit)
#   v3_tuned                v3_tuned.py                      -> v3_tuned.json
#   sample_throughput       sample_throughput.py             -> sample_throughput.json
#   write_tuning_yaml       write_tuning_yaml.py             (final merge; prints the nulls)
# Filling the plan table is a separate, explicit step (it edits ANALYSIS_PLAN.md):
#   python scripts/tuning/fill_tuning_table.py --plan ROOT/ANALYSIS_PLAN.md --tuning ROOT/config/tuning.yaml
#
# Env: PAIRS_THREADS (default 2 on the shared node; 48 on a free node), STORE / FOLDCOMP (structure store for the
# boot_sub pairs; default ROOT/data/interim/structures/{store,fc}), V3_ARGS / E5_ARGS / E4_ARGS / TAU_ARGS (extra
# flags), PAIRS_RUN_THREADS (threads of the run that wrote data/interim/pairs/*.done; default 48).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT="$(cd "$HERE/../.." && pwd)"
ROOT="$PROJECT"; ONLY=""; FROM=""; FORCE=0; DRY=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --root) ROOT="$(cd "$2" && pwd)"; shift 2;;
    --only) ONLY="$2"; shift 2;;
    --from) FROM="$2"; shift 2;;
    --force) FORCE=1; shift;;
    --dry-run) DRY=1; shift;;
    *) echo "unknown arg $1" >&2; exit 2;;
  esac
done
if [[ -z "${CONDA_PREFIX:-}" || "$CONDA_PREFIX" != "$PROJECT/conda_envs/main" ]]; then
  source /opt/conda/etc/profile.d/conda.sh && conda activate "$PROJECT/conda_envs/main"
fi
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-2} OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-2} \
       MKL_NUM_THREADS=${MKL_NUM_THREADS:-2} POLARS_MAX_THREADS=${POLARS_MAX_THREADS:-2}
PAIRS_THREADS=${PAIRS_THREADS:-2}
STORE=${STORE:-$ROOT/data/interim/structures/store}
FOLDCOMP=${FOLDCOMP:-$ROOT/data/interim/structures/fc}
PY="python"
S="$PROJECT/scripts/tuning"
LOG="$ROOT/logs/tuning"; mkdir -p "$LOG"
started=0; [[ -z "$FROM" ]] && started=1

step() {   # step NAME CMD...
  local name="$1"; shift
  [[ "$name" == "$FROM" ]] && started=1
  if [[ $started -eq 0 ]]; then return 0; fi
  if [[ -n "$ONLY" && " $ONLY " != *" $name "* ]]; then return 0; fi
  if [[ -f "$LOG/$name.done" && $FORCE -eq 0 ]]; then echo "[skip] $name (sentinel $LOG/$name.done)"; return 0; fi
  echo "[run ] $name: $*"
  if [[ $DRY -eq 1 ]]; then return 0; fi
  local t0=$(date +%s)
  if ! "$@" > "$LOG/$name.log" 2>&1; then
    echo "[FAIL] $name -- see $LOG/$name.log"; tail -20 "$LOG/$name.log"; exit 1
  fi
  local t1=$(date +%s)
  printf '{"step": "%s", "finished": "%s", "seconds": %d, "root": "%s"}\n' "$name" "$(date -Iseconds)" $((t1 - t0)) "$ROOT" > "$LOG/$name.done"
  echo "[done] $name ($((t1 - t0)) s)"
}

cd "$PROJECT"
step validate_profile_alpha $PY "$S/profile_alpha_validation.py" --root "$ROOT"
step tune_alpha            $PY "$S/tune_alpha.py" --root "$ROOT"
step yaml_alpha            $PY "$S/write_tuning_yaml.py" --root "$ROOT"
step floors_tuned          $PY "$S/floors_tuned.py" --root "$ROOT"
step select_bootsub        $PY "$S/select_bootsub.py" --root "$ROOT"
step bootsub_pairs         $PY "$PROJECT/scripts/pairs/run_pairs.py" --pairs "$ROOT/data/interim/pairs_bootsub/bootsub_pairs.tsv" \
                               --foldcomp "$FOLDCOMP" --store "$STORE" --out "$ROOT/data/interim/pairs_bootsub" \
                               --dump-columns "$ROOT/data/interim/pairs_bootsub/columns" --n-shards 8 --threads "$PAIRS_THREADS"
step se_calibration        $PY "$S/se_calibration.py" --root "$ROOT"
step yaml_se               $PY "$S/write_tuning_yaml.py" --root "$ROOT"
step apply_tuning          $PY "$S/apply_tuning.py" --root "$ROOT"
step tau_misfit            $PY "$S/tau_misfit.py" --root "$ROOT" ${TAU_ARGS:-}
step b_e4_timing           $PY "$S/b_e4_timing.py" --root "$ROOT" ${E4_ARGS:-}
step e5_power              $PY "$S/e5_power.py" --root "$ROOT" ${E5_ARGS:-}
step v3_tuned              $PY "$S/v3_tuned.py" --root "$ROOT" ${V3_ARGS:-}
step sample_throughput     $PY "$S/sample_throughput.py" --root "$ROOT" --threads-used "${PAIRS_RUN_THREADS:-48}"
step write_tuning_yaml     $PY "$S/write_tuning_yaml.py" --root "$ROOT"
if [[ $DRY -eq 0 && ( -z "$ONLY" || " $ONLY " == *" write_tuning_yaml "* ) ]]; then
  echo "--- nulls in $ROOT/config/tuning.yaml (from the last write):"; sed -n '/^NULLS/,$p' "$LOG/write_tuning_yaml.log" 2>/dev/null || true
fi
echo "next (explicit): python scripts/tuning/fill_tuning_table.py --plan $ROOT/ANALYSIS_PLAN.md --tuning $ROOT/config/tuning.yaml"
