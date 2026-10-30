#!/usr/bin/env bash
# Post-Ancestra pipeline (ANALYSIS_PLAN §1.6 steps 2-7, §1.2; X4/H1e). Resumable: every step writes a sentinel
# $ROOT/logs/pipeline/<step>.done and is skipped on rerun; pair shards additionally have their own .done sentinels
# (a step interrupted mid-way resumes at the first unfinished shard).
#
#   scripts/run_after_ancestra.sh                 # production (euk180), 48 threads + GPU
#   ROOT=data/interim/itest_euk20 ORTHOXML=results/ancestra/runs/euk20_accurate/ancestra/hogs.orthoxml \
#     FASTA_DIR=data/interim/ancestra/euk20 THREADS=8 scripts/run_after_ancestra.sh       # euk20 test root
#   STEPS="b1 b2" scripts/run_after_ancestra.sh   # run a subset (default: all, in order)
#   ALLOW_BUILD=0 scripts/run_after_ancestra.sh   # never run build_hogs here: adopt an existing build (hogs.parquet,
#                                                 # split.tsv, pairlists/, hogs_build_summary.json) or stop (exit 3)
#   scripts/run_after_ancestra.sh status          # progress + ETA of the pair steps
#
# Steps
#   a   labels: NOT the O(Σn²) all-pairs table. build_hogs labels primary/E6 pairs on demand (HogIndex, proven
#       identical to hog_pair_labels.label_pairs); H1e labels only its candidate pairs (step d). Step a only
#       checks the OrthoXML and records its size.
#   b1  build_hogs --projection-only (counts per stratum, primary/E6 pairs; C25 caps)     -> ETA for step c
#   b2  build_hogs (real; writes split.tsv ONCE, hogs.parquet, pairlists/)                 [C25: E6 caps 5/5]
#   c0  per-split pair lists + AFDB v4 FoldComp fetch (store is built by the first run_pairs call)
#   c1  tune:    seqonly+tdiaa per shard  (seqonly must precede the others in each shard: it supplies aa_ml_seq)
#   c2  tune:    tmalign
#   c3  confirm: seqonly+tdiaa
#   c4  confirm: tmalign
#   d   H1e: SIFTS ortholog pairs with qualifying PDB chains on both sides (Ancestra 1:1 labels) -> AF-AF (h1e_af),
#       PDB-PDB (h1e_pdb) [+ AF-PDB floor] through the pipeline into data/interim/pairs_h1e/
# Alignment order: the PI asked for tdiaa -> seqonly -> tmalign. tdiaa tables join aa_ml_seq from the seqonly shard,
# so each shard runs seqonly immediately before tdiaa (seqonly is the cheapest stage); tmalign (3x slower) runs last.
# Tune shards: ids 0..N_TUNE-1; confirm shards: 1000.. (--shard-offset), all in data/interim/pairs/{aln}/.
set -euo pipefail
PROJ="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJ"
ROOT="${ROOT:-$PROJ}"
case "$ROOT" in /*) ;; *) ROOT="$PROJ/$ROOT" ;; esac
ORTHOXML="${ORTHOXML:-$PROJ/data/interim/ancestra_euk180/ancestra/hogs.orthoxml}"
FASTA_DIR="${FASTA_DIR:-$PROJ/data/interim/ancestra_input}"
THREADS="${THREADS:-48}"
E6CAPS="${E6CAPS:---e6-max-o 5 --e6-max-per-dcg 5}"     # C25
PAIRS_PER_SHARD="${PAIRS_PER_SHARD:-15000}"
STEPS="${STEPS:-a b1 b2 c0 c1 c2 c3 c4 d}"
CONFIRM_OFFSET=1000
H1E_THREADS="${H1E_THREADS:-$THREADS}"

source /opt/conda/etc/profile.d/conda.sh
conda activate "$PROJ/conda_envs/main"
export PYTHONUNBUFFERED=1

DI="$ROOT/data/interim"
S="$ROOT/logs/pipeline"
mkdir -p "$S"
LOG="$S/run_after_ancestra.log"
log() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }

run_step() {   # run_step NAME CMD...
  local name="$1"; shift
  if [[ " $STEPS " != *" $name "* ]]; then return 0; fi
  if [[ -f "$S/$name.done" ]]; then log "skip $name (sentinel $(cat "$S/$name.done"))"; return 0; fi
  log "START $name"
  local t0; t0=$(date +%s)
  set +e
  ( set -e; "$@" ) 2>&1 | tee -a "$S/$name.log"
  local rc=${PIPESTATUS[0]}
  set -e
  if [[ $rc -ne 0 ]]; then log "FAILED $name (rc $rc); see $S/$name.log"; exit $rc; fi
  local dt=$(( $(date +%s) - t0 ))
  echo "{\"step\": \"$name\", \"wall_s\": $dt, \"finished\": \"$(date -Is)\", \"threads\": $THREADS}" > "$S/$name.done"
  log "DONE $name in ${dt}s"
}

npairs() { python -c "import json;print(json.load(open('$DI/pairruns/counts.json'))['$1']['pairs'])"; }
nshards() { python -c "import math;print(max(1, math.ceil($1/$PAIRS_PER_SHARD)))"; }

status() {
  [[ -f "$DI/pairruns/counts.json" ]] || { echo "no pair lists yet"; return 0; }
  local nt nc; nt=$(npairs tune); nc=$(npairs confirm)
  python scripts/pipeline/eta.py progress --out "$DI/pairs" --alns seqonly tdiaa tmalign --offset-range 0 999 --expected "$nt"
  python scripts/pipeline/eta.py progress --out "$DI/pairs" --alns seqonly tdiaa tmalign --offset-range 1000 1999 --expected "$nc"
  ls "$S"/*.done 2>/dev/null | xargs -n1 basename
}
if [[ "${1:-}" == "status" ]]; then status; exit 0; fi

log "ROOT=$ROOT ORTHOXML=$ORTHOXML THREADS=$THREADS STEPS='$STEPS'"
[[ "$ROOT" != "$PROJ" ]] && log "NOTE: test root; production files under $PROJ/data/interim are not touched"

# ---------------------------------------------------------------- a
step_a() {
  [[ -s "$ORTHOXML" ]] || { echo "OrthoXML missing: $ORTHOXML"; return 1; }
  echo "orthoxml $(du -h "$ORTHOXML" | cut -f1) $(date -r "$ORTHOXML" -Is)"
  echo "labels: on demand inside build_hogs (HogIndex) and H1e (candidate pairs only); no all-pairs table"
}
run_step a step_a

# ---------------------------------------------------------------- b
BH="python scripts/hogs/build_hogs.py --orthoxml $ORTHOXML --outdir $DI --fasta-dir $FASTA_DIR --threads $THREADS $E6CAPS"
[[ "$ROOT" != "$PROJ" ]] && BH="$BH --logs $ROOT/logs/filters"
step_b1() {
  $BH --projection-only
  python - "$DI/hogs_projection_summary.json" "$THREADS" <<'EOF'
import json, sys, subprocess
d = json.load(open(sys.argv[1])); f = d.get("final", {})
print("[projection]", json.dumps({k: f.get(k) for k in ("n_candidates", "n_selected", "per_stratum_eligible",
      "per_stratum_sampled", "primary_pairs", "e6_pairs", "total_pairs")}))
n = f.get("total_pairs") or 0
for lab, frac in (("tune(~20%)", 0.2), ("confirm(~80%)", 0.8)):
    subprocess.run(["python", "scripts/pipeline/eta.py", "predict", "--pairs", str(int(n * frac)), "--alns",
                    "seqonly", "tdiaa", "tmalign", "--threads", sys.argv[2], "--label", lab], check=True)
EOF
}
# A build made outside this script (e.g. by hand) is detected and adopted: no rebuild, sentinel written.
# Never start build_hogs while another build_hogs process is running (it would race on split.tsv).
build_running() { pgrep -f "scripts/hogs/build_hogs.py" >/dev/null 2>&1; }
adopt() {   # adopt STEP FILE...  -> write the sentinel if all files exist
  local name="$1"; shift
  [[ -f "$S/$name.done" ]] && return 0
  for f in "$@"; do [[ -s "$f" ]] || return 1; done
  echo "{\"step\": \"$name\", \"adopted_existing\": \"$*\", \"finished\": \"$(date -Is)\"}" > "$S/$name.done"
  log "adopted existing outputs for $name: $*"
}
BUILT=("$DI/hogs.parquet" "$DI/split.tsv" "$DI/pairlists/primary.tsv" "$DI/pairlists/e6.tsv" "$DI/hogs_build_summary.json")
if [[ " $STEPS " == *" b1 "* || " $STEPS " == *" b2 "* ]]; then
  if build_running; then
    if adopt b2 "${BUILT[@]}"; then :; else
      log "build_hogs is running (pid $(pgrep -f scripts/hogs/build_hogs.py | tr '\n' ' ')); waiting for it to finish"
      while build_running; do sleep 30; done
    fi
  fi
  adopt b2 "${BUILT[@]}" && { adopt b1 "${BUILT[@]}" || true; }
  adopt b1 "$DI/hogs_projection_summary.json" || true
  if [[ "${ALLOW_BUILD:-1}" == 0 && ! -f "$S/b2.done" ]]; then
    log "ALLOW_BUILD=0 and no complete build_hogs output in $DI (need: ${BUILT[*]})"; exit 3
  fi
fi
run_step b1 step_b1
run_step b2 $BH

# ---------------------------------------------------------------- c
step_c0() {
  python scripts/pipeline/prep_pair_lists.py --pairlists "$DI/pairlists" --outdir "$DI/pairruns"
  if [[ ! -f "$DI/structures/store/BUILD.json" ]]; then
    mkdir -p "$DI/structures"
    python scripts/pairs/fetch_structures.py "$DI/pairruns/accessions.tsv" "$DI/structures/fc" \
      --procs 8 --threads "$(( THREADS > 8 ? 64 : 16 ))"
  fi
}
run_step c0 step_c0

RP="python scripts/pairs/run_pairs.py --foldcomp $DI/structures/fc --store $DI/structures/store --out $DI/pairs --threads $THREADS"
pairstep() {   # pairstep SPLIT OFFSET ALNS...
  local split="$1" off="$2"; shift 2
  local n ns hi; n=$(npairs "$split")
  ns="${N_SHARDS_OVERRIDE:-}"
  if [[ "$split" == tune ]]; then ns="${N_SHARDS_TUNE:-$(nshards "$n")}"; else ns="${N_SHARDS_CONFIRM:-$(nshards "$n")}"; fi
  hi=$(( off + 999 ))
  # ETA: measured rate from earlier shards of this run (tune measures confirm), else the benchmark
  python scripts/pipeline/eta.py predict --pairs "$n" --alns "$@" --threads "$THREADS" --measured "$DI/pairs" \
    --label "$split [$*] $ns shards"
  $RP --pairs "$DI/pairruns/pairs_${split}.tsv" --n-shards "$ns" --shard-offset "$off" --alns "$@"
  python scripts/pipeline/eta.py rates --out "$DI/pairs" --offset-range "$off" "$hi"
}
run_step c1 pairstep tune 0 seqonly tdiaa
run_step c2 pairstep tune 0 tmalign
run_step c3 pairstep confirm $CONFIRM_OFFSET seqonly tdiaa
run_step c4 pairstep confirm $CONFIRM_OFFSET tmalign

# ---------------------------------------------------------------- d
step_d() {
  [[ -f scripts/h1e/run_h1e.py ]] || { echo "scripts/h1e/run_h1e.py missing"; return 1; }
  python scripts/h1e/run_h1e.py --orthoxml "$ORTHOXML" --root "$ROOT" --threads "$H1E_THREADS" ${H1E_ARGS:-}
}
run_step d step_d

log "pipeline complete; next: scripts/tuning/run_tuning.sh (tune split only), then freeze (ANALYSIS_PLAN §3.4)"
status
