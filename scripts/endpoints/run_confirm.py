#!/usr/bin/env python3
"""Confirmation run (ANALYSIS_PLAN §3 step 5, §4, §5). Run ONCE. Any rerun is a DEVIATION.

  CONFIRM_RUN=1 python scripts/endpoints/run_confirm.py [--root ROOT]
  # testing only (ROOT must differ from the project root; guard failures are reported, not fatal):
  CONFIRM_RUN=1 python scripts/endpoints/run_confirm.py --root TEST_ROOT --test-root [--provisional] [--quick]

Guards, all evaluated BEFORE any data file is read (failures are listed together):
  G1 tuning_no_nulls        config/tuning.yaml has no null value (recursive)
  G2 git_clean              `git status --porcelain --untracked-files=no` is empty
  G3 tuning_tracked         config/tuning.yaml is tracked by git
  G4 plan_md5_ledger        md5(ANALYSIS_PLAN.md) equals the value of the LAST ledger row whose quantity is
                            'freeze: ANALYSIS_PLAN.md md5'
  G5 tuning_md5_ledger      if a 'freeze: config/tuning.yaml md5' row exists: it equals md5(config/tuning.yaml)
  G6 split_md5_ledger       if a 'freeze: data/interim/split.tsv md5' row exists: it equals md5(split.tsv)
  G7 commit_ledger          if a 'freeze: commit' row exists: it equals `git rev-parse HEAD`
  G8 tuning_split_md5       tuning.yaml split_md5 equals md5(data/interim/split.tsv)
  G9 confirm_run_env        CONFIRM_RUN=1
  G10 not_a_rerun           results/endpoints/ has no CONFIRM_RUN.json / E*.json, unless
                            --rerun-logged-deviation <ID> and <ID> appears in DEVIATIONS.md

Ledger row format (RESULTS_LEDGER.md, markdown table | date | quantity | value | script | commit | input md5 | reviewed by |):
  | 2026-10-31 | freeze: ANALYSIS_PLAN.md md5 | <32 hex> | scripts/endpoints/run_confirm.py --print-freeze-rows | <commit> | - | PI |
`--print-freeze-rows` prints these rows for the current state (read-only; paste them into the ledger at the freeze).

Writes ROOT/results/endpoints/{E1..E6}.json (§4 key set), holm.json, interpretation.json, CONFIRM_RUN.json.
"""
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "2")
import argparse, json, re, subprocess, sys, time
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT / "src"))
from strucevo.stats.common import md5_file, write_json  # noqa: E402  (no data access)
from strucevo.stats.driver import tuning_nulls  # noqa: E402

FREEZE = {"plan": "freeze: ANALYSIS_PLAN.md md5", "tuning": "freeze: config/tuning.yaml md5",
          "split": "freeze: data/interim/split.tsv md5", "commit": "freeze: commit"}


def _git(root, *args):
    try:
        r = subprocess.run(["git", "-C", str(root), *args], capture_output=True, text=True, timeout=30)
        return r.returncode, r.stdout.strip(), r.stderr.strip()
    except Exception as e:
        return 1, "", str(e)


def ledger_values(root):
    """{quantity: value} from the LAST ledger row per quantity (markdown table rows)."""
    p = Path(root) / "RESULTS_LEDGER.md"
    out = {}
    if not p.exists():
        return out
    for line in p.read_text().splitlines():
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) >= 3 and cells[1].lower().startswith("freeze:"):
            out[cells[1]] = cells[2].strip("` ")
    return out


def evaluate_guards(root, rerun_id=None):
    root = Path(root)
    g = {}
    nulls = tuning_nulls(root)
    g["G1_tuning_no_nulls"] = (not nulls, f"{len(nulls)} null value(s): {nulls[:8]}{' ...' if len(nulls) > 8 else ''}" if nulls else "ok")
    rc, out, err = _git(root, "status", "--porcelain", "--untracked-files=no")
    g["G2_git_clean"] = (rc == 0 and out == "", "ok" if (rc == 0 and out == "") else (err or f"dirty: {out.splitlines()[:5]}"))
    rc, out, err = _git(root, "ls-files", "--error-unmatch", "config/tuning.yaml")
    g["G3_tuning_tracked"] = (rc == 0, "ok" if rc == 0 else (err or "config/tuning.yaml not tracked"))
    led = ledger_values(root)
    pm = md5_file(root / "ANALYSIS_PLAN.md")
    lv = led.get(FREEZE["plan"])
    g["G4_plan_md5_ledger"] = (pm is not None and lv == pm,
                               "ok" if (pm is not None and lv == pm) else f"plan md5 {pm} vs ledger {lv!r} (row '{FREEZE['plan']}')")
    for key, gid, path in (("tuning", "G5_tuning_md5_ledger", "config/tuning.yaml"),
                           ("split", "G6_split_md5_ledger", "data/interim/split.tsv")):
        lv = led.get(FREEZE[key]); m = md5_file(root / path)
        g[gid] = (lv is None or lv == m, "no ledger row (not checked)" if lv is None else ("ok" if lv == m else f"{path} md5 {m} vs ledger {lv}"))
    lv = led.get(FREEZE["commit"])
    rc, head, _ = _git(root, "rev-parse", "HEAD")
    g["G7_commit_ledger"] = (lv is None or (rc == 0 and head.startswith(lv)),
                             "no ledger row (not checked)" if lv is None else ("ok" if head.startswith(lv) else f"HEAD {head} vs ledger {lv}"))
    try:
        import yaml
        tsm = (yaml.safe_load((root / "config" / "tuning.yaml").read_text()) or {}).get("split_md5")
    except Exception:
        tsm = None
    sm = md5_file(root / "data" / "interim" / "split.tsv")
    g["G8_tuning_split_md5"] = (tsm is not None and tsm == sm, "ok" if (tsm is not None and tsm == sm) else f"tuning.yaml split_md5 {tsm} vs {sm}")
    g["G9_confirm_run_env"] = (os.environ.get("CONFIRM_RUN") == "1", "ok" if os.environ.get("CONFIRM_RUN") == "1" else "CONFIRM_RUN != 1")
    ep = root / "results" / "endpoints"
    prev = sorted(p.name for p in ep.glob("E*.json")) + (["CONFIRM_RUN.json"] if (ep / "CONFIRM_RUN.json").exists() else [])
    if not prev:
        g["G10_not_a_rerun"] = (True, "ok (first run)")
    else:
        dev = (root / "DEVIATIONS.md").read_text() if (root / "DEVIATIONS.md").exists() else ""
        ok = bool(rerun_id) and re.search(r"(^|\W)" + re.escape(rerun_id) + r"(\W|$)", dev) is not None
        g["G10_not_a_rerun"] = (ok, f"rerun under logged deviation {rerun_id}" if ok else
                                f"previous confirmation outputs exist ({prev[:4]}); a rerun requires --rerun-logged-deviation <ID> with <ID> in DEVIATIONS.md")
    return {k: {"pass": bool(v[0]), "detail": v[1]} for k, v in g.items()}


def freeze_rows(root):
    root = Path(root)
    rc, head, _ = _git(root, "rev-parse", "HEAD")
    d = time.strftime("%Y-%m-%d")
    rows = [(FREEZE["plan"], md5_file(root / "ANALYSIS_PLAN.md")), (FREEZE["tuning"], md5_file(root / "config/tuning.yaml")),
            (FREEZE["split"], md5_file(root / "data/interim/split.tsv")), (FREEZE["commit"], head if rc == 0 else None)]
    return [f"| {d} | {q} | {v} | scripts/endpoints/run_confirm.py --print-freeze-rows | {head if rc == 0 else '-'} | - | |" for q, v in rows]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=str(PROJECT))
    ap.add_argument("--test-root", action="store_true", help="TESTING ONLY (ROOT != project root): guards reported, not fatal")
    ap.add_argument("--provisional", action="store_true", help="TESTING ONLY: provisional tuning values for nulls")
    ap.add_argument("--quick", action="store_true", help="TESTING ONLY: small bootstrap sizes")
    ap.add_argument("--rerun-logged-deviation", default=None, metavar="ID")
    ap.add_argument("--print-freeze-rows", action="store_true")
    ap.add_argument("--no-alt", action="store_true", help="TESTING ONLY: skip alternative alignments")
    ap.add_argument("--device", default=None)
    ap.add_argument("--size", action="append", default=[], metavar="KEY=INT",
                    help="override one bootstrap size (keys: see driver.PLAN_SIZES) -- TESTING ONLY (--test-root)")
    a = ap.parse_args()
    root = Path(a.root).resolve()
    if a.print_freeze_rows:
        print("\n".join(freeze_rows(root)))
        return 0
    is_project = root == PROJECT.resolve()
    if a.test_root and is_project:
        sys.exit("run_confirm.py: --test-root is not permitted on the project root")
    for opt in ("provisional", "quick", "no_alt", "size"):
        if getattr(a, opt) and not a.test_root:
            sys.exit(f"run_confirm.py: --{opt.replace('_', '-')} is permitted only with --test-root")
    guards = evaluate_guards(root, a.rerun_logged_deviation)
    failed = [k for k, v in guards.items() if not v["pass"]]
    print("Guards:")
    for k, v in guards.items():
        print(f"  [{'PASS' if v['pass'] else 'FAIL'}] {k}: {v['detail']}")
    if failed and not a.test_root:
        print(f"\nREFUSED: {len(failed)} guard(s) failed: {', '.join(failed)}. No data was read.", file=sys.stderr)
        return 2
    if failed:
        print(f"\nTEST ROOT: {len(failed)} guard(s) failed ({', '.join(failed)}) -- continuing because --test-root.")
    if os.environ.get("CONFIRM_RUN") != "1":       # data loaders raise anyway; fail early and clearly
        print("REFUSED: CONFIRM_RUN=1 is required to read split == 'confirm'", file=sys.stderr)
        return 2
    # ------------------------------------------------------------------ run
    import torch
    torch.set_num_threads(2)
    from strucevo.stats.driver import EndpointDriver, QUICK_SIZES, finalize_E5
    from strucevo.stats.endpoints_misc import apply_holm
    sys.path.insert(0, str(PROJECT / "scripts" / "endpoints"))
    from interpret import interpret_dir
    out = root / "results" / "endpoints"
    out.mkdir(parents=True, exist_ok=True)
    rc, head, _ = _git(root, "rev-parse", "HEAD")
    log = {"started": time.strftime("%Y-%m-%dT%H:%M:%S"), "root": str(root), "test_root": a.test_root,
           "guards": guards, "guards_failed": failed, "rerun_logged_deviation": a.rerun_logged_deviation,
           "commit": head if rc == 0 else None, "provisional": a.provisional, "quick": a.quick,
           "plan_md5": md5_file(root / "ANALYSIS_PLAN.md"), "tuning_md5": md5_file(root / "config/tuning.yaml"),
           "split_md5": md5_file(root / "data/interim/split.tsv"), "status": "running"}
    write_json(log, out / "CONFIRM_RUN.json")
    t0 = time.time()
    label = "TEST_ROOT_not_a_confirmation_result" if a.test_root else None
    sizes = dict(QUICK_SIZES) if a.quick else {}
    sizes.update({k: int(v) for k, v in (kv.split("=", 1) for kv in a.size)})
    drv = EndpointDriver(root, split="confirm", provisional=a.provisional, sizes=sizes or None,
                         outdir=out, label=label, device=a.device, run_alt_alignments=not a.no_alt,
                         log=lambda s: print(s, flush=True))
    recs = drv.run_all(compute_power=False)
    holm = apply_holm(out)
    finalize_E5(out)
    interp = interpret_dir(out)
    for e in ("E1", "E2", "E3", "E4", "E5", "E6"):      # stamp guards + commit into every record
        p = out / f"{e}.json"
        if p.exists():
            r = json.loads(p.read_text())
            r["commit"] = log["commit"]; r["confirm_guards_failed"] = failed
            if a.rerun_logged_deviation:
                r["flags"] = list(r.get("flags") or []) + [f"rerun_deviation:{a.rerun_logged_deviation}"]
            write_json(r, p)
    log.update({"status": "done", "finished": time.strftime("%Y-%m-%dT%H:%M:%S"), "runtime_s": round(time.time() - t0, 1),
                "timing_s": drv.timing,
                "decisions": {e: json.loads((out / f"{e}.json").read_text()).get("decision") for e in recs if (out / f"{e}.json").exists()},
                "holm": holm["p_holm"], "interpretation": interp["row_ids"]})
    write_json(log, out / "CONFIRM_RUN.json")
    print(json.dumps({k: log[k] for k in ("decisions", "holm", "interpretation", "runtime_s", "timing_s")}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
