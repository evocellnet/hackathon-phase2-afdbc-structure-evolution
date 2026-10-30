#!/usr/bin/env python3
"""Dry run of all endpoints on the TUNING split (ANALYSIS_PLAN §3 step 3). Outputs are EXPLORATORY and are never
pooled with confirmation results. Refuses --split confirm.

  python scripts/endpoints/run_dryrun.py [--root ROOT] [--split tune] [--out ROOT/results/tuning/dryrun]
         [--provisional] [--quick] [--no-alt]

--provisional: fill null tuning.yaml values with kappa=1, rho=0.3, alpha=1, eta/S0 from results/floors/floors.json,
tau_misfit=3, small B_E4; every output records which values were provisional.
"""
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "2")
import argparse, json, sys, time
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT / "src"))
import torch  # noqa: E402
torch.set_num_threads(2)
from strucevo.stats.driver import EndpointDriver, QUICK_SIZES, finalize_E5  # noqa: E402
from strucevo.stats.endpoints_misc import apply_holm  # noqa: E402

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--root", default=str(PROJECT))
ap.add_argument("--split", default="tune")
ap.add_argument("--out", default=None, help="default ROOT/results/tuning/dryrun")
ap.add_argument("--provisional", action="store_true")
ap.add_argument("--quick", action="store_true", help="small bootstrap sizes (testing only)")
ap.add_argument("--no-alt", action="store_true", help="skip seqonly/tmalign sensitivity alignments")
ap.add_argument("--device", default=None)
ap.add_argument("--size", action="append", default=[], metavar="KEY=INT", help="override one bootstrap size (keys: see driver.PLAN_SIZES)")
a = ap.parse_args()
if a.split != "tune":
    sys.exit("run_dryrun.py: refusing split != 'tune' (dry runs are exploratory and tuning-only; use run_confirm.py)")
root = Path(a.root).resolve()
out = Path(a.out) if a.out else root / "results" / "tuning" / "dryrun"
if (root / "results" / "endpoints").resolve() == out.resolve():
    sys.exit("run_dryrun.py: refusing to write dry-run output into results/endpoints/")
t0 = time.time()
sizes = dict(QUICK_SIZES) if a.quick else {}
sizes.update({k: int(v) for k, v in (kv.split("=", 1) for kv in a.size)})
drv = EndpointDriver(root, split="tune", provisional=a.provisional, sizes=sizes or None,
                     outdir=out, label="exploratory_dryrun_tune", device=a.device, run_alt_alignments=not a.no_alt,
                     log=lambda s: print(s, flush=True))
drv.run_all(compute_power=True)
holm = apply_holm(out)
finalize_E5(out)
sys.path.insert(0, str(PROJECT / "scripts" / "endpoints"))
from interpret import interpret_dir  # noqa: E402
interp = interpret_dir(out)
print(json.dumps({"holm": holm, "interpretation": interp, "timing": drv.timing, "runtime_s": time.time() - t0},
                 indent=1, default=str))
