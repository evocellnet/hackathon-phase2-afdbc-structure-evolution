#!/usr/bin/env python3
"""Apply the tuned parameters to the pair tables (recompute from counts_*.npz; no realignment).

  python scripts/tuning/apply_tuning.py --root R                                   # in place, alpha from tuning.yaml
  python scripts/tuning/apply_tuning.py --root R --alpha-aa A --alpha-tdi B         # explicit alpha
  python scripts/tuning/apply_tuning.py --root R --alpha-scale 1.5 \
      --out-dir data/interim/pairs_variants/alpha_x1.5                              # sensitivity variant
  python scripts/tuning/apply_tuning.py --root R --alpha-scale 0.6666667 --out-dir data/interim/pairs_variants/alpha_div1.5

For every alignment (default tdiaa seqonly tmalign) and every shard of ROOT/data/interim/<--pairs-subdir>/{aln}/
(all splits; the transform is deterministic and reads no outcome): aa_*/tdi_* (core + halves) at the tuned
alpha (x --alpha-scale), lddt_D / lddt_sat with the C28 S0 table, lddt_D_se from the §2.3 SE model.
aa_ml_seq / aa_pid_seq are kept at alpha = 1 (C26). boot_sub = pair_id in
ROOT/data/interim/pairs_bootsub/bootsub_pairs.tsv (if present). Sources, first found wins:
  alpha    : --alpha-aa/--alpha-tdi > ROOT/config/tuning.yaml alpha > ROOT/results/tuning/alpha.json
  S0       : tuning.yaml s0.{aln} > results/tuning/floors_tuned.json > results/floors/floors.json (C28 built)
  SE model : tuning.yaml lddt_se_model.{aln} > results/tuning/se_calibration.json (else lddt_D_se stays NaN)
In place (default) updates the .done md5; --out-dir (ROOT-relative or absolute) writes DIR/{aln}/shard_k.parquet,
counts_shard_k.npz (copy) and .done, plus DIR/apply_tuning.json (log; never written inside a {aln} dir).
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import time
from pathlib import Path

from se_common import ALNS, PROJECT, add_root_arg, out_json, read_alpha, read_s0, write_json

spec = importlib.util.spec_from_file_location("recompute_ml", PROJECT / "scripts" / "pairs" / "recompute_ml.py")
rml = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rml)


def se_model_for(root, aln):
    import yaml
    y = Path(root) / "config" / "tuning.yaml"
    if y.exists():
        m = ((yaml.safe_load(open(y)) or {}).get("lddt_se_model") or {}).get(aln) or {}
        if m.get("sd_S") is not None:
            return m, str(y)
    j = out_json(root, "se_calibration")
    if j.exists():
        m = (json.load(open(j)).get("lddt_se_model") or {}).get(aln)
        if m and m.get("sd_S") is not None:
            return m, str(j)
    return None, None


def main():
    ap = add_root_arg(argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter))
    ap.add_argument("--alpha-aa", type=float)
    ap.add_argument("--alpha-tdi", type=float)
    ap.add_argument("--alpha-scale", type=float, default=1.0, help="multiply both alphas (sensitivity variants)")
    ap.add_argument("--out-dir", help="variant dir (ROOT-relative or absolute); default: in place")
    ap.add_argument("--alns", nargs="+", default=list(ALNS))
    ap.add_argument("--pairs-subdir", default="pairs", help="data/interim/<subdir> (e.g. pairs_floors)")
    ap.add_argument("--device", default=None)
    a = ap.parse_args()
    t0 = time.time()
    root = Path(a.root)
    alpha0, asrc = read_alpha(root, a.alpha_aa, a.alpha_tdi)
    alpha = {k: v * a.alpha_scale for k, v in alpha0.items()}
    src = root / "data" / "interim" / a.pairs_subdir
    dst = None
    if a.out_dir:
        dst = Path(a.out_dir) if Path(a.out_dir).is_absolute() else root / a.out_dir
        if dst.resolve() == src.resolve():
            raise SystemExit("--out-dir equals the source pairs dir; omit it for in place")
        dst.mkdir(parents=True, exist_ok=True)
    log = {"alpha": alpha, "alpha_tuned": alpha0, "alpha_scale": a.alpha_scale, "alpha_source": asrc,
           "source": str(src), "dest": str(dst) if dst else "in place", "selection": "kept at alpha=1 (C26)",
           "alignments": {}}
    bs_path = root / "data" / "interim" / "pairs_bootsub" / "bootsub_pairs.tsv"
    boot = None
    if bs_path.exists():
        import polars as pl
        b = pl.read_csv(bs_path, separator="\t", infer_schema_length=0)
        boot = set((b["gene_a"] + "__" + b["gene_b"]).to_list())
    log["boot_sub_list"] = str(bs_path) if boot is not None else None
    for aln in a.alns:
        s0, s0src = read_s0(root, aln)
        se, sesrc = se_model_for(root, aln)
        files = rml.shard_files(str(src), aln)
        n = 0
        for pq in files:
            n += rml.recompute_shard(pq, aln, alpha, str(dst) if dst else None, s0, se, False, 1.0, a.device,
                                     boot_sub=boot)["rows"]
        log["alignments"][aln] = {"shards": len(files), "rows": n, "s0_source": s0src, "se_model_source": sesrc}
        print(f"{aln}: {len(files)} shards, {n} rows; alpha {alpha}; S0 {s0src}; SE model {sesrc}", flush=True)
    log["seconds"] = round(time.time() - t0, 1)
    logp = (dst if dst else src) / "apply_tuning.json"
    write_json(log, logp)
    print(logp)


if __name__ == "__main__":
    main()
