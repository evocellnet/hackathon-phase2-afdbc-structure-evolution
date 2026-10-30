"""Merge the tuning fragments ROOT/results/tuning/*.json into ROOT/config/tuning.yaml (ANALYSIS_PLAN §3.4).

The existing tuning.yaml is used as the template: every line, key order and comment is kept; only the VALUE text of
known leaves is replaced (flow style). Keys are never renamed; two keys are appended (allowed by the schema rule):
`sources` (yaml path -> fragment file + date, or 'computed') and `provisional` (paths whose fragment was computed with
provisional noise parameters; such values are written as null unless --allow-provisional).

Fragment -> yaml mapping (fragment keys expected; anything missing stays null):
  alpha.json            alpha.{aa,tdi}
  floors_tuned.json     eta.{aln}.{tdi_ml,lddt_D,lddt_raw}; s0.{aln}.{bin_lo,median,q95,source}
  se_calibration.json   kappa.{aln}.{bin_edges,aa_ml,tdi_ml}; rho.{aln}.{bin_edges,aa_tdi,aa_lddt};
                        lddt_se_model.{aln}.{len_bin_lo,s_bin_lo,sd_S}
  tau_misfit.json       tau_misfit
  b_e4_timing.json      B_E4
  e5_power.json         yaml.{power_f1.5,power_f1.25,power_f2,null_rate_f1,p_accept_f1.5,min_t_snr} -> e5.*
  v3_tuned.json         yaml.{false_concave_rate,min_detectable_c,min_detectable_beta_e4} -> v3.*
  sample_throughput.json sample.{hogs_per_stratum,n_primary_pairs,n_e6_pairs};
                        throughput.{aln}.pairs_per_s_48threads -> throughput_pairs_per_s.{aln}
  split_md5             md5 of ROOT/data/interim/split.tsv (computed here; every fragment's split_md5 must match)
Prints the list of null leaves at the end (exit 0; --strict -> exit 1 if any null).

  python scripts/tuning/write_tuning_yaml.py --root R [--allow-provisional] [--out PATH] [--strict]
"""
import argparse
import json
import re
import sys
import time
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
PROJECT = HERE.parents[1]
sys.path.insert(0, str(PROJECT / "src"))
from strucevo.stats.common import md5_file  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--root", default=str(PROJECT))
ap.add_argument("--template", default=None, help="default: ROOT/config/tuning.yaml")
ap.add_argument("--out", default=None, help="default: ROOT/config/tuning.yaml")
ap.add_argument("--allow-provisional", action="store_true")
ap.add_argument("--strict", action="store_true")
a = ap.parse_args()
root = Path(a.root)
tpl_path = Path(a.template) if a.template else root / "config" / "tuning.yaml"
out_path = Path(a.out) if a.out else root / "config" / "tuning.yaml"
tdir = root / "results" / "tuning"
lines = tpl_path.read_text().splitlines()
# strip previously appended blocks (idempotent re-runs)
for i, l in enumerate(lines):
    if l.startswith("# ---- appended by write_tuning_yaml.py"):
        lines = lines[:i]
        break
while lines and not lines[-1].strip():
    lines.pop()
template = yaml.safe_load("\n".join(lines))


def frag(name):
    p = tdir / f"{name}.json"
    if not p.exists():
        return None, None
    d = json.loads(p.read_text())
    date = time.strftime("%Y-%m-%d", time.localtime(p.stat().st_mtime))
    return d, f"results/tuning/{name}.json ({date})"


def is_prov(d):
    if not isinstance(d, dict):
        return []
    n = d.get("noise", {})
    return list(n.get("provisional", [])) if isinstance(n, dict) else []


values, sources, provisional, mismatch = {}, {}, {}, []
smd5 = md5_file(root / "data" / "interim" / "split.tsv")


def put(path, v, src, prov=()):
    if v is None:
        return
    values[path] = v
    sources[path] = src
    if prov:
        provisional[path] = list(prov)


def leaves(d, pre=()):
    for k, v in d.items():
        if isinstance(v, dict) and v:
            yield from leaves(v, pre + (str(k),))
        else:
            yield pre + (str(k),), v


put(("split_md5",), smd5, "computed: md5(data/interim/split.tsv)")
# generic blocks whose nesting equals the yaml
for name, blocks in (("alpha", ("alpha",)), ("floors_tuned", ("eta", "s0")),
                     ("se_calibration", ("kappa", "rho", "lddt_se_model"))):
    d, src = frag(name)
    if d is None:
        continue
    if d.get("split_md5") not in (None, smd5):
        mismatch.append(name)
    for b in blocks:
        blk = d.get(b)
        if not isinstance(blk, dict):
            continue
        for p, v in leaves(blk, (b,)):
            put(p, v, src, is_prov(d))
d, src = frag("tau_misfit")
if d:
    mismatch += ["tau_misfit"] if d.get("split_md5") not in (None, smd5) else []
    put(("tau_misfit",), d.get("tau_misfit"), src, is_prov(d))
d, src = frag("b_e4_timing")
if d:
    mismatch += ["b_e4_timing"] if d.get("split_md5") not in (None, smd5) else []
    put(("B_E4",), d.get("B_E4"), src)          # timing does not depend on the noise values
for name, blk in (("e5_power", "e5"), ("v3_tuned", "v3")):
    d, src = frag(name)
    if d and isinstance(d.get("yaml"), dict):
        mismatch += [name] if d.get("split_md5") not in (None, smd5) else []
        extra = ""
        if name == "e5_power" and d.get("power_source"):
            extra = f" [{d['power_source']}]"
        for k, v in d["yaml"].items():
            if (blk, k) == ("e5", "min_t_snr"):
                continue          # fixed setting (D8), not a tuned value: the template value is kept
            put((blk, k), v, src + extra, is_prov(d))
d, src = frag("sample_throughput")
if d:
    s = d.get("sample", {})
    hps = s.get("hogs_per_stratum")
    put(("sample", "hogs_per_stratum"), {int(k): v for k, v in hps.items()} if isinstance(hps, dict) else hps, src)
    put(("sample", "n_primary_pairs"), s.get("n_primary_pairs"), src)
    put(("sample", "n_e6_pairs"), s.get("n_e6_pairs"), src)
    for aln, t in d.get("throughput", {}).items():
        if t.get("pairs_per_s_48threads") is not None:
            put(("throughput_pairs_per_s", aln), round(float(t["pairs_per_s_48threads"]), 1),
                src + f" (measured at {t.get('threads_used')} threads, scaled to 48)")

if not a.allow_provisional:
    for p in list(provisional):
        values[p] = None          # explicit null (overrides any value already in the template)
        sources[p] = sources[p] + " -- PROVISIONAL noise, left null (use --allow-provisional)"


def fmt(v):
    if isinstance(v, float):
        return repr(round(v, 6)) if abs(v) >= 1e-4 or v == 0 else repr(v)
    s = yaml.safe_dump(v, default_flow_style=True, width=10 ** 9, sort_keys=False).strip()
    if s.endswith("\n..."):
        s = s[:-4].strip()
    if s.endswith("..."):
        s = s[:-3].strip()
    return s


def round_floats(v):
    if isinstance(v, float):
        return round(v, 6) if abs(v) >= 1e-4 or v == 0 else v
    if isinstance(v, list):
        return [round_floats(x) for x in v]
    if isinstance(v, dict):
        return {k: round_floats(x) for k, x in v.items()}
    return v


# ---- render: walk template lines, replace leaf values (or whole flow maps) in place
pat = re.compile(r"^(\s*)([A-Za-z0-9_.]+):(\s*)([^#]*?)(\s*#.*)?$")
stack = []          # (indent, key)
out = []
used = set()
for l in lines:
    m = pat.match(l)
    if not m or l.lstrip().startswith("#"):
        out.append(l); continue
    ind, key, sp, val, com = len(m.group(1)), m.group(2), m.group(3), m.group(4), m.group(5) or ""
    while stack and stack[-1][0] >= ind:
        stack.pop()
    path = tuple(k for _, k in stack) + (key,)
    if val.strip() == "":
        stack.append((ind, key)); out.append(l); continue
    tv = val.strip()
    if path in values:                           # whole value (scalar, list or map) from a fragment
        newv = fmt(round_floats(values[path])); used.add(path)
    elif tv.startswith("{"):                     # flow map: merge sub-leaves
        sub = yaml.safe_load(tv) or {}
        changed = False
        for k in list(sub) + [p[-1] for p in values if p[:-1] == path and p[-1] not in sub]:
            if path + (k,) in values:
                sub[k] = round_floats(values[path + (k,)]); used.add(path + (k,)); changed = True
        newv = fmt(sub) if changed else tv
    elif path in values:
        newv = fmt(round_floats(values[path])); used.add(path)
    else:
        newv = tv
    pad = sp if sp else " "
    if com:
        # keep the comment column where possible
        base = f"{m.group(1)}{key}:{pad}{newv}"
        col = len(l) - len(com.lstrip())          # original comment column
        gap = col - len(base) if len(base) < col else 1
        out.append(base + " " * gap + com.lstrip())
    else:
        out.append(f"{m.group(1)}{key}:{pad}{newv}")
unused = [p for p in values if p not in used]
out.append("")
out.append("# ---- appended by write_tuning_yaml.py (new keys; schema keys above are unchanged) ----")
out.append(f"written: {fmt(time.strftime('%Y-%m-%d %H:%M:%S %Z'))}")
out.append("sources:                        # yaml path -> fragment (date); PROVISIONAL = computed with untuned noise")
for p in sorted(sources):
    out.append(f"  {fmt('.'.join(p))}: {fmt(sources[p])}")
out.append("provisional:                    # paths whose fragment used provisional kappa/rho/eta")
if provisional:
    for p in sorted(provisional):
        out.append(f"  {fmt('.'.join(p))}: {fmt(provisional[p])}")
else:
    out[-1] = "provisional: {}                 # paths whose fragment used provisional kappa/rho/eta"
text = "\n".join(out) + "\n"
parsed = yaml.safe_load(text)           # must parse
# check: every merged value round-trips
bad = []
for p in used:
    cur = parsed
    for k in p:
        cur = cur[k]
    if json.dumps(round_floats(values[p]), default=str) != json.dumps(cur, default=str):
        bad.append(".".join(p))
if bad:
    sys.exit(f"round-trip mismatch for {bad}")
out_path.parent.mkdir(parents=True, exist_ok=True)
out_path.write_text(text)
nulls = [".".join(p) for p, v in leaves({k: v for k, v in parsed.items() if k not in ("sources", "provisional")}) if v is None]
print(f"wrote {out_path}  (split_md5 {smd5})")
if mismatch:
    print("WARNING split_md5 mismatch in fragments:", mismatch)
if unused:
    print("fragment values with no slot in the schema (not written):", [".".join(p) for p in unused][:20])
if provisional:
    print(("provisional values WRITTEN: " if a.allow_provisional else "provisional values left null: ")
          + ", ".join(".".join(p) for p in sorted(provisional)))
print(f"NULLS ({len(nulls)}):")
for n in nulls:
    print("  ", n)
sys.exit(1 if (a.strict and nulls) else 0)
