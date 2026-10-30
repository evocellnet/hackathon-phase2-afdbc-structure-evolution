"""Fill the '## §Tuning' table of ANALYSIS_PLAN.md from config/tuning.yaml (ANALYSIS_PLAN §3.4).

Replaces ONLY the Value and Source/date cells of the table rows below the '## §Tuning' heading (rows are matched by
their Parameter cell; a row whose yaml values are all missing is left untouched). Every other line of the file is
verified byte-identical before writing (the script aborts otherwise) and the changed lines are printed as a diff.

  python scripts/tuning/fill_tuning_table.py --plan ANALYSIS_PLAN.md --tuning config/tuning.yaml [--dry-run]
"""
import argparse
import difflib
import sys
import time
from pathlib import Path

import yaml

ap = argparse.ArgumentParser()
ap.add_argument("--plan", required=True)
ap.add_argument("--tuning", required=True)
ap.add_argument("--aln", default="tdiaa", help="primary alignment shown in the table (D22)")
ap.add_argument("--dry-run", action="store_true")
a = ap.parse_args()
plan = Path(a.plan)
text = plan.read_text()
lines = text.split("\n")
Y = yaml.safe_load(Path(a.tuning).read_text())
SRC = Y.get("sources", {}) or {}
PROV = Y.get("provisional", {}) or {}
aln = a.aln
today = time.strftime("%Y-%m-%d")


def g(*path):
    cur = Y
    for k in path:
        if not isinstance(cur, dict) or k not in cur:
            return None
        cur = cur[k]
    return cur


def f(x, nd=3):
    if x is None:
        return "n/a"
    if isinstance(x, float):
        if x.is_integer():
            return str(int(x))
        return f"{x:.{nd}g}" if 1e-3 <= abs(x) < 10 ** nd else (f"{x:.0f}" if abs(x) >= 10 ** nd else f"{x:.2e}")
    return str(x)


def src(*paths):
    s = sorted({str(SRC.get(".".join(p), "")).split(" (")[0] for p in paths if SRC.get(".".join(p))})
    prov = any(".".join(p) in PROV for p in paths)
    files = ", ".join(x.replace("results/tuning/", "") for x in s if x) or "config/tuning.yaml"
    return f"{files}{' PROVISIONAL' if prov else ''}; {today}"


def rng(v, nd=3):
    v = [x for x in (v or []) if x is not None]
    if not v:
        return None
    return f"{min(v):.{nd}g}–{max(v):.{nd}g}"


def row_alpha_aa():
    v = g("alpha", "aa")
    return (f(v), src(("alpha", "aa"))) if v is not None else None


def row_alpha_tdi():
    v = g("alpha", "tdi")
    return (f(v), src(("alpha", "tdi"))) if v is not None else None


def row_s0():
    s = g("s0", aln)
    if not s or s.get("median") is None:
        return None
    med = ", ".join(f"{x:.3f}" for x in s["median"]); q = ", ".join(f"{x:.3f}" for x in s["q95"])
    bins = ", ".join(str(b) for b in s["bin_lo"])
    return (f"{aln} lddt_D floor (C28 {'/'.join(sorted(set(s.get('source') or [])))}), len bins ≥ [{bins}]: "
            f"median [{med}]; q95 [{q}] (all alignments in tuning.yaml#s0)", src(("s0", aln, "median")))


def row_tau():
    v = g("tau_misfit")
    return (f(v), src(("tau_misfit",))) if v is not None else None


def row_eta():
    e = g("eta", aln) or {}
    if e.get("tdi_ml") is None and e.get("lddt_D") is None:
        return None
    return (f"{aln}: η_tdi = {f(e.get('tdi_ml'))}, η_lddtD = {f(e.get('lddt_D'))} (ν = η/2)",
            src(("eta", aln, "tdi_ml"), ("eta", aln, "lddt_D")))


def row_kappa():
    k = g("kappa", aln) or {}
    r = g("rho", aln) or {}
    if k.get("aa_ml") is None and r.get("aa_tdi") is None:
        return None
    parts = []
    if k.get("aa_ml") is not None:
        parts.append(f"κ_AA {rng(k['aa_ml'])}, κ_3Di {rng(k.get('tdi_ml'))}")
    if r.get("aa_tdi") is not None:
        parts.append(f"ρ_AA,3Di {rng(r['aa_tdi'])}" + (f", ρ_AA,lddt {rng(r.get('aa_lddt'))}" if r.get("aa_lddt") else ""))
    return (f"{aln}, range over aa_ml bins: " + "; ".join(parts) + " (per bin in tuning.yaml)",
            src(("kappa", aln, "aa_ml"), ("rho", aln, "aa_tdi")))


def row_b():
    v = g("B_E4")
    return (f(v), src(("B_E4",))) if v is not None else None


def row_sample():
    h = g("sample", "hogs_per_stratum"); p = g("sample", "n_primary_pairs"); e = g("sample", "n_e6_pairs")
    if h is None and p is None:
        return None
    hs = ", ".join(f"{k}: {v}" for k, v in sorted((h or {}).items(), key=lambda kv: str(kv[0])))
    tot = (p or 0) + (e or 0)
    return (f"strata {{{hs}}}; primary {p:,} + E6 {e:,} = {tot:,} pairs".replace(",", ",") if p is not None and e is not None
            else f"strata {{{hs}}}", src(("sample", "hogs_per_stratum"), ("sample", "n_primary_pairs")))


def row_thr():
    t = g("throughput_pairs_per_s") or {}
    if t.get("tmalign") is None:
        return None
    others = ", ".join(f"{k} {f(v, 4)}" for k, v in t.items() if k != "tmalign" and v is not None)
    s = SRC.get("throughput_pairs_per_s.tmalign", "")
    note = s[s.find("(measured"):].strip("()") if "(measured" in s else "48 threads"
    return (f"{f(t['tmalign'], 4)} (48 threads; {others}) [{note}]", src(("throughput_pairs_per_s", "tmalign")))


def row_e5():
    e = g("e5") or {}
    if e.get("power_f1.5") is None:
        return None
    s = SRC.get("e5.power_f1.5", "")
    basis = s[s.find("[") + 1:s.find("]")] if "[" in s else ""
    return (f"{f(e['power_f1.5'])} (f = 1.25: {f(e.get('power_f1.25'))}; f = 2: {f(e.get('power_f2'))}; "
            f"null rate f = 1: {f(e.get('null_rate_f1'))}; P(accept \\| f = 1.5) = {f(e.get('p_accept_f1.5'))}; "
            f"min_t_snr = {f(e.get('min_t_snr'))}{'; ' + basis if basis else ''})", src(("e5", "power_f1.5")))


def row_v3():
    v = g("v3") or {}
    if v.get("false_concave_rate") is None and v.get("min_detectable_c") is None:
        return None
    return (f"false-concave {f(v.get('false_concave_rate'))}; min. detectable \\|c−1\\| = {f(v.get('min_detectable_c'))}; "
            f"\\|β_E4\\| = {f(v.get('min_detectable_beta_e4'))}",
            src(("v3", "false_concave_rate"), ("v3", "min_detectable_c"), ("v3", "min_detectable_beta_e4")))


ROWS = [("α_AA", row_alpha_aa), ("α_3Di", row_alpha_tdi), ("S0(L)", row_s0), ("HOG-misfit threshold", row_tau),
        ("η_tdi", row_eta), ("κ_m(bin)", row_kappa), ("B_E4", row_b), ("HOGs per stratum", row_sample),
        ("tmalign pairs/s", row_thr), ("E5 spike-in power", row_e5), ("V3 false-concave", row_v3)]

# locate the table
try:
    h = next(i for i, l in enumerate(lines) if l.startswith("## §Tuning"))
except StopIteration:
    sys.exit("no '## §Tuning' heading")
i = h + 1
while i < len(lines) and not lines[i].startswith("|"):
    i += 1
start = i
while i < len(lines) and lines[i].startswith("|"):
    i += 1
end = i
new = list(lines)
changed, filled, missing = [], [], []
for j in range(start + 2, end):           # skip header + separator
    l = lines[j]
    # split on unescaped pipes
    cells, buf, k = [], "", 0
    while k < len(l):
        if l[k] == "\\" and k + 1 < len(l):
            buf += l[k:k + 2]; k += 2; continue
        if l[k] == "|":
            cells.append(buf); buf = ""
        else:
            buf += l[k]
        k += 1
    cells.append(buf)
    # cells: ['', ' Parameter ', ' Value ', ' Source ', '']
    if len(cells) != 5:
        sys.exit(f"unexpected table row: {l}")
    param = cells[1].strip()
    fn = next((fn for key, fn in ROWS if param.startswith(key)), None)
    if fn is None:
        missing.append(param); continue
    res = fn()
    if res is None:
        missing.append(param); continue
    val, srcs = res
    # keep the plan's annotation for tau (e.g. '(≥ 3)')
    nl = f"| {param} | {val} | {srcs} |"
    if nl != l:
        new[j] = nl; changed.append(j); filled.append(param)
# verification: every other line unchanged
assert len(new) == len(lines)
others = [k for k in range(len(lines)) if k not in changed]
assert all(new[k] == lines[k] for k in others), "a line outside the §Tuning rows changed"
assert all(start + 2 <= k < end for k in changed)
out = "\n".join(new)
diff = list(difflib.unified_diff(lines, new, "plan(before)", "plan(after)", lineterm="", n=0))
print("\n".join(diff))
print(f"\nrows filled: {len(filled)}; rows left unchanged (no value): {missing}")
print(f"lines changed: {len(changed)} (all inside the §Tuning table, lines {start + 3}-{end}); other {len(others)} lines identical")
if not a.dry_run:
    plan.write_text(out)
    # re-verify from disk
    back = plan.read_text().split("\n")
    assert len(back) == len(lines) and all(back[k] == lines[k] for k in others)
    print(f"wrote {plan}")
