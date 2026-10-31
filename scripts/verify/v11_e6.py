"""V11: independent E6 point estimates (§4.6). Reference curve f_h(x) = eta + a_h x^c, a_h = pass-2 profile a_h at the
E1 c-hat (own pass-2 weights). P1 choice: the plan's seeded draw (stage 'e6_p1_choice'), plus orientation-free variant."""
import sys
import numpy as np
from v11_common import *
from strucevo.stats.common import stage_rng   # seed helper only (to reproduce the seeded P1 choice)
import v11_curve as vc

c_list = [float(a) for a in sys.argv[1:]] or [endpoint("E1")["estimate"]]
eta = TUN["eta"][PRIMARY]["tdi_ml"]
d, cnt, _ = vc.prepare()
hog = d["hog"].to_numpy(); x = d["aa_ml"].to_numpy().astype(float); y = d["tdi_ml"].to_numpy().astype(float)
sy = d["tdi_ml_se_k"].to_numpy().astype(float)
F = vc.Fit(hog, x, y, eta)
c1, _ = F.chat(1 / sy ** 2); _, _, r1 = F.prof(c1, 1 / sy ** 2)
s2 = max(0.0, np.mean(r1 ** 2) - np.mean(sy ** 2)); w = 1 / (sy ** 2 + s2)
hog_ids = np.unique(hog)

P = pairs().filter(pl.col("pair_set").str.starts_with("e6_"))
assert set(P["split"].unique()) == {"confirm"}
P = P.filter(pl.col("status") == "ok")
look = {}
for i, (a, b) in enumerate(zip(P["gene_a"].to_list(), P["gene_b"].to_list())):
    look[(a, b)] = i; look[(b, a)] = i
XA, XB = P["aa_ml_hA"].to_numpy().astype(float), P["aa_ml_hB"].to_numpy().astype(float)
YB = P["tdi_ml_hB"].to_numpy().astype(float)
para = load_pairs(PRIMARY, SPLIT, pair_set="e6_para").select("gene_a", "gene_b", "aa_pid_seq", "aa_ml_seq")
pk = {}
for a, b, pid, ds in para.iter_rows():
    pk[(a, b)] = pk[(b, a)] = (pid, ds)
T = pl.read_csv(ROOT / "data/interim/pairlists/e6_triplets.tsv", separator="\t", infer_schema_length=0)
T = T.filter(pl.col("split") == SPLIT)
assert set(T["split"].unique()) == {"confirm"}
proj = endpoint("E6")
res = {"n_triplets_confirm": T.height}


def run(c, flipmode):
    _, a_h, _ = F.prof(c, w)
    ah = dict(zip(hog_ids, a_h))
    def r(g1, g2, h, xa):
        i = look.get((g1, g2))
        if i is None: return np.nan
        xv = (XA if xa == "A" else XB)[i]; yv = YB[i]
        if not (np.isfinite(xv) and np.isfinite(yv)): return np.nan
        return yv - (eta + ah[h] * max(xv, 0.0) ** c)
    def xv(g1, g2):
        i = look.get((g1, g2)); return XA[i] if i is not None else np.nan
    n = T.height
    if flipmode == "seeded":
        flip = stage_rng("e6_p1_choice", 0).random(n) < 0.5
    else:
        flip = np.zeros(n, bool)
    excl = dict(leaf_level=0, pid_or_dseq=0, no_E1_fit=0, no_O1=0, no_O2=0, missing_a=0, missing_b=0)
    da, db = {}, {}
    for i, row in enumerate(T.iter_rows(named=True)):
        h = row["hog"]
        if row["event_level"] == "leaf": excl["leaf_level"] += 1; continue
        pid, ds = pk.get((row["P1"], row["P2"]), (np.nan, np.nan))
        if not (pid <= 0.95 and ds <= 2): excl["pid_or_dseq"] += 1; continue
        if h not in ah or not np.isfinite(ah[h]): excl["no_E1_fit"] += 1; continue
        o1 = [g for g in (row["O1"] or "").split(",") if g]; o2 = [g for g in (row["O2"] or "").split(",") if g]
        assert not (set(o1) & set(o2))
        if not o1: excl["no_O1"] += 1; continue
        c1_, c2_ = row["P1"], row["P2"]
        def delta(p1, p2):
            rpp = r(p1, p2, h, "A"); ro = [v for v in (r(p1, o, h, "A") for o in o1) if np.isfinite(v)]
            return rpp - np.mean(ro) if (np.isfinite(rpp) and ro) else np.nan
        if flipmode == "both":
            v = np.nanmean([delta(c1_, c2_), delta(c2_, c1_)]) if np.isfinite([delta(c1_, c2_), delta(c2_, c1_)]).any() else np.nan
        else:
            p1, p2 = (c2_, c1_) if flip[i] else (c1_, c2_)
            v = delta(p1, p2)
        if np.isfinite(v): da.setdefault((h, row["dcg_id"]), []).append(v)
        else: excl["missing_a"] += 1
        if not o2: excl["no_O2"] += 1; continue
        m1 = np.nanmean([xv(c1_, o) for o in o1]); m2 = np.nanmean([xv(c2_, o) for o in o1])
        if not (np.isfinite(m1) and np.isfinite(m2)): excl["missing_b"] += 1; continue
        mdo, ldo = (c1_, c2_) if m1 > m2 else (c2_, c1_)
        ex = [v for v in (r(mdo, o, h, "B") - r(ldo, o, h, "B") for o in o2) if np.isfinite(v)]
        if ex: db.setdefault((h, row["dcg_id"]), []).append(np.mean(ex))
        else: excl["missing_b"] += 1
    va = np.array([np.mean(v) for v in da.values()]); vb = np.array([np.mean(v) for v in db.values()])
    nh = len(set(k[0] for k in da) | set(k[0] for k in db))
    return {"c": c, "E6a": float(va.mean()), "E6b": float(vb.mean()), "n_a": len(va), "n_b": len(vb), "n_hogs": nh, "excl": excl}


for c in c_list:
    for fm in ("seeded", "both"):
        o = run(c, fm)
        o["E6a_absdiff"] = abs(o["E6a"] - proj["E6a"]); o["E6b_absdiff"] = abs(o["E6b"] - proj["E6b"])
        res[f"c={c:.6f}|P1={fm}"] = o
        print(o, flush=True)
res["project"] = {"E6a": proj["E6a"], "E6b": proj["E6b"], "n_units": proj["n_units"], "n_hogs": proj["n_hogs"], "n_excluded": proj["n_excluded"]}
wjson(res, OUT / "V11_E6.json")
