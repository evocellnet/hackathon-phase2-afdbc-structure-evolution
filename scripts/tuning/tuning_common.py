"""Shared helpers for the §3 tuning scripts (scripts/tuning/*.py).

Every tuning script takes --root (default: the project root) and resolves ALL paths relative to it; it reads only
split == 'tune' through strucevo.stats.common.load_pairs(..., 'tune', root=ROOT).

Tuned-parameter fragments written by the SE-calibration / floors agent are read when present (else provisional
values are used and marked `provisional`). Expected keys (same nesting as config/tuning.yaml; extra keys ignored):
  results/tuning/alpha.json          {"alpha": {"aa": float, "tdi": float}, ...}
  results/tuning/se_calibration.json {"kappa": {aln: {"bin_edges": [..], "aa_ml": [..], "tdi_ml": [..]}},
                                      "rho":   {aln: {"bin_edges": [..], "aa_tdi": [..], "aa_lddt": [..]}},
                                      "lddt_se_model": {aln: {...}}}
  results/tuning/floors_tuned.json   {"eta": {aln: {"tdi_ml": float, "lddt_D": float, "lddt_raw": float}},
                                      "s0": {aln: {...}}}
Lenient readers also accept the same blocks nested one level deeper under "values"/"tuning", and a top-level
alignment key ({aln: {"kappa": ...}}).

Provisional values (until the fragments exist): kappa = 1, rho = 0.3, eta_tdi = floors.json[aln].eta.tdi_ml.median
(nu = eta/2), tip SD = robust SD of floor_ident tdi_ml / sqrt(2) (an upper bound: includes estimation noise).
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS", "POLARS_MAX_THREADS"):
    os.environ.setdefault(_v, "2")

import json
import sys
import time
from pathlib import Path

import numpy as np

PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT / "src"))

from strucevo.stats.common import GLOBAL_SEED, load_pairs, md5_file, write_json  # noqa: E402
from strucevo.stats.trees import PanelTree  # noqa: E402

try:
    import torch
    torch.set_num_threads(2)
except Exception:  # pragma: no cover
    torch = None

# tuning-only random streams: SeedSequence(20261030, spawn_key=(1000 + k, ...)), disjoint from config/seeds.yaml
# stage indices (0..15); keys fixed here, never change.
TUNING_STREAMS = {"tau_misfit_sim": 1000, "b_e4_world": 1001, "e5_world": 1002, "v3_tuned": 1003, "v3_beta": 1004}


def tuning_rng(name: str, *extra: int) -> np.random.Generator:
    return np.random.default_rng(np.random.SeedSequence(GLOBAL_SEED, spawn_key=(TUNING_STREAMS[name], *extra)))


def add_common_args(ap):
    ap.add_argument("--root", default=str(PROJECT), help="project root (all paths resolved relative to it)")
    ap.add_argument("--aln", default="tdiaa", help="alignment (primary = tdiaa, D22)")
    return ap


def out_path(root, step) -> Path:
    p = Path(root) / "results" / "tuning" / f"{step}.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    return p


def split_md5(root):
    return md5_file(Path(root) / "data" / "interim" / "split.tsv")


def header(root, step, args) -> dict:
    return {"step": step, "root": str(root), "split_md5": split_md5(root), "split_used": "tune",
            "date": time.strftime("%Y-%m-%d %H:%M:%S %Z"), "args": {k: v for k, v in vars(args).items()},
            "threads": {"omp": os.environ.get("OMP_NUM_THREADS"),
                        "torch": (torch.get_num_threads() if torch is not None else None)}}


# ------------------------------------------------------------------------------------------ panel tree
def load_panel(root) -> tuple[PanelTree, dict]:
    """Panel tree from ROOT/data/interim/timetree/panel_tree.nwk with ages (median, lo95, hi95) from node_ages.tsv.
    Leaf labels 'UPxxxx_<taxid>' -> taxid. Undated internal nodes (median NA) are interpolated between the parent
    age and the oldest dated descendant (evenly along the longest undated chain) and listed in info['age_filled'];
    nodes without a CI keep lo/hi = NaN, i.e. they are not redrawn in the E4 age bootstrap."""
    import dendropy
    import pandas as pd
    tdir = Path(root) / "data" / "interim" / "timetree"
    t = dendropy.Tree.get(path=str(tdir / "panel_tree.nwk"), schema="newick", preserve_underscores=True,
                          suppress_internal_node_taxa=True, suppress_leaf_node_taxa=False)
    for lf in t.leaf_node_iter():
        lf.taxon.label = lf.taxon.label.rsplit("_", 1)[-1]
    na = pd.read_csv(tdir / "node_ages.tsv", sep="\t")
    na = na[~na["is_leaf"].astype(bool)]
    ages = {str(r.node): (float(r.median), float(r.lo95), float(r.hi95)) for r in na.itertuples()}
    panel = PanelTree.from_dendropy(t, ages)
    age = panel.age.copy()
    N = panel.n_nodes
    base = np.zeros(N); chain = np.zeros(N, int)
    for i in range(N - 1, -1, -1):
        if panel.is_leaf[i]:
            continue
        b, k = 0.0, 0
        for c in panel.children[i]:
            if np.isfinite(age[c]):
                b = max(b, age[c])
            else:
                b = max(b, base[c]); k = max(k, chain[c])
        base[i] = b; chain[i] = k + 1
    filled = []
    for i in range(N):
        if not panel.is_leaf[i] and not np.isfinite(age[i]):
            p = panel.parent[i]
            pa = age[p] if p >= 0 else base[i] * 1.1 + 1
            age[i] = base[i] + (pa - base[i]) * chain[i] / (chain[i] + 1)
            filled.append(panel.label[i])
    panel = panel.with_ages(age)
    info = {"n_nodes": N, "n_leaves": int(panel.is_leaf.sum()), "age_filled": filled,
            "n_no_ci": int(np.sum(~panel.is_leaf & ~(np.isfinite(panel.lo) & np.isfinite(panel.hi) & (panel.lo > 0))))}
    return panel, info


# ------------------------------------------------------------------------------------------ tuned params
def _read(root, name):
    p = Path(root) / "results" / "tuning" / f"{name}.json"
    if p.exists():
        try:
            return json.loads(p.read_text()), str(p.relative_to(root))
        except Exception:
            return None, None
    return None, None


def dig(d, key, aln=None):
    """Find block `key` (optionally its [aln] sub-block) in a fragment with lenient nesting."""
    if not isinstance(d, dict):
        return None
    cands = [d, d.get("values", {}), d.get("tuning", {})]
    if aln is not None and isinstance(d.get(aln), dict):
        cands.append(d[aln])
    for c in cands:
        if isinstance(c, dict) and key in c:
            v = c[key]
            if aln is not None and isinstance(v, dict) and aln in v:
                return v[aln]
            return v
    return None


class BinLookup:
    """Piecewise-constant lookup by aa_ml bin: bin_edges (len k+1) + values (len k); scalar allowed."""

    def __init__(self, edges, values, default):
        self.scalar = None
        if values is None or (edges is None and not np.isscalar(values)):
            self.scalar = float(default)
        elif np.isscalar(values):
            self.scalar = float(values)
        else:
            v = np.array([np.nan if x is None else float(x) for x in values])
            v = np.where(np.isfinite(v), v, np.nanmedian(v) if np.isfinite(v).any() else default)
            self.edges = np.asarray([np.inf if x is None else float(x) for x in edges])
            self.values = v

    def __call__(self, x):
        x = np.asarray(x, float)
        if self.scalar is not None:
            return np.full(x.shape, self.scalar)
        k = np.clip(np.searchsorted(self.edges, np.where(np.isfinite(x), x, 0.0), side="right") - 1, 0,
                    len(self.values) - 1)
        return self.values[k]

    def summary(self):
        return self.scalar if self.scalar is not None else float(np.median(self.values))


def load_tuned(root, aln="tdiaa") -> dict:
    """Tuned noise parameters for `aln` (fragments if present, else provisional)."""
    out = {"provisional": [], "sources": {}}
    se, se_src = _read(root, "se_calibration")
    fl, fl_src = _read(root, "floors_tuned")
    al, al_src = _read(root, "alpha")
    k = dig(se, "kappa", aln) if se else None
    if isinstance(k, dict) and (k.get("aa_ml") is not None or k.get("tdi_ml") is not None):
        out["kappa_aa"] = BinLookup(k.get("bin_edges"), k.get("aa_ml"), 1.0)
        out["kappa_tdi"] = BinLookup(k.get("bin_edges"), k.get("tdi_ml"), 1.0)
        out["sources"]["kappa"] = se_src
    else:
        out["kappa_aa"] = BinLookup(None, None, 1.0); out["kappa_tdi"] = BinLookup(None, None, 1.0)
        out["provisional"].append("kappa=1"); out["sources"]["kappa"] = "provisional"
    r = dig(se, "rho", aln) if se else None
    if isinstance(r, dict) and r.get("aa_tdi") is not None:
        out["rho_tdi"] = BinLookup(r.get("bin_edges"), r.get("aa_tdi"), 0.3)
        out["sources"]["rho"] = se_src
    else:
        out["rho_tdi"] = BinLookup(None, None, 0.3)
        out["provisional"].append("rho=0.3"); out["sources"]["rho"] = "provisional"
    e = dig(fl, "eta", aln) if fl else None
    floors = json.loads((Path(root) / "results" / "floors" / "floors.json").read_text())
    feta = floors.get(aln, {}).get("eta", {}).get("tdi_ml", {})
    if isinstance(e, dict) and e.get("tdi_ml") is not None:
        out["eta_tdi"] = float(e["tdi_ml"]); out["sources"]["eta"] = fl_src
    else:
        out["eta_tdi"] = float(feta.get("median", 0.04))
        out["provisional"].append("eta_tdi=floors.json median (untuned alpha)")
        out["sources"]["eta"] = "results/floors/floors.json (provisional)"
    out["nu_tdi"] = out["eta_tdi"] / 2
    # per-protein tip noise SD: robust SD of floor_ident tdi_ml (q05..q95 range / 3.29) / sqrt(2)
    if feta.get("q95") is not None:
        out["tip_sd"] = float((feta["q95"] - feta["q05"]) / 3.29 / np.sqrt(2))
        out["sources"]["tip_sd"] = "floors.json floor_ident tdi_ml (q95-q05)/3.29/sqrt2 (upper bound)"
    else:
        out["tip_sd"] = 0.01; out["provisional"].append("tip_sd=0.01")
    a = dig(al, "alpha") if al else None
    out["alpha"] = a if isinstance(a, dict) else None
    out["sources"]["alpha"] = al_src or "absent (pairs computed at alpha=1)"
    out["summary"] = {"kappa_aa": out["kappa_aa"].summary(), "kappa_tdi": out["kappa_tdi"].summary(),
                      "rho_tdi": out["rho_tdi"].summary(), "eta_tdi": out["eta_tdi"], "nu_tdi": out["nu_tdi"],
                      "tip_sd": out["tip_sd"], "alpha": out["alpha"]}
    return out


# ------------------------------------------------------------------------------------------ tune data
def hog_species_tune(root) -> dict:
    import polars as pl
    h = pl.read_parquet(Path(root) / "data" / "interim" / "hogs.parquet", columns=["hog", "split", "in_sample", "selected"])
    h = h.filter((pl.col("split") == "tune") & pl.col("in_sample"))
    return {r[0]: np.asarray(r[1], dtype=np.int64) for r in h.select(["hog", "selected"]).iter_rows()}


def load_tune_primary(root, aln, tuned) -> tuple:
    """Primary pairs of split=='tune' with calibrated SEs (kappa by aa_ml bin applied to aa_ml_se/tdi_ml_se).
    Returns (pandas-free dict table, hog_species, metrics)."""
    df = load_pairs(aln, "tune", root=root, pair_set="primary")
    tab = {c: df.get_column(c).to_numpy() for c in df.columns}
    x = tab["aa_ml"].astype(float)
    tab["aa_ml_se_fisher"] = tab["aa_ml_se"].astype(float)
    tab["tdi_ml_se_fisher"] = tab["tdi_ml_se"].astype(float)
    tab["aa_ml_se"] = tab["aa_ml_se_fisher"] * tuned["kappa_aa"](x)
    tab["tdi_ml_se"] = tab["tdi_ml_se_fisher"] * tuned["kappa_tdi"](x)
    for c in ("aa_sat", "tdi_sat", "lddt_sat"):
        if c in tab:
            v = tab[c]
            tab[c] = np.array([bool(z) if z is not None else False for z in v]) if v.dtype == object else v.astype(bool)
    hs = hog_species_tune(root)
    # S_h must contain every species that appears in the pairs
    hog = tab["hog"].astype(str)
    for h in np.unique(hog):
        m = hog == h
        sp = np.unique(np.r_[tab["species_a"][m], tab["species_b"][m]]).astype(np.int64)
        if h not in hs or not np.all(np.isin(sp, hs[h])):
            hs[h] = np.union1d(hs.get(h, np.zeros(0, np.int64)), sp)
    metrics = ["aa_ml", "tdi_ml"]
    if "lddt_D" in tab and np.isfinite(tab["lddt_D"].astype(float)).any():
        metrics.append("lddt_D")
    return tab, hs, metrics


def fit_tune_branches(root, aln, tuned, tau_misfit=None, metrics=None, misfit=True):
    from strucevo.stats.branch import fit_branches
    panel, pinfo = load_panel(root)
    tab, hs, mets = load_tune_primary(root, aln, tuned)
    mets = metrics or mets
    nu = {"tdi_ml": tuned["nu_tdi"]}
    t = time.time()
    bf = fit_branches(tab, panel, metrics=mets, nu=nu, aln=aln, hog_species={h: hs[h] for h in np.unique(tab["hog"].astype(str))},
                      tau_misfit=tau_misfit, misfit=misfit, se_method="model")
    return bf, panel, tab, {"panel": pinfo, "branch_fit_s": time.time() - t, "metrics": mets,
                            "n_hogs": len(bf.hogs), "n_pairs": int(len(tab["hog"]))}


def tau_from_fragment(root, default=None):
    d, _ = _read(root, "tau_misfit")
    if d and d.get("tau_misfit") is not None:
        return float(d["tau_misfit"])
    return default


def observe_world_tuned(w, tuned, rng, **kw):
    """observe_gauss with the tuned noise (scalar summaries of kappa / rho)."""
    from strucevo.stats.synth import observe_gauss
    k = tuned["summary"]["kappa_tdi"]
    t = observe_gauss(w, rng=rng, kappa=k, rho=tuned["summary"]["rho_tdi"], nu_tdi=tuned["nu_tdi"],
                      tip_sd=tuned["tip_sd"], **kw)
    # observe_gauss returns Fisher-like SEs (true noise = kappa x SE); return calibrated SEs, as the branch engine
    # receives on real data (load_tune_primary applies kappa)
    t["aa_ml_se"] = t["aa_ml_se"] * k
    t["tdi_ml_se"] = t["tdi_ml_se"] * k
    return t


def dump(rec, path):
    write_json(rec, path)
    print(f"wrote {path}")
    return path
