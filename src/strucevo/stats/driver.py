"""Endpoint driver (ANALYSIS_PLAN §3 step 4-5, §4, §5): real-data glue between the pair/branch tables and the
strucevo.stats engines. Used by scripts/endpoints/run_confirm.py and scripts/endpoints/run_dryrun.py.

Every path is relative to `root` (default: project root), so the same code runs on test roots.

Mapping choices (documented here; PI decisions listed in the endpoint-engineer hand-off):
  * kappa_m(bin) (tuning.yaml kappa[aln], bins on the pair's own aa_ml) multiplies the Fisher SE columns
    aa_ml_se / tdi_ml_se of every pair BEFORE any engine sees them ("SE used everywhere = kappa x Fisher SE",
    §2.3). The engines are then called with kappa = 1 (kx = ky = 1); SIMEX therefore adds noise with the
    calibrated SE. lddt_D_se comes from the §2.3 lddt SE model (already bootstrap-calibrated): kappa = 1.
  * rho_m(bin) -> one scalar per endpoint = the pair-count-weighted mean of rho(bin) over the pairs that
    enter the fit (engines take a scalar rho).
  * Primary alignment = tuning.yaml primary_alignment (tdiaa, D22); seqonly and tmalign are sensitivities.
  * misfit_ratio / excl_misfit come from the primary alignment (D22), AA metric.
  * hogs-level columns are written to a companion table data/interim/hogs_endpoint_cols_{split}.parquet;
    hogs.parquet itself is never modified.
"""
from __future__ import annotations

import copy
import json
import os
import time
import traceback
from pathlib import Path

import numpy as np

from .common import (ROOT as PROJECT_ROOT, ConfirmSplitGuard, confirm_allowed, guard_split, load_pairs, md5_file,
                     stage_rng, write_json, endpoint_record, md5_arrays)

ALNS = ("tdiaa", "seqonly", "tmalign")
PAIR_COLS_SE = {"aa_ml": "aa_ml_se", "tdi_ml": "tdi_ml_se"}

# ---------------------------------------------------------------------------------------- bootstrap sizes
PLAN_SIZES = dict(B1=200, Bs=50, Bs_boot=10, B_naive=1000, B_sens=50, B_alt=200, B3=1000, B4=None, B4_sens=None,
                  B4_alt=200, B5=1000, P=100, B_cv=200, spike_reps=50, spike_B=200, B6=1000, neg_perms=20,
                  neg_B=200, B3_alt=200)
QUICK_SIZES = dict(B1=20, Bs=10, Bs_boot=3, B_naive=20, B_sens=8, B_alt=10, B3=50, B4=10, B4_sens=5, B4_alt=5,
                   B5=50, P=10, B_cv=20, spike_reps=3, spike_B=20, B6=100, neg_perms=20, neg_B=10, B3_alt=10)
PROVISIONAL_B_E4 = 50


def _now():
    return time.strftime("%Y-%m-%dT%H:%M:%S")


# ======================================================================================== tuning
def _nulls(o, path=""):
    out = []
    if o is None:
        return [path or "<root>"]
    if isinstance(o, dict):
        for k, v in o.items():
            out += _nulls(v, f"{path}.{k}" if path else str(k))
    elif isinstance(o, list):
        for i, v in enumerate(o):
            out += _nulls(v, f"{path}[{i}]")
    return out


def tuning_nulls(root) -> list:
    import yaml
    p = Path(root) / "config" / "tuning.yaml"
    if not p.exists():
        return ["<config/tuning.yaml missing>"]
    return _nulls(yaml.safe_load(p.read_text()))


def _floors(root):
    p = Path(root) / "results" / "floors" / "floors.json"
    return json.loads(p.read_text()) if p.exists() else None


def _s0_from_floors(fl, aln):
    """C28: S0_ok (status-ok random pairs) per length bin, fallback to unconditional S0 if n < 200."""
    f = fl[aln]
    lo = [b[0] for b in fl["bin_edges_len_a"]]
    med, q95, src = [], [], []
    for b in fl["bins"]:
        ok = f["S0_ok"][b]["LDDT_global"]
        if ok["n"] >= 200:
            med.append(ok["median"]); q95.append(ok["q95"]); src.append("S0_ok")
        else:
            u = f["S0"][b]["LDDT_global"]
            med.append(u["median"]); q95.append(u["q95"]); src.append("S0_fallback")
    return {"bin_lo": lo, "median": med, "q95": q95, "source": src}


class Tuning:
    """tuning.yaml values with provisional fall-backs (only where allowed); records what was provisional."""

    def __init__(self, root, provisional=False):
        import yaml
        self.root = Path(root)
        p = self.root / "config" / "tuning.yaml"
        self.raw = yaml.safe_load(p.read_text()) if p.exists() else {}
        self.md5 = md5_file(p)
        self.provisional = provisional
        self.used_provisional = []
        self.floors = _floors(self.root)
        self.primary = self.raw.get("primary_alignment") or "tdiaa"

    def _get(self, path, prov_value, why=""):
        o = self.raw
        for k in (path if isinstance(path, (list, tuple)) else path.split(".")):
            o = o.get(k) if isinstance(o, dict) else None
        if o is not None and not (isinstance(o, dict) and _nulls(o)):
            return o
        if not self.provisional:
            if prov_value is None:      # optional quantity (e.g. power when computed here)
                return None
            raise ValueError(f"tuning.yaml: '{path}' is null and --provisional is not set")
        self.used_provisional.append(f"{path}={why or prov_value}")
        return prov_value

    # -- scalar / table accessors
    def alpha(self):
        return {"aa": float(self._get("alpha.aa", 1.0)), "tdi": float(self._get("alpha.tdi", 1.0))}

    def eta(self, aln):
        fl = self.floors[aln] if self.floors else {}
        e_t = fl.get("eta_tdi_ml", 0.0)
        e_raw = fl.get("eta_lddt_raw", 0.0)
        s0 = self.s0(aln)
        # provisional eta_lddtD: first-order D ~ (1 - S)/(1 - S0) near S = 1, S0 at the 200-300 bin
        e_D = e_raw / (1 - s0["median"][min(2, len(s0["median"]) - 1)])
        return {"tdi_ml": float(self._get(f"eta.{aln}.tdi_ml", e_t, f"floors.json[{aln}].eta_tdi_ml")),
                "lddt_D": float(self._get(f"eta.{aln}.lddt_D", e_D, "eta_lddt_raw/(1-S0_med[200-300])")),
                "lddt_raw": float(self._get(f"eta.{aln}.lddt_raw", e_raw, f"floors.json[{aln}].eta_lddt_raw"))}

    def s0(self, aln):
        prov = _s0_from_floors(self.floors, aln) if self.floors else None
        return self._get(f"s0.{aln}", prov, "floors.json S0_ok (C28)")

    def kappa(self, aln):
        return self._get(f"kappa.{aln}", {"bin_edges": [0.0], "aa_ml": [1.0], "tdi_ml": [1.0]}, "kappa=1")

    def rho(self, aln):
        return self._get(f"rho.{aln}", {"bin_edges": [0.0], "aa_tdi": [0.3], "aa_lddt": [0.3]}, "rho=0.3")

    def lddt_se_model(self, aln):
        return self._get(f"lddt_se_model.{aln}", None, "proxy sd_S=sqrt(S(1-S)/n_core)")

    def tau_misfit(self):
        return float(self._get("tau_misfit", 3.0))

    def B_E4(self):
        return int(self._get("B_E4", PROVISIONAL_B_E4))

    def min_t_snr(self):
        return self._get("e5.min_t_snr", 8)

    def power(self):
        return self._get(("e5", "power_f1.5"), None, "not tuned -> computed by spike-in on this split")

    def tip_sd(self, aln):
        """Per-protein 3Di tip-noise SD (E4 cell rule D30): tuning.yaml tip_sd.{aln} if present, else the floor_ident
        tdi_ml spread of floors.json, (q95 - q05)/3.29/sqrt(2) (an upper bound; same formula as the V3 tuning)."""
        o = (self.raw.get("tip_sd") or {}).get(aln) if isinstance(self.raw.get("tip_sd"), dict) else None
        if o is not None:
            return float(o)
        e = ((self.floors or {}).get(aln, {}).get("eta", {}) or {}).get("tdi_ml", {})
        if e.get("q95") is None:
            return 0.0
        return float((e["q95"] - e["q05"]) / 3.29 / np.sqrt(2))

    def v3_false_concave(self):
        return self._get("v3.false_concave_rate", None, "not tuned")

    def snapshot(self):
        return {"tuning_md5": self.md5, "provisional": self.provisional, "used_provisional": self.used_provisional}


# ======================================================================================== tables
def _bin_lookup(edges, vals, x):
    edges = np.asarray(edges, float); vals = np.asarray(vals, float)
    b = np.clip(np.searchsorted(edges, np.nan_to_num(np.asarray(x, float), nan=0.0), side="right") - 1, 0, len(vals) - 1)
    return vals[b]


def _with_cols(df, cols: dict):
    import polars as pl
    return df.with_columns([pl.Series(k, np.asarray(v)) for k, v in cols.items()])


def lddt_columns(df, s0, se_model, flags):
    """Fill lddt_D / lddt_D_se / lddt_sat when absent (tables not yet recomputed at tuned S0)."""
    from strucevo.pairs.ml import lddt_D as _lddtD
    S = df["lddt_global"].to_numpy().astype(float)
    have_D = np.isfinite(df["lddt_D"].to_numpy().astype(float)).any() if "lddt_D" in df.columns else False
    have_se = np.isfinite(df["lddt_D_se"].to_numpy().astype(float)).any() if "lddt_D_se" in df.columns else False
    cols = {}
    if not have_D:
        r = _lddtD(S, df["len_a"].to_numpy(), s0)
        cols.update({"lddt_D": r["lddt_D"], "lddt_sat": r["lddt_sat"]})
        flags.append("lddt_D_computed_in_driver")
    if not have_se:
        lo = np.asarray(s0["bin_lo"], float)
        b = np.clip(np.searchsorted(lo, df["len_a"].to_numpy().astype(float), side="right") - 1, 0, len(lo) - 1)
        s0m = np.asarray(s0["median"], float)[b]
        if se_model is not None:
            lb = np.clip(np.searchsorted(np.asarray(se_model["len_bin_lo"], float), df["len_a"].to_numpy(), side="right") - 1, 0, None)
            sb = np.clip(np.searchsorted(np.asarray(se_model["s_bin_lo"], float), S, side="right") - 1, 0, None)
            sd = np.asarray(se_model["sd_S"], float)
            sdS = sd[np.minimum(lb, sd.shape[0] - 1), np.minimum(sb, sd.shape[1] - 1)]
        else:
            n = np.maximum(df["n_core"].to_numpy().astype(float), 1.0)
            sdS = np.sqrt(np.clip(S * (1 - S), 1e-4, None) / n)
            flags.append("lddt_D_se_provisional_proxy")
        with np.errstate(invalid="ignore", divide="ignore"):
            cols["lddt_D_se"] = sdS / np.abs(S - s0m)
    return _with_cols(df, cols) if cols else df


def apply_kappa(df, kap):
    """SE used everywhere = kappa(bin of aa_ml) x Fisher SE (aa_ml_se, tdi_ml_se)."""
    x = df["aa_ml"].to_numpy().astype(float)
    cols = {}
    for m, sc in PAIR_COLS_SE.items():
        k = _bin_lookup(kap["bin_edges"], kap[m], x)
        cols[sc] = df[sc].to_numpy().astype(float) * k
    return _with_cols(df, cols)


def scalar_rho(rho_tab, key, x):
    """Pair-count-weighted mean of rho(bin) over the pairs entering a fit (engines take a scalar rho)."""
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if not len(x):
        return float(np.mean(rho_tab[key]))
    return float(np.mean(_bin_lookup(rho_tab["bin_edges"], rho_tab[key], x)))


def table_alpha(root, aln):
    """alpha the shard tables were computed at (from the .done sentinels; default 1.0)."""
    a = set()
    for d in sorted((Path(root) / "data" / "interim" / "pairs" / aln).glob("shard_*.done")):
        try:
            al = json.loads(d.read_text()).get("alpha") or {"aa": 1.0, "tdi": 1.0}
            a.add((float(al["aa"]), float(al["tdi"])))
        except Exception:
            pass
    return sorted(a)


def variant_tables(root, name, aln, split, alpha, log=print):
    """Recompute aa_*/tdi_* (core + split halves) from counts_shard_*.npz at `alpha` for rows of `split`, into
    ROOT/data/interim/pairs_variants/<name>/<aln>/shard_XXXX__<split>.parquet (cached by alpha). aa_ml_seq /
    aa_pid_seq are NOT changed (C26: selection variable fixed at alpha = 1). Never writes into pairs/<aln>."""
    import polars as pl
    from strucevo.pairs.ml import ml_from_counts
    if split == "confirm" and not confirm_allowed():
        raise ConfirmSplitGuard("variant_tables(split='confirm') requires CONFIRM_RUN=1")
    src = Path(root) / "data" / "interim" / "pairs" / aln
    # if apply_tuning.py already produced the variant, use it
    out = Path(root) / "data" / "interim" / "pairs_variants" / name / aln
    out.mkdir(parents=True, exist_ok=True)
    meta_p = out / f"_meta__{split}.json"
    sig = [[f.name, f.stat().st_size, int(f.stat().st_mtime)] for f in sorted(src.glob("*shard_[0-9][0-9][0-9][0-9]*"))]
    if meta_p.exists():
        meta = json.loads(meta_p.read_text())
        if meta.get("alpha") == alpha and meta.get("done") and meta.get("source_sig") == sig:
            return out
        for f in out.glob(f"shard_*__{split}.parquet"):     # stale: source tables or alpha changed
            f.unlink()
    for pq in sorted(src.glob("shard_[0-9][0-9][0-9][0-9].parquet")):
        k = pq.stem[6:10]
        df = pl.scan_parquet(pq).with_row_index("_r").filter(pl.col("split") == split).collect()
        if df.height == 0:
            continue
        guard_split(df["split"].to_numpy())
        z = np.load(src / f"counts_shard_{k}.npz")
        rows = df["_r"].to_numpy()
        assert np.all(z["pair_id"][rows] == df["pair_id"].to_numpy()), f"pair_id mismatch shard {k}"
        halves = df["pair_set"].str.starts_with("e6_").to_numpy()
        n = df.height
        cols = {}
        for key in ("aa", "tdi"):
            r = ml_from_counts(z[key][rows], key, alpha[key])
            cols.update({f"{key}_p": r["p"], f"{key}_ml": r["ml"], f"{key}_ml_se": r["se"], f"{key}_sat": r["sat"]})
            for h in ("hA", "hB"):
                v = np.full(n, np.nan); s = np.full(n, np.nan)
                if halves.any():
                    rh = ml_from_counts(z[f"{key}_{h}"][rows][halves], key, alpha[key])
                    v[halves] = rh["ml"]; s[halves] = rh["se"]
                cols[f"{key}_ml_{h}"] = v; cols[f"{key}_ml_{h}_se"] = s
        df = _with_cols(df.drop("_r"), cols)
        tmp = out / f"shard_{k}__{split}.parquet.tmp"
        df.write_parquet(tmp); os.replace(tmp, out / f"shard_{k}__{split}.parquet")
    meta_p.write_text(json.dumps({"alpha": alpha, "split": split, "done": True, "created": _now(), "source_sig": sig,
                                  "note": "aa_ml_seq unchanged (C26)"}))
    log(f"  variant {name}/{aln}/{split}: recomputed at alpha {alpha}")
    return out


def load_variant(path, split):
    import polars as pl
    if split == "confirm" and not confirm_allowed():
        raise ConfirmSplitGuard("load_variant(split='confirm') requires CONFIRM_RUN=1")
    files = sorted(Path(path).glob(f"shard_*__{split}.parquet"))
    if not files:
        raise FileNotFoundError(f"no variant shards in {path} for {split}")
    df = pl.read_parquet([str(f) for f in files]).filter(pl.col("split") == split)
    guard_split(df["split"].to_numpy())
    return df


# ======================================================================================== panel tree
def load_panel(root, flags, age_override=None):
    """PanelTree from data/interim/timetree/panel_tree.nwk (+ node_ages.tsv median/lo95/hi95). Leaf labels
    'UPxxxx_<taxid>' -> taxid. Nodes without an age get (parent + oldest child)/2 (flag age_imputed);
    nodes without a CI keep sigma = 0 in age draws (flag age_ci_missing).
    age_override: {node: (median, lo95, hi95)} replacing node_ages.tsv values BEFORE imputation (E4 deep-age
    sensitivities, §4.4); followed by the top-down repair child = min(child, parent - 0.1) (flag age_override_repair)."""
    import dendropy
    import polars as pl
    from .trees import PanelTree
    tt = Path(root) / "data" / "interim" / "timetree"
    t = dendropy.Tree.get(path=str(tt / "panel_tree.nwk"), schema="newick", preserve_underscores=True,
                          suppress_internal_node_taxa=True, suppress_leaf_node_taxa=False)
    for lf in t.leaf_node_iter():
        lf.taxon.label = str(int(str(lf.taxon.label).rsplit("_", 1)[-1]))
    na = pl.read_csv(tt / "node_ages.tsv", separator="\t", infer_schema_length=0)
    f = lambda v: float(v) if v not in (None, "", "NA") else np.nan
    ages = {r["node"]: (f(r["median"]), f(r["lo95"]), f(r["hi95"])) for r in na.iter_rows(named=True)}
    if age_override:
        missing = sorted(set(age_override) - set(ages))
        if missing:
            raise KeyError(f"age_override nodes not in node_ages.tsv: {missing}")
        ages.update(age_override)
    panel = PanelTree.from_dendropy(t, ages)
    age = panel.age.copy()
    miss = np.flatnonzero(~np.isfinite(age))
    if len(miss):
        oldest = np.zeros(len(age))
        for i in range(len(age) - 1, -1, -1):        # postorder
            if not panel.is_leaf[i]:
                ch = [age[c] if np.isfinite(age[c]) else oldest[c] for c in panel.children[i]]
                oldest[i] = max(ch) if ch else 0.0
        for i in range(len(age)):                    # preorder
            if not np.isfinite(age[i]):
                p = panel.parent[i]
                pa = age[p] if p >= 0 else oldest[i] * 1.1
                age[i] = 0.5 * (pa + oldest[i])
        panel = panel.with_ages(age)
        flags.append(f"age_imputed:{len(miss)}")
    if age_override:
        age = panel.age.copy(); nfix = 0
        for i in range(len(age)):                    # preorder (parents first)
            p = panel.parent[i]
            if p >= 0 and not panel.is_leaf[i] and age[i] > age[p] - 0.1:
                age[i] = age[p] - 0.1; nfix += 1
        if nfix:
            panel = panel.with_ages(age, panel.lo, panel.hi)
        flags.append(f"age_override_repair:{nfix}")
    nci = int(np.sum(~panel.is_leaf & ~(np.isfinite(panel.lo) & np.isfinite(panel.hi))))
    if nci:
        flags.append(f"age_ci_missing:{nci}")
    return panel


def deep_age_overrides(root):
    """{'ages_strassert2021': {node: (med, lo, hi)}, 'ages_betts2018': {...}} from
    data/interim/timetree/deep_ages_panel.tsv (scripts/tuning/map_deep_ages.py), or None if absent."""
    import polars as pl
    p = Path(root) / "data" / "interim" / "timetree" / "deep_ages_panel.tsv"
    if not p.exists():
        return None
    d = pl.read_csv(p, separator="\t", infer_schema_length=0, null_values=["NA", ""])
    out = {}
    for name, pre in (("ages_strassert2021", "s21"), ("ages_betts2018", "b18")):
        ov = {}
        for r in d.iter_rows(named=True):
            if r[f"{pre}_median"] is not None:
                ov[r["node"]] = (float(r[f"{pre}_median"]), float(r[f"{pre}_lo"]), float(r[f"{pre}_hi"]))
        out[name] = ov
    return out


# ======================================================================================== annotations
def cath_classes(root, hog_genes: dict, flags):
    """C24: HOG class = the TED/CATH class covering >= 60% of CATH-labelled TED-domain residues of the HOG's
    selected proteins ('mixed' otherwise; 'none' if no labelled residue). Literal plan rule (unlabelled domain
    residues in the denominator) returned as the sensitivity map."""
    import polars as pl
    p = Path(root) / "data" / "processed" / "annotations" / "protein_annot.parquet"
    if not p.exists():
        flags.append("cath_unavailable")
        return None, None
    genes = sorted({g for v in hog_genes.values() for g in v})
    a = (pl.scan_parquet(p).select("acc", "ted_dom_residues", "ted_frac_class1", "ted_frac_class2",
                                   "ted_frac_class3", "ted_frac_class4", "ted_frac_unlabeled")
         .filter(pl.col("acc").is_in(genes)).collect())
    d = {r[0]: r[1:] for r in a.iter_rows()}
    names = {0: "alpha", 1: "beta", 2: "alphabeta", 3: "few_ss"}
    lab, lit = {}, {}
    for h, gs in hog_genes.items():
        num = np.zeros(4); den_lab = 0.0; den_all = 0.0
        for g in gs:
            r = d.get(g)
            if r is None or r[0] in (None, 0):
                continue
            n = float(r[0]); fr = np.array([v or 0.0 for v in r[1:5]]); un = float(r[5] or 0.0)
            num += n * fr; den_lab += n * (1 - un); den_all += n
        for out, den in ((lab, den_lab), (lit, den_all)):
            if den <= 0:
                out[h] = "none"; continue
            k = int(np.argmax(num)); fr = num[k] / den
            out[h] = names[k] if (fr >= 0.6 and k < 3) else "mixed"
    if not d:
        flags.append("cath_unavailable")
        return None, None
    return lab, lit


# ======================================================================================== the driver
class EndpointDriver:
    def __init__(self, root=None, split="tune", provisional=False, sizes=None, outdir=None, log=print,
                 device=None, label=None, run_alt_alignments=True):
        self.root = Path(root or PROJECT_ROOT).resolve()
        self.split = split
        self.tun = Tuning(self.root, provisional)
        self.sizes = dict(PLAN_SIZES) | (sizes or {})
        if self.sizes.get("B4") is None:
            self.sizes["B4"] = self.tun.B_E4()
        if self.sizes.get("B4_sens") is None:
            self.sizes["B4_sens"] = max(20, self.sizes["B4"] // 5)
        self.outdir = Path(outdir) if outdir else self.root / "results" / "endpoints"
        self.log = log
        self.device = device
        self.label = label
        self.flags = []
        self.timing = {}
        self.run_alt = run_alt_alignments
        self.primary = self.tun.primary

    # ------------------------------------------------------------------ data
    def _t(self, k, t0):
        self.timing[k] = round(time.time() - t0, 2)

    def load(self):
        import polars as pl
        t0 = time.time()
        r = self.root
        self.panel = load_panel(r, self.flags)
        hogs = pl.read_parquet(r / "data" / "interim" / "hogs.parquet")
        hogs = hogs.filter(pl.col("split") == self.split)
        guard_split(hogs["split"].to_numpy())
        self.hogs = hogs
        self.hog_species = {h: np.asarray(s, np.int64) for h, s in zip(hogs["hog"].to_list(), hogs["selected"].to_list())}
        self.hog_genes = {h: ([] if g is None else list(g)) for h, g in zip(hogs["hog"].to_list(), hogs["selected_genes"].to_list())}
        self.single_copy = set(hogs.filter(pl.col("n_genes_filtered") == pl.col("n_species_1to1"))["hog"].to_list())
        self.alpha = self.tun.alpha()
        self.tabs, self.e6tabs, self.pflags = {}, {}, {}
        for aln in ALNS:
            if aln != self.primary and not self.run_alt:
                continue
            fl = []
            ta = table_alpha(r, aln)
            df = load_pairs(aln, self.split, root=r)
            if ta and (ta != [(self.alpha["aa"], self.alpha["tdi"])]):
                fl.append(f"pairs_alpha_{ta}_ne_tuned_{self.alpha}->variant_alpha_tuned")
                vp = variant_tables(r, "alpha_tuned", aln, self.split, self.alpha, self.log)
                df = load_variant(vp, self.split)
            df = lddt_columns(df, self.tun.s0(aln), self.tun.lddt_se_model(aln), fl)
            df = apply_kappa(df, self.tun.kappa(aln))
            self.tabs[aln] = df.filter(pl.col("pair_set") == "primary")
            self.e6tabs[aln] = df.filter(pl.col("pair_set").str.starts_with("e6_"))
            self.pflags[aln] = sorted(set(fl))
            self.log(f"  loaded {aln}: {self.tabs[aln].height} primary, {self.e6tabs[aln].height} e6 pairs; {fl}")
        p = self.tabs[self.primary]
        sp = np.r_[p["species_a"].to_numpy(), p["species_b"].to_numpy()]
        pv = np.r_[p["mean_plddt_a"].to_numpy(), p["mean_plddt_b"].to_numpy()]
        sp_tsv = r / "results" / "panel" / "species.tsv"
        if sp_tsv.exists():          # §1.5 results/panel/species.tsv (taxid, mean_plddt, n_proteins)
            st = pl.read_csv(sp_tsv, separator="\t")
            self.species_plddt = dict(zip(st["taxid"].cast(pl.Int64).to_list(), st["mean_plddt"].cast(pl.Float64).to_list()))
        else:
            self.species_plddt = {int(s): float(np.nanmean(pv[sp == s])) for s in np.unique(sp)}
            self.flags.append("species_plddt_from_pair_tables")
        self.cath, self.cath_literal = cath_classes(r, self.hog_genes, self.flags)
        self._t("load", t0)

    # ------------------------------------------------------------------ branches
    def fit(self, tab, aln, hog_subset=None, panel=None):
        from .branch import fit_branches
        import polars as pl
        if hog_subset is not None:
            tab = tab.filter(pl.col("hog").is_in(list(hog_subset)))
        eta = self.tun.eta(aln)
        hs = {h: self.hog_species[h] for h in np.unique(tab["hog"].to_numpy()) if h in self.hog_species}
        return fit_branches(tab, self.panel if panel is None else panel, metrics=("aa_ml", "tdi_ml", "lddt_D"),
                            nu={"tdi_ml": eta["tdi_ml"] / 2, "lddt_D": eta["lddt_D"] / 2}, aln=aln,
                            hog_species=hs, tau_misfit=self.tau, se_method="model", misfit=True)

    def save_branches(self, bf, aln):
        out = self.root / "data" / "interim" / "branches"
        if self.split != "confirm":
            out = out / self.split
        out.mkdir(parents=True, exist_ok=True)
        bf.edge_frame(self.panel).write_parquet(out / f"{aln}.parquet")
        bf.fit_frame().write_parquet(out / f"{aln}_fit.parquet")
        import polars as pl
        pf = bf.pair_flags
        pl.DataFrame({"row": pf["row"].astype(np.int64), "hidden_para": pf["hidden_para"].astype(bool),
                      "resid_aa_z": pf["resid_aa_z"].astype(float)}).write_parquet(out / f"{aln}_pairflags.parquet")

    def branches(self):
        import polars as pl
        t0 = time.time()
        self.tau = self.tun.tau_misfit()
        self.bf = {}
        for aln, tab in self.tabs.items():
            self.bf[aln] = self.fit(tab, aln)
            self.save_branches(self.bf[aln], aln)
        # hidden paralogy (residual rule §2.5.7) OR out-paralog label (C27), primary alignment rows
        tab = self.tabs[self.primary]
        pf = self.bf[self.primary].pair_flags
        hp = np.zeros(tab.height, bool)
        if len(pf["row"]):
            hp[pf["row"].astype(int)[pf["hidden_para"].astype(bool)]] = True
        hp |= tab["relation"].to_numpy().astype(str) != "ortholog"
        self.hp_rows = np.flatnonzero(hp)
        self.bf_hp = self.fit(tab.filter(pl.Series(~hp)), self.primary)
        self.bf_sc = self.fit(tab, self.primary, hog_subset=self.single_copy) if self.single_copy else None
        self.hog_cols()
        self._t("branches", t0)

    def hog_cols(self):
        """hogs-level columns (§1.1) -> data/interim/hogs_endpoint_cols_{split}.parquet."""
        import polars as pl
        mk = {"aa_ml": "aa", "tdi_ml": "tdi", "lddt_D": "lddt"}
        cols = {"hog": self.hogs["hog"].to_list()}
        idx = {h: i for i, h in enumerate(cols["hog"])}
        H = len(idx)
        self.excl = {}
        for aln, bf in self.bf.items():
            F = bf.fits
            for m, s in mk.items():
                v = np.full(H, np.nan)
                sel = F["metric"] == m
                for h, ff in zip(F["hog"][sel], F["fail_frac"][sel]):
                    v[idx[h]] = ff
                cols[f"fail_frac_{s}_{aln}"] = v
                ex = v > 0.10
                cols[f"excl_{s}_{aln}"] = ex
                self.excl[(m, aln)] = {h for h, e in zip(cols["hog"], ex) if e}
        F = self.bf[self.primary].fits
        sel = F["metric"] == "aa_ml"
        mr = np.full(H, np.nan)
        for h, v in zip(F["hog"][sel], F["misfit_ratio"][sel]):
            mr[idx[h]] = v
        cols["misfit_ratio"] = mr
        cols["excl_misfit"] = np.where(np.isfinite(mr), mr > self.tau, False)
        cols["cath_class"] = [self.cath.get(h, "none") if self.cath else "cath_unavailable" for h in cols["hog"]]
        cols["cath_class_literal"] = [self.cath_literal.get(h, "none") if self.cath_literal else "cath_unavailable"
                                      for h in cols["hog"]]
        df = pl.DataFrame(cols).with_columns(pl.lit(self.split).alias("split"))
        p = self.root / "data" / "interim" / f"hogs_endpoint_cols_{self.split}.parquet"
        df.write_parquet(p)
        self.hogcols = df

    # ------------------------------------------------------------------ helpers
    def _rho(self, aln, key, tab):
        return scalar_rho(self.tun.rho(aln), key, tab["aa_ml"].to_numpy())

    def _rho_pair(self, aln, key, tab):
        """Per-pair rho(bin of the pair's aa_ml) (per-pair Sigma_ij in SIMEX; D29)."""
        rt = self.tun.rho(aln)
        x = tab["aa_ml"].to_numpy().astype(float)
        if np.isscalar(rt.get(key)) or len(rt.get("bin_edges", [])) <= 1:
            return np.full(len(x), scalar_rho(rt, key, x))
        return _bin_lookup(rt["bin_edges"], rt[key], x)

    def _impute_c14(self, tab, ycol, ysecol, ysat):
        """C14: failed / saturated (y or AA) pairs -> y = the HOG's maximum eligible y, SE = median eligible
        SE; x = aa_ml if finite and > 0, else aa_ml_seq."""
        hog = tab["hog"].to_numpy().astype(str)
        x = tab["aa_ml"].to_numpy().astype(float); xs = tab["aa_ml_seq"].to_numpy().astype(float)
        y = tab[ycol].to_numpy().astype(float); sy = tab[ysecol].to_numpy().astype(float)
        sx = tab["aa_ml_se"].to_numpy().astype(float)
        st = tab["status"].to_numpy().astype(str)
        sat = tab[ysat].fill_null(False).to_numpy().astype(bool) | tab["aa_sat"].fill_null(False).to_numpy().astype(bool)
        elig = (st == "ok") & ~sat & np.isfinite(x) & np.isfinite(y) & np.isfinite(sy) & (x > 0)
        bad = ~elig
        hmax = {}
        for h, v in zip(hog[elig], y[elig]):
            hmax[h] = max(hmax.get(h, -np.inf), v)
        med_sy = np.nanmedian(sy[elig]); med_sx = np.nanmedian(sx[elig])
        y2 = np.where(bad, [hmax.get(h, np.nan) for h in hog], y)
        x2 = np.where(bad & ~(np.isfinite(x) & (x > 0)), xs, x)
        sx2 = np.where(bad & ~np.isfinite(sx), med_sx, sx)
        return _with_cols(tab, {ycol: y2, ysecol: np.where(bad, med_sy, sy), "aa_ml": x2, "aa_ml_se": sx2,
                                "status": np.where(bad & np.isfinite(y2) & np.isfinite(x2) & (x2 > 0), "ok", st),
                                ysat: np.zeros(len(y), bool), "aa_sat": np.zeros(len(y), bool)})

    def _alpha_variants(self, aln):
        out = {}
        for name, f in (("alpha_x1.5", 1.5), ("alpha_div1.5", 1 / 1.5)):
            al = {"aa": self.alpha["aa"] * f, "tdi": self.alpha["tdi"] * f}
            p = variant_tables(self.root, name, aln, self.split, al, self.log)
            df = load_variant(p, self.split)
            fl = []
            df = lddt_columns(df, self.tun.s0(aln), self.tun.lddt_se_model(aln), fl)
            df = apply_kappa(df, self.tun.kappa(aln))
            import polars as pl
            out[name] = df.filter(pl.col("pair_set") == "primary")
        return out

    def _all_columns(self, aln):
        """E1/E2 sensitivity 'All aligned columns' (§4.1): aa_ml / tdi_ml (+ SE, sat, p) recomputed at the tuned alpha
        from AA/3Di counts over EVERY aligned column of the same Foldseek alignment (strucevo.pairs.allcols; the
        alignment is re-run, core counts verified identical). kappa(bin) applied as for the core columns; status,
        aa_ml_seq (selection, C26), lddt_* unchanged. Branch engine refitted on the all-column table (PB design).
        Cached per alignment. Returns (primary-pair table, BranchFit, meta)."""
        cache = getattr(self, "_allcols_cache", {})
        if aln in cache:
            return cache[aln]
        import polars as pl
        from strucevo.pairs.allcols import all_column_counts
        from strucevo.pairs.ml import ml_from_counts
        f, meta = all_column_counts(self.root, aln, self.split, threads=int(os.environ.get("ALLCOLS_THREADS", "40")),
                                    log=self.log)
        z = np.load(f)
        tab = self.tabs[aln]
        pos = {p: i for i, p in enumerate(z["pair_id"].tolist())}
        rows = np.array([pos.get(p, -1) for p in tab["pair_id"].to_list()])
        good = rows >= 0
        good[good] = z["core_match"][rows[good]] & z["aligned"][rows[good]]
        cols = {}
        for key in ("aa", "tdi"):
            cnt = np.zeros((tab.height, 20, 20), np.int64)
            cnt[good] = z[f"{key}_all"][rows[good]]
            r = ml_from_counts(cnt, key, self.alpha[key], self.device)
            nan = lambda v: np.where(good, v, np.nan)
            cols.update({f"{key}_p": nan(r["p"]), f"{key}_ml": nan(r["ml"]), f"{key}_ml_se": nan(r["se"]),
                         f"{key}_sat": np.where(good, r["sat"], False)})
        cols["status"] = np.where(good, tab["status"].to_numpy().astype(str), "fail_allcols_recompute")
        df = apply_kappa(_with_cols(tab, cols), self.tun.kappa(aln))
        bfa = self.fit(df, aln)
        meta = {k: v for k, v in meta.items() if k != "source_sig"} | {"n_rows_not_matched": int((~good).sum())}
        cache[aln] = (df, bfa, meta)
        self._allcols_cache = cache
        return cache[aln]

    def _finish(self, rec, name):
        bad = []
        for k, v in (rec.get("sensitivity") or {}).items():      # graceful degradation: mark empty sensitivities
            est = v.get("estimate") if isinstance(v, dict) else None
            fin = (isinstance(est, dict) and any(isinstance(x, (int, float)) and np.isfinite(x) for x in est.values())) or \
                  (isinstance(est, (int, float, np.floating)) and np.isfinite(est))
            if not fin and isinstance(v, dict) and v.get("decision") not in ("unavailable", "not_estimable"):
                v["decision"] = "not_estimable"; bad.append(k)
        if bad:
            rec["flags"] = list(rec.get("flags") or []) + ["sensitivity_not_estimable:" + ",".join(bad)]
        rec["plan_md5"] = md5_file(self.root / "ANALYSIS_PLAN.md") or rec.get("plan_md5")
        rec["tuning_md5"] = self.tun.md5
        rec["run_meta"] = {"root": str(self.root), "split": self.split, "primary_alignment": self.primary,
                           "label": self.label, "created": _now(), "sizes": self.sizes,
                           "tuning": self.tun.snapshot(), "pair_table_flags": self.pflags,
                           "split_md5": md5_file(self.root / "data" / "interim" / "split.tsv"),
                           "driver_flags": self.flags}
        if self.tun.provisional and self.tun.used_provisional:
            rec["flags"] = list(rec.get("flags") or []) + ["provisional_tuning"]
        if self.label:
            rec["flags"] = list(rec.get("flags") or []) + [self.label]
        write_json(rec, self.outdir / f"{name}.json")
        return rec

    def _failed(self, name, e):
        rec = endpoint_record(endpoint=name, estimand=None, decision="not_estimable",
                              flags=[f"error:{type(e).__name__}:{str(e)[:300]}"])
        rec["traceback"] = traceback.format_exc()[-3000:]
        self.log(f"  {name} FAILED: {e}")
        return self._finish(rec, name)

    # ------------------------------------------------------------------ endpoints
    def run_curve(self, ep):
        from .endpoints_curve import run_curve_endpoint, CURVE_COLS
        S = self.sizes
        t0 = time.time()
        ycol, ysecol, ysat = CURVE_COLS[ep]
        mkey, rkey, ekey = ("tdi_ml", "aa_tdi", "tdi_ml") if ep == "E1" else ("lddt_D", "aa_lddt", "lddt_D")
        aln = self.primary
        tab = self.tabs[aln]
        eta = self.tun.eta(aln)
        rho = self._rho(aln, rkey, tab)
        common = dict(kx=1.0, ky=1.0, B_s=S["Bs"], B_s_boot=S["Bs_boot"], device=self.device,
                      v3_false_concave=self.tun.v3_false_concave())
        # alignment sensitivities (full endpoint at B_alt)
        alt = {}
        for a2 in ALNS:
            if a2 == aln or a2 not in self.tabs:
                continue
            try:
                r2, cv2, _ = run_curve_endpoint(self.tabs[a2], ep, eta=self.tun.eta(a2)[ekey], rho=self._rho(a2, rkey, self.tabs[a2]),
                                                excluded_hogs=self.excl[(mkey, a2)], B=S["B_alt"], B_naive=S["B_alt"],
                                                rho_pair=self._rho_pair(a2, rkey, self.tabs[a2]),
                                                bf=self.bf.get(a2), **common)
                alt[a2] = (r2["estimate"], r2["se_boot"])
                self.alt_curves = getattr(self, "alt_curves", {}); self.alt_curves[(ep, a2)] = cv2
            except Exception as e:
                self.log(f"  {ep} alignment {a2} failed: {e}")
                alt[a2] = (np.nan, np.nan)
        # sensitivities
        hp_mask = np.ones(tab.height, bool); hp_mask[self.hp_rows] = False
        hogarr = tab["hog"].to_numpy().astype(str)
        sens = {"eta_free": {"eta_free": True},
                "d_aa_le_1": {"x_range": (0.05, 1.0)},
                "failures_imputed_C14": {"tab": self._impute_c14(tab, ycol, ysecol, ysat), "bf": self.bf[aln]},
                "hidden_paralogy_removed": {"extra_mask": hp_mask},
                "single_copy_hogs": {"extra_mask": np.isin(hogarr, np.asarray(sorted(self.single_copy), dtype=str))}}
        try:
            for k, v in self._alpha_variants(aln).items():
                sens[k] = {"tab": v}
        except Exception as e:
            self.flags.append(f"alpha_variant_failed:{e}")
        allcols_meta, allcols_err = None, None
        try:
            t_all, bf_all, allcols_meta = self._all_columns(aln)
            sens["all_aligned_columns"] = {"tab": t_all, "bf": bf_all,
                                           "note": "aa_ml/tdi_ml over all aligned columns (core = pLDDT>=70 both)"
                                                   + ("; y = lddt_D unchanged (lddt_global is already all-column)" if ep == "E2" else "")}
        except Exception as e:
            allcols_err = f"{type(e).__name__}: {str(e)[:300]}"
            self.log(f"  {ep} all-columns sensitivity failed: {allcols_err}")
        rec, curve, _ = run_curve_endpoint(tab, ep, eta=eta[ekey], rho=rho, excluded_hogs=self.excl[(mkey, aln)],
                                           rho_pair=self._rho_pair(aln, rkey, tab), bf=self.bf[aln],
                                           B=S["B1"], B_naive=S["B_naive"], sensitivities=sens, B_sens=S["B_sens"],
                                           eta_raw=eta["lddt_raw"] if ep == "E2" else None,
                                           alignment_results={a: v for a, v in alt.items() if np.all(np.isfinite(v))}, **common)
        for a2, v in alt.items():
            if not np.all(np.isfinite(v)):
                rec["sensitivity"][f"aln_{a2}"] = {"estimate": None, "ci98": None, "decision": "not_estimable"}
                rec["alignment_robust"] = None
                rec["flags"].append(f"alignment_{a2}_not_estimable")
        if allcols_err is not None:
            rec["sensitivity"]["all_aligned_columns"] = {"estimate": None, "ci98": None, "decision": "not_estimable",
                                                         "error": allcols_err}
            rec["flags"].append("all_columns_sensitivity_not_estimable")
        else:
            rec["sensitivity"]["all_aligned_columns"]["recompute_meta"] = allcols_meta
        rec["rho_scalar"] = rho
        rec["sensitivity_failures_imputed_rule"] = "C14: HOG max eligible y, SE = median eligible SE"
        self._t(ep, t0)
        self.curves = getattr(self, "curves", {}); self.curves[ep] = curve
        return self._finish(rec, ep)

    def run_E3(self):
        from .endpoints_branch import run_E3
        t0 = time.time()
        S = self.sizes
        alt = {}
        for a2, bf in self.bf.items():
            if a2 != self.primary:
                try:
                    alt[a2] = run_E3(bf, self.panel, self.tau, B=S["B3_alt"])["estimate"]
                except Exception as e:
                    self.log(f"  E3 aln {a2} failed: {e}")
        hogs = self.hogs
        # groups: length tertile (mean len over primary pairs), pLDDT tertile, stratum
        import polars as pl
        t = self.tabs[self.primary].group_by("hog").agg(((pl.col("len_a") + pl.col("len_b")) / 2).mean().alias("L"))
        Lm = dict(zip(t["hog"], t["L"]))
        def tert(d):
            v = np.array(list(d.values()), float)
            q = np.nanquantile(v, [1 / 3, 2 / 3]) if len(v) else [0, 0]
            return {h: f"T{1 + int(x > q[0]) + int(x > q[1])}" for h, x in d.items()}
        groups = {"length_tertile": tert(Lm), "plddt_tertile": tert(dict(zip(hogs["hog"], hogs["mean_plddt"]))),
                  "stratum": {h: str(s) for h, s in zip(hogs["hog"], hogs["stratum"])}}
        rec = run_E3(self.bf[self.primary], self.panel, self.tau, cath_class=self.cath, hog_groups=groups,
                     B=S["B3"], alignment_k=alt or None)
        if self.cath_literal:
            r2 = run_E3(self.bf[self.primary], self.panel, self.tau, cath_class=self.cath_literal, B=S["B3_alt"])
            rec["sensitivity"]["cath_literal_rule"] = {"estimate": (r2.get("cath_contrast") or {}).get("delta"),
                                                       "ci98": None, "decision": None, "detail": r2.get("cath_contrast")}
        if self.cath is None:
            rec["flags"].append("cath_unavailable")
        elif rec.get("cath_contrast") is None:
            rec["flags"].append("cath_contrast_not_estimable")
        rec["cath_rule"] = "C24 ted_class_labeled-equivalent (>= 60% of CATH-labelled TED residues)"
        self._t("E3", t0)
        return self._finish(rec, "E3")

    def run_E4(self):
        from .endpoints_branch import run_E4, edge_view
        from .endpoints_misc import e4_robust_alpha_timescale
        import polars as pl
        t0 = time.time()
        S = self.sizes
        bf = self.bf[self.primary]
        extra = {}
        tv = self.tun.tip_sd(self.primary) ** 2      # D30 cell rule: per-protein tip-noise variance
        def point(bfv, B):
            r = run_E4(bfv, self.panel, self.tau, B=B, sensitivities=False, tip_var=tv)
            return {"estimate": r["estimate"], "ci98": r["ci98"], "p_unadjusted": r["p_value"], "B": B,
                    "decision": None, "flags": r["flags"]}
        try:
            for name, tab in self._alpha_variants(self.primary).items():
                extra[name] = point(self.fit(tab, self.primary), S["B4_sens"])
        except Exception as e:
            self.flags.append(f"E4_alpha_variant_failed:{e}")
        tab = self.tabs[self.primary]
        mx = tab.group_by("hog").agg(pl.col("aa_ml_seq").max().alias("m"))
        sub = set(mx.filter(pl.col("m") <= 1.0)["hog"].to_list())
        try:
            extra["aa_ml_seq_le_1_hogs"] = point(self.fit(tab, self.primary, hog_subset=sub), S["B4_sens"]) | {
                "note": "HOGs whose primary pairs all have aa_ml_seq <= 1", "n_hogs": len(sub)}
        except Exception as e:
            extra["aa_ml_seq_le_1_hogs"] = {"estimate": None, "ci98": None, "decision": "not_estimable", "error": str(e)[:200]}
        # deep-age sensitivities (§4.4): literature ages on the deep panel nodes (TimeTree elsewhere), branch engine
        # refitted on the re-aged panel (D1 root-merged shares use ages), then E4 at B4_sens (as the alpha variants)
        deep = deep_age_overrides(self.root)
        for name in ("ages_strassert2021", "ages_betts2018"):
            ov = (deep or {}).get(name)
            if not ov:
                extra[name] = {"estimate": None, "ci98": None, "decision": "unavailable",
                               "note": "no deep-age table data/interim/timetree/deep_ages_panel.tsv (scripts/tuning/map_deep_ages.py)"}
                continue
            try:
                pfl = []
                pan = load_panel(self.root, pfl, age_override=ov)
                bfd = self.fit(tab, self.primary, panel=pan)
                r = run_E4(bfd, pan, self.tau, B=S["B4_sens"], sensitivities=False, tip_var=tv)
                extra[name] = {"estimate": r["estimate"], "ci98": r["ci98"], "p_unadjusted": r["p_value"],
                               "B": S["B4_sens"], "decision": None, "flags": r["flags"], "n_cells": r["n_units"],
                               "nodes_overridden": {k: list(v) for k, v in ov.items()}, "panel_flags": pfl,
                               "note": "literature ages on deep panel nodes (TimeTree median >= 500 My) only; "
                                       "TimeTree elsewhere; branch engine refitted on the re-aged panel"}
            except Exception as e:
                extra[name] = {"estimate": None, "ci98": None, "decision": "not_estimable", "error": str(e)[:300]}
        signs = {}
        for a2, bfa in self.bf.items():
            if a2 != self.primary:
                try:
                    extra[f"aln_{a2}"] = point(bfa, S["B4_alt"]); signs[a2] = np.sign(extra[f"aln_{a2}"]["estimate"])
                except Exception as e:
                    extra[f"aln_{a2}"] = {"estimate": None, "ci98": None, "decision": "not_estimable", "error": str(e)[:200]}
        rec = run_E4(bf, self.panel, self.tau, B=S["B4"], sensitivities=True, B_sens=S["B4_sens"], extra_sens=extra,
                     tip_var=tv)
        need = ("alpha_x1.5", "alpha_div1.5", "ages_strassert2021", "ages_betts2018")
        sv = rec["sensitivity"]
        avail = [n for n in need if isinstance(sv.get(n), dict) and sv[n].get("estimate") is not None
                 and np.isfinite(sv[n]["estimate"])]
        rec["robust_alpha_timescale_detail"] = {
            n: ({"sign_same": bool(np.sign(sv[n]["estimate"]) == np.sign(rec["estimate"])),
                 "p_unadjusted": sv[n].get("p_unadjusted")} if n in avail else "unavailable") for n in need}
        if len(avail) < len(need):
            rec["flags"].append("robust_alpha_timescale_incomplete:" + ",".join(n for n in need if n not in avail))
            rec["robust_alpha_timescale"] = False
            rec["robust_alpha_timescale_note"] = "a required alpha/deep-age sensitivity is unavailable -> cannot be established"
        else:
            rec["robust_alpha_timescale"] = bool(e4_robust_alpha_timescale(rec, need[:2], need[2:]))
            rec["robust_alpha_timescale_note"] = ("§4.4: same sign as beta-hat and unadjusted p < 0.05 in alpha x1.5, "
                                                  "alpha /1.5, Strassert 2021 and Betts 2018 deep ages")
        if signs:
            rec["alignment_robust"] = bool(all(s == np.sign(rec["estimate"]) for s in signs.values()))
            rec["alignment_robust_rule"] = "same sign of beta across tdiaa/seqonly/tmalign"
        if self.split == "tune":       # D13 negative control: tau permuted across HOG edges
            nE = len(edge_view(bf, self.panel, tau_misfit=self.tau).child)
            cover = []
            for k in range(S["neg_perms"]):
                rng = stage_rng("age_perm", k)
                r = run_E4(bf, self.panel, self.tau, B=S["neg_B"], sensitivities=False, tau_perm=rng.permutation(nE),
                           tip_var=tv)
                cover.append(bool(r["ci98"][0] <= 0 <= r["ci98"][1]))
            rec["negative_control_D13"] = {"n_perm": len(cover), "n_cover_0": int(sum(cover)),
                                           "pass": bool(sum(cover) >= 0.9 * len(cover)), "B_per_perm": S["neg_B"]}
        self._t("E4", t0)
        self.rec4 = rec
        return self._finish(rec, "E4")

    def run_E6(self):
        from .endpoints_misc import run_E6
        import polars as pl
        t0 = time.time()
        aln = self.primary
        trip, pairs = self.e6_inputs(aln)
        if trip.height == 0:
            raise RuntimeError("no E6 triplets for this split")
        rec = run_E6(trip, pairs, self.curves["E1"], B=self.sizes["B6"])
        # alignment sensitivities (needs the E1 curve of the same alignment: use the primary curve, flagged)
        for a2 in ALNS:
            if a2 != aln and a2 in self.e6tabs:
                try:
                    t2, p2 = self.e6_inputs(a2)
                    cv = getattr(self, "alt_curves", {}).get(("E1", a2))
                    r2 = run_E6(t2, p2, cv or self.curves["E1"], B=self.sizes["B6"] // 5)
                    rec["sensitivity"][f"aln_{a2}"] = {"estimate": r2["estimate"], "ci98": r2["ci98"], "E6b": r2["E6b"],
                                                       "p_value": r2["p_value"], "decision": None,
                                                       "note": f"reference curve = E1 curve of {a2 if cv else 'the primary alignment (alt fit missing)'}"}
                except Exception as e:
                    rec["sensitivity"][f"aln_{a2}"] = {"estimate": None, "ci98": None, "decision": "not_estimable", "error": str(e)[:200]}
        if "E2" in self.curves:      # secondary: lddt_D on full columns (not split-half valid)
            try:
                p2 = pairs.with_columns(pl.col("aa_ml").alias("aa_ml_hA"), pl.col("aa_ml").alias("aa_ml_hB"),
                                        pl.col("lddt_D").alias("lddt_D_hA"), pl.col("lddt_D").alias("lddt_D_hB"))
                r2 = run_E6(trip, p2, self.curves["E2"], B=self.sizes["B6"] // 5, metric=("aa_ml", "lddt_D"), secondary=True)
                rec["secondary_lddtD"] = {k: r2[k] for k in ("E6a", "E6b", "ci98", "ci98_b", "p_value", "n_units", "flags")}
            except Exception as e:
                rec["secondary_lddtD"] = {"error": str(e)[:200]}
        rec["replication_direction"]["expected"] = "positive (C19-C21: MDO shows MORE structural divergence)"
        self.e6_rec = rec
        self._t("E6", t0)
        return self._finish(rec, "E6")

    def e6_inputs(self, aln):
        """pairlists/e6_triplets.tsv (P1,P2,O1,O2 comma lists) + e6_* pair rows -> run_E6 schema."""
        import polars as pl
        p = self.root / "data" / "interim" / "pairlists" / "e6_triplets.tsv"
        if self.split == "confirm" and not confirm_allowed():
            raise ConfirmSplitGuard("e6 triplets of split confirm require CONFIRM_RUN=1")
        t = pl.read_csv(p, separator="\t", infer_schema_length=0).filter(pl.col("split") == self.split)
        guard_split(t["split"].to_numpy())
        pairs = self.e6tabs[aln]
        pp = pairs.filter(pl.col("pair_set") == "e6_para").select("gene_a", "gene_b", "aa_pid_seq", "aa_ml_seq")
        key = {}
        for a, b, pid, d in pp.iter_rows():
            key[(a, b)] = key[(b, a)] = (pid, d)
        sl = lambda s: [g for g in (s or "").split(",") if g]
        rows = []
        for r in t.iter_rows(named=True):
            pid, d = key.get((r["P1"], r["P2"]), (np.nan, np.nan))
            rows.append({"hog": r["hog"], "dcg_id": r["dcg_id"], "species": r["species"], "copy1": r["P1"],
                         "copy2": r["P2"], "event_level": r["event_level"], "aa_pid_seq": pid, "aa_ml_seq": d,
                         "o1": sl(r["O1"]), "o2": sl(r["O2"])})
        trip = pl.DataFrame(rows, schema_overrides={"o1": pl.List(pl.String), "o2": pl.List(pl.String)}) if rows else pl.DataFrame()
        return trip, pairs

    def run_E5(self, compute_power):
        from .endpoints_lineage import run_E5
        t0 = time.time()
        S = self.sizes
        aln = self.primary
        tab = self.tabs[aln]
        rho = self._rho(aln, "aa_tdi", tab)
        tuned_power = self.tun.power()
        from .endpoints_lineage import E5Data
        from .endpoints_branch import edge_view
        guards, dropped = {}, {}
        fe = {k: self.rec4[k] for k in ("mu", "gamma", "x_center", "z_center")} | {"beta": self.rec4["estimate"]}
        for gname, bfg in (("hidden_paralogy_removed", self.bf_hp), ("single_copy_hogs", self.bf_sc)):
            if bfg is None:
                dropped[gname] = "no HOGs"; continue
            try:
                evg = edge_view(bfg, self.panel, tau_misfit=self.tau)
                nl = len(E5Data(bfg, self.panel, evg, fe, self.species_plddt, min_t_snr=self.tun.min_t_snr()).lineages)
            except Exception as e:
                nl = 0
            if nl >= 3:
                guards[gname] = bfg
            else:
                dropped[gname] = f"n_lineages={nl} < 3"
        kw = dict(spike_reps=S["spike_reps"], spike_B=S["spike_B"]) if (compute_power or tuned_power is None) else \
            dict(spike_reps=0, spike_fs=(), power_override=float(tuned_power))
        rec, dat, core = run_E5(self.bf[aln], self.panel, self.rec4, tau_misfit=self.tau, species_plddt=self.species_plddt,
                                rho=rho, nu_t=self.tun.eta(aln)["tdi_ml"] / 2, P=S["P"], B=S["B5"], B_cv=S["B_cv"],
                                guards=guards, e6_significant=None, min_t_snr=self.tun.min_t_snr(), **kw)
        if "power_override" in kw:
            rec["power"] = {"f=1.5": {"detect_rate": float(tuned_power), "accept_rate": None, "source": "tuning.yaml"}}
        rec["rho_scalar"] = rho
        if dropped:
            rec["guards_not_estimable"] = dropped
            rec["flags"].append("guards_not_estimable:" + ",".join(dropped))
            for gname, why in dropped.items():
                rec["sensitivity"][gname] = {"estimate": None, "ci98": None, "decision": "not_estimable", "note": why}
            if rec["decision"] == "N0_rejected":   # rejection must hold under EVERY guard (§4.5)
                rec["decision"] = "artefact_suspect"
                rec["flags"].append("rejection_unverified_guard->artefact_suspect")
        rec["e6_gate"] = "pending_holm (finalised by finalize_E5 after Holm)"
        self._t("E5", t0)
        return self._finish(rec, "E5")

    # ------------------------------------------------------------------ run all
    def run_all(self, compute_power=None):
        t0 = time.time()
        self.outdir.mkdir(parents=True, exist_ok=True)
        compute_power = (self.split == "tune") if compute_power is None else compute_power
        self.load()
        self.branches()
        self.curves = {}
        done = {}
        for name, fn in (("E1", lambda: self.run_curve("E1")), ("E2", lambda: self.run_curve("E2")),
                         ("E3", self.run_E3), ("E4", self.run_E4), ("E6", self.run_E6),
                         ("E5", lambda: self.run_E5(compute_power))):
            self.log(f"[{time.time() - t0:7.1f}s] {name} ...")
            if name == "E5" and not hasattr(self, "rec4"):
                done[name] = self._failed(name, RuntimeError("E4 fixed effects unavailable"))
                continue
            try:
                done[name] = fn()
            except Exception as e:
                done[name] = self._failed(name, e)
            self.log(f"[{time.time() - t0:7.1f}s] {name}: estimate={done[name].get('estimate')} decision={done[name].get('decision')}")
        self.timing["total"] = round(time.time() - t0, 1)
        write_json({"timing_s": self.timing, "sizes": self.sizes, "flags": self.flags, "root": str(self.root),
                    "split": self.split, "tuning": self.tun.snapshot(), "label": self.label,
                    "n_pairs": {a: t.height for a, t in self.tabs.items()}, "n_hogs": self.hogs.height},
                   self.outdir / "run_meta.json")
        return done


def finalize_E5(outdir):
    """After Holm: apply the E6-significance gate of the E5 acceptance rule (p_E5 is unaffected: acceptance and
    'inconclusive' both use the IUT p)."""
    d = Path(outdir)
    p5, p6 = d / "E5.json", d / "E6.json"
    if not (p5.exists() and p6.exists()):
        return None
    r5 = json.loads(p5.read_text()); r6 = json.loads(p6.read_text())
    det = r5.get("decision_detail") or {}
    e6_sig = r6.get("decision") == "detected"
    det["e6_significant"] = e6_sig
    pw = det.get("power")
    if r5.get("decision") == "inconclusive" and det.get("margins_met") and e6_sig and pw is not None and pw >= 0.80:
        r5["decision"] = "N0_accepted"
    r5["decision_detail"] = det
    r5["e6_gate"] = f"applied after Holm: E6 {'detected' if e6_sig else 'not detected'}"
    write_json(r5, p5)
    return r5["decision"]
