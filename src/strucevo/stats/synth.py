"""Synthetic worlds for verification (ANALYSIS_PLAN §7 V2-V4, V6).

* `random_panel`      : dated ultrametric species tree (~100 species) with collapsed short internodes.
* `make_world`        : HOGs (clade-closed S_h of 8-30 species), per-HOG AA rates, shared lineage rate
                         multipliers, b_3Di = r_h * b_AA (N0), optional age effect and lineage spike-in.
* `observe_gauss`     : pair table with realistic SEs (tabulated from the ML engine), correlated AA/3Di
                         errors (kappa, rho) and per-protein tip noise.
* `observe_counts`    : pair table from `strucevo.distances.simulate_counts` + `ml_distance` (LG, Q.3Di.AF).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

import numpy as np

from .trees import PanelTree, prune

ROOT = Path(__file__).resolve().parents[3]
Q3DI_PATH = ROOT / "data" / "raw" / "models" / "Q.3Di.AF"


# ------------------------------------------------------------------------------------------ tree
def random_panel(n_species=100, root_age=1500.0, rng=None, min_internode=5.0, ci_sigma=0.1) -> PanelTree:
    """Random dated tree by recursive random splitting; internal ages = parent * U(0.35, 0.95) with tips at 0.
    Internodes shorter than `min_internode` My are collapsed into polytomies (ANALYSIS_PLAN §2.2).
    ci_sigma: log-scale SD of the age uncertainty (lo/hi = median * exp(-/+ 1.96 sigma))."""
    rng = np.random.default_rng(rng)
    parent, age, leaf = [-1], [root_age], [False]
    taxa = [None]
    stack = [(0, list(range(n_species)))]
    while stack:
        node, sp = stack.pop()
        k = rng.integers(1, len(sp))
        perm = rng.permutation(sp)
        for part in (perm[:k], perm[k:]):
            idx = len(parent)
            parent.append(node)
            if len(part) == 1:
                age.append(0.0); leaf.append(True); taxa.append(1000 + int(part[0]))
            else:
                # deeper clades keep more of the parent's age; small clades are young
                frac = rng.uniform(0.35, 0.95) * min(1.0, 0.25 + 0.75 * np.log(len(part)) / np.log(n_species))
                age.append(min(max(age[node] * frac, 2.0), age[node])); leaf.append(False); taxa.append(None)
                stack.append((idx, list(part)))
    parent = np.array(parent); age = np.array(age); leaf = np.array(leaf)
    # collapse short internal internodes
    changed = True
    alive = np.ones(len(parent), bool)
    while changed:
        changed = False
        for i in range(1, len(parent)):
            if alive[i] and not leaf[i] and age[parent[i]] - age[i] < min_internode:
                for j in np.flatnonzero((parent == i) & alive):
                    parent[j] = parent[i]
                alive[i] = False
                changed = True
    # re-index in preorder
    kids = {i: [] for i in range(len(parent))}
    for i in range(1, len(parent)):
        if alive[i]:
            kids[parent[i]].append(i)
    order = []
    st = [0]
    while st:
        i = st.pop()
        order.append(i)
        st.extend(reversed(kids[i]))
    new = {o: k for k, o in enumerate(order)}
    par = np.array([new[parent[o]] if o != 0 else -1 for o in order])
    ag = age[order]
    lf = leaf[order]
    tx = np.array([taxa[o] if lf[k] else -1 for k, o in enumerate(order)])
    lab = [str(tx[k]) if lf[k] else f"N{k}" for k in range(len(order))]
    lo = np.where(lf, 0.0, ag * np.exp(-1.96 * ci_sigma))
    hi = np.where(lf, 0.0, ag * np.exp(1.96 * ci_sigma))
    return PanelTree(parent=par, label=lab, is_leaf=lf, taxid=tx, age=ag, lo=lo, hi=hi)


# ------------------------------------------------------------------------------------------ world
@dataclass
class World:
    panel: PanelTree
    hogs: list                     # hog ids
    species: dict                  # hog -> taxids (S_h)
    b_aa: dict                     # hog -> per-panel-edge AA length (dict child_idx -> b) restricted to the HOG's span
    r_h: dict                      # hog -> 3Di/AA ratio
    lineage_mult: np.ndarray       # (N,) AA lineage rate multipliers (shared, cancel in the ratio under N0)
    spiked: np.ndarray             # (N,) 3Di multiplier per panel edge (1 = none)
    age_beta: float = 0.0
    pairs_true: dict = field(default_factory=dict)  # hog -> (ia, ib, D_aa, D_tdi) on S_h order


def _panel_edge_lengths_for_hog(panel, S, rate, lmult, noise_sd, rng):
    """AA branch lengths on every panel edge below MRCA(S) that leads to S."""
    ht = prune(panel, S)
    b = {}
    for e in range(len(ht.edge_child)):
        c = ht.edge_child[e]
        while c != ht.edge_parent[e]:
            p = panel.parent[c]
            dur = panel.age[p] - panel.age[c]
            b[int(c)] = rate * lmult[c] * dur * np.exp(noise_sd * rng.standard_normal() - noise_sd ** 2 / 2)
            c = p
    return ht, b


def make_world(panel: PanelTree, n_hogs=400, rng=None, min_sp=8, max_sp=30, lineage_sd=0.3, hog_edge_sd=0.2,
               r_median=0.35, r_sd=0.3, max_dist=(0.5, 1.9), age_beta=0.0, spike_frac=0.0, spike_f=1.5,
               spike_min_hogs=0) -> World:
    """b_3Di = r_h * b_AA * exp(age_beta*(log tau - mean_h log tau)) * spike(edge).

    age_beta is on y = log(3Di/AA) (beta, not beta_k). spike multiplies 3Di on a random `spike_frac` of panel
    edges (lineages)."""
    rng = np.random.default_rng(rng)
    N = panel.n_nodes
    lmult = np.exp(lineage_sd * rng.standard_normal(N) - lineage_sd ** 2 / 2)
    internal = np.flatnonzero(~panel.is_leaf)
    nleaf = np.array([bin(panel.clade[i]).count("1") for i in range(N)])
    cand = internal[nleaf[internal] >= min_sp]
    hogs, species, baa, rh = [], {}, {}, {}
    # weight candidates so that young and old roots are both represented (uniform over log root age)
    la = np.log(panel.age[cand])
    wts = 1.0 / np.maximum(1, np.histogram(la, bins=8)[0][np.clip(np.digitize(la, np.histogram(la, bins=8)[1][1:-1]), 0, 7)])
    wts = wts / wts.sum()
    for h in range(n_hogs):
        node = rng.choice(cand, p=wts)
        leaves = [panel.leaf_order[k] for k in range(len(panel.leaf_order)) if (panel.clade[node] >> k) & 1]
        tax = panel.taxid[leaves]
        n = len(tax)
        keep_n = int(min(max_sp, n, max(min_sp, rng.integers(min_sp, max_sp + 1))))
        S = np.sort(rng.choice(tax, size=keep_n, replace=False))
        ht0 = prune(panel, S)
        root_age = panel.age[ht0.root]
        target = rng.uniform(*max_dist)
        rate = target / (2 * root_age)
        hid = f"H{h:05d}"
        ht, b = _panel_edge_lengths_for_hog(panel, S, rate, lmult, hog_edge_sd, rng)
        hogs.append(hid); species[hid] = S; baa[hid] = b
        rh[hid] = r_median * np.exp(r_sd * rng.standard_normal())
    spiked = np.ones(N)
    if spike_frac > 0:
        edges = panel.panel_edge_ids()
        nsp = max(1, int(round(spike_frac * len(edges))))
        spiked[rng.choice(edges, size=nsp, replace=False)] = spike_f
    w = World(panel=panel, hogs=hogs, species=species, b_aa=baa, r_h=rh, lineage_mult=lmult, spiked=spiked,
              age_beta=age_beta)
    for hid in hogs:
        w.pairs_true[hid] = true_pair_distances(w, hid)
    return w


def true_pair_distances(w: World, hid):
    panel = w.panel
    S = w.species[hid]
    ht = prune(panel, S)
    b = w.b_aa[hid]
    pe_aa, pe_3 = {}, {}
    taus = {c: 0.5 * (panel.age[c] + panel.age[panel.parent[c]]) for c in b}
    lt = np.log(np.array(list(taus.values())))
    mlt = lt.mean()
    for c, v in b.items():
        mult = np.exp(w.age_beta * (np.log(taus[c]) - mlt)) if w.age_beta else 1.0
        pe_aa[c] = v
        pe_3[c] = w.r_h[hid] * v * mult * w.spiked[c]
    E = len(ht.edge_child)
    be_aa = np.zeros(E); be_3 = np.zeros(E)
    for e in range(E):
        c = ht.edge_child[e]
        while c != ht.edge_parent[e]:
            be_aa[e] += pe_aa[int(c)]; be_3[e] += pe_3[int(c)]
            c = panel.parent[c]
    n = len(S)
    ia, ib = np.triu_indices(n, 1)
    Xe = (ht.member[ia] ^ ht.member[ib]).astype(float)
    return ia, ib, Xe @ be_aa, Xe @ be_3, be_aa, be_3, ht


# ------------------------------------------------------------------------------------------ observation models
@lru_cache(maxsize=8)
def se_profile(model_name: str, alpha=1.0, L=1000, n=400):
    """Tabulate sqrt(L) * SE(d) of the ML engine on simulated counts (realistic SE shape)."""
    from strucevo.distances import load_model, load_paml, AA_ORDER, simulate_counts, ml_distance
    model = load_model("LG") if model_name == "LG" else load_paml(str(Q3DI_PATH), AA_ORDER)
    grid = np.array([0.01, 0.03, 0.06, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0])
    rng = np.random.default_rng(12345)
    d = np.repeat(grid, n)
    c = simulate_counts(model, d, L, alpha=alpha, rng=rng)
    res = ml_distance(c, model, alpha=alpha, d_max=20)
    sd = np.array([np.nanstd(res.d[d == g]) for g in grid]) * np.sqrt(L)
    return grid, sd


def se_of(model_name, d, L, alpha=1.0):
    g, s = se_profile(model_name, alpha)
    return np.interp(np.maximum(d, 0.0), g, s) / np.sqrt(L) + 1.0 / L


def observe_gauss(w: World, rng=None, L_range=(150, 500), kappa=1.0, rho=0.3, nu_tdi=0.02, tip_sd=0.01,
                  split="tune", hogs=None, c_true=None):
    """Pair table: x = D_aa + e_x, y = D_tdi + u_i + u_j + e_y with corr(e_x, e_y) = rho, Var = kappa^2 SE^2.
    Returned SEs are the *Fisher-like* SEs (true noise = kappa * SE), as in the real tables."""
    rng = np.random.default_rng(rng)
    cols = {k: [] for k in ("hog", "split", "species_a", "species_b", "status", "aa_ml", "aa_ml_se", "aa_sat",
                            "tdi_ml", "tdi_ml_se", "tdi_sat", "aa_ml_seq", "one_to_one", "len_a",
                            "aa_true", "tdi_true")}
    for hid in (hogs or w.hogs):
        ia, ib, Daa, D3, *_ = w.pairs_true[hid]
        S = w.species[hid]
        n = len(ia)
        u = nu_tdi + tip_sd * rng.standard_normal(len(S))
        L = rng.integers(*L_range, size=n)
        if c_true is not None:      # pairwise power-law truth (V3 power): D_3Di = r_h * D_AA^c
            D3 = w.r_h[hid] * Daa ** c_true
        sx = se_of("LG", Daa, L)
        y_true = D3 + u[ia] + u[ib]
        sy = se_of("Q3Di", y_true, L)
        z1 = rng.standard_normal(n)
        z2 = rho * z1 + np.sqrt(1 - rho ** 2) * rng.standard_normal(n)
        x = Daa + kappa * sx * z1
        y = y_true + kappa * sy * z2
        cols["hog"] += [hid] * n; cols["split"] += [split] * n
        cols["species_a"].append(S[ia]); cols["species_b"].append(S[ib]); cols["status"] += ["ok"] * n
        cols["aa_ml"].append(x); cols["aa_ml_se"].append(sx); cols["aa_sat"].append(np.zeros(n, bool))
        cols["tdi_ml"].append(y); cols["tdi_ml_se"].append(sy); cols["tdi_sat"].append(np.zeros(n, bool))
        cols["aa_ml_seq"].append(Daa + 1.2 * sx * rng.standard_normal(n)); cols["one_to_one"].append(np.ones(n, bool))
        cols["len_a"].append(L); cols["aa_true"].append(Daa); cols["tdi_true"].append(y_true)
    return {k: (np.concatenate(v) if isinstance(v[0], np.ndarray) else np.asarray(v)) for k, v in cols.items()}


def observe_counts(w: World, rng=None, L_range=(150, 500), alpha_aa=1.0, alpha_tdi=1.0, nu_tdi=0.02, tip_sd=0.01,
                   split="tune", hogs=None, device=None):
    """Pair table from simulated 20x20 count matrices (LG+G4 for AA, Q.3Di.AF+G4 for 3Di), ML distances with
    d_max = 20, SE floor 1/n_core (trap #3). A separate, independent count draw gives the `aa_ml_seq` proxy."""
    from strucevo.distances import load_model, load_paml, AA_ORDER, simulate_counts, ml_distance
    rng = np.random.default_rng(rng)
    lg = load_model("LG")
    q3 = load_paml(str(Q3DI_PATH), AA_ORDER)
    hs, sa, sb, Ls, Da, D3 = [], [], [], [], [], []
    for hid in (hogs or w.hogs):
        ia, ib, Daa, Dt, *_ = w.pairs_true[hid]
        S = w.species[hid]
        u = nu_tdi + tip_sd * rng.standard_normal(len(S))
        hs += [hid] * len(ia); sa.append(S[ia]); sb.append(S[ib])
        Ls.append(rng.integers(*L_range, size=len(ia))); Da.append(Daa); D3.append(Dt + u[ia] + u[ib])
    sa, sb, L, Da, D3 = map(np.concatenate, (sa, sb, Ls, Da, D3))
    kw = {} if device is None else {"device": device}
    def sim(model, d, alpha, chunk=50000):
        return np.concatenate([simulate_counts(model, d[i:i + chunk], L[i:i + chunk], alpha=alpha, rng=rng)
                               for i in range(0, len(d), chunk)])
    ca = sim(lg, Da, alpha_aa)
    ct = sim(q3, np.maximum(D3, 0.0), alpha_tdi)
    cs = sim(lg, Da, alpha_aa)
    ra = ml_distance(ca, lg, alpha=alpha_aa, d_max=20, **kw)
    rt = ml_distance(ct, q3, alpha=alpha_tdi, d_max=20, **kw)
    rs = ml_distance(cs, lg, alpha=alpha_aa, d_max=20, **kw)
    floor = 1.0 / L
    se_a = np.where(np.isfinite(ra.se), np.maximum(ra.se, floor), floor)
    se_t = np.where(np.isfinite(rt.se), np.maximum(rt.se, floor), floor)
    n = len(L)
    return {"hog": np.asarray(hs), "split": np.full(n, split), "species_a": sa, "species_b": sb,
            "status": np.full(n, "ok"), "aa_ml": ra.d, "aa_ml_se": se_a, "aa_sat": np.asarray(ra.saturated, bool),
            "tdi_ml": rt.d, "tdi_ml_se": se_t, "tdi_sat": np.asarray(rt.saturated, bool),
            "aa_ml_seq": rs.d, "one_to_one": np.ones(n, bool), "len_a": L, "aa_true": Da, "tdi_true": D3,
            "_counts_aa": ca, "_counts_tdi": ct}


def subset(tab: dict, mask) -> dict:
    return {k: v[mask] for k, v in tab.items()}


# ------------------------------------------------------------------------------------------ E6 paralog triplets
def make_triplets(w: World, n=300, rng=None, m_aa=1.0, m_3=1.0, nu_tdi=0.02, tip_sd=0.01, L_range=(150, 500),
                  hogs=None):
    """Synthetic E6 data. Duplication on the stem of internal panel node P (midpoint); copies P1, P2 in species s.
    P2's post-duplication AA rate is x m_aa and its 3Di per AA x m_3 (N0 world: m_aa = m_3 = 1).
    O1 = HOG species in the sister clade(s) of P; O2 = species in the sister clade(s) of parent(P).
    Split-half ML distances come from independent count draws for halves A and B (L/2 columns each)."""
    from strucevo.distances import load_model, load_paml, AA_ORDER, simulate_counts, ml_distance
    rng = np.random.default_rng(rng)
    panel = w.panel
    T = {k: [] for k in ("hog", "dcg_id", "species", "copy1", "copy2", "event_level", "aa_pid_seq", "aa_ml_seq", "o1", "o2")}
    P_ga, P_gb, P_daa, P_d3, P_L = [], [], [], [], []
    hogs = list(hogs or w.hogs)
    tries = 0
    while len(T["hog"]) < n and tries < 50 * n:
        tries += 1
        h = hogs[rng.integers(len(hogs))]
        S = w.species[h]
        ia, ib, Daa, D3, be_aa, be_3, ht = w.pairs_true[h]
        b = w.b_aa[h]
        root = ht.root
        s = S[rng.integers(len(S))]
        leaf = panel.leaf_of_taxid[int(s)]
        anc = []
        c = panel.parent[leaf]
        while c >= 0 and c != root:
            anc.append(c); c = panel.parent[c]
        if not anc:
            continue
        Pn = anc[rng.integers(len(anc))]
        par = panel.parent[Pn]
        if par == root:
            continue       # need a sister of parent(P) inside the HOG clade
        gp = panel.parent[par]
        sp_in = lambda node: [t for t in S if (panel.clade[node] >> panel.leaf_pos[panel.leaf_of_taxid[int(t)]]) & 1]
        o1 = [t for k in panel.children[par] if k != Pn for t in sp_in(k)]
        o2 = [t for k in panel.children[gp] if k != par for t in sp_in(k)]
        if not o1 or not o2:
            continue
        # AA path s -> P plus half of P's stem
        path = 0.0
        c = leaf
        while c != Pn:
            path += b[int(c)]; c = panel.parent[c]
        Ddup = path + 0.5 * b[int(Pn)]
        r = w.r_h[h]
        d_pp_aa = Ddup * (1 + m_aa)
        if not (0.06 <= d_pp_aa <= 2.0):
            continue
        d_pp_3 = r * Ddup * (1 + m_aa * m_3)
        pos = {int(t): k for k, t in enumerate(S)}
        def dpair(t1, t2):
            i, j = sorted((pos[int(t1)], pos[int(t2)]))
            k = np.flatnonzero((ia == i) & (ib == j))[0]
            return Daa[k], D3[k]
        k = len(T["hog"])
        g1, g2 = f"{h}|{s}|c1|{k}", f"{h}|{s}|c2|{k}"
        u = {g1: nu_tdi + tip_sd * rng.standard_normal(), g2: nu_tdi + tip_sd * rng.standard_normal()}
        oset1 = [f"{h}|{o}|o|{k}" for o in o1]
        oset2 = [f"{h}|{o}|o|{k}" for o in o2]
        for o, og in zip(o1 + o2, oset1 + oset2):
            u[og] = nu_tdi + tip_sd * rng.standard_normal()
            da, d3 = dpair(s, o)
            # P1 is the ancestral-rate copy; P2 has the accelerated post-dup segment
            for g, extra_a, extra_3 in ((g1, 0.0, 0.0), (g2, (m_aa - 1) * Ddup, r * (m_aa * m_3 - 1) * Ddup)):
                P_ga.append(g); P_gb.append(og); P_daa.append(da + extra_a); P_d3.append(d3 + extra_3 + u[g] + u[og])
        P_ga.append(g1); P_gb.append(g2); P_daa.append(d_pp_aa); P_d3.append(d_pp_3 + u[g1] + u[g2])
        T["hog"].append(h); T["dcg_id"].append(f"dcg{k}"); T["species"].append(int(s)); T["copy1"].append(g1)
        T["copy2"].append(g2); T["event_level"].append(panel.label[Pn])
        T["aa_ml_seq"].append(d_pp_aa); T["aa_pid_seq"].append(float(np.exp(-d_pp_aa)))
        T["o1"].append(oset1); T["o2"].append(oset2)
    n_p = len(P_ga)
    L = rng.integers(*L_range, size=n_p)
    lg = load_model("LG"); q3 = load_paml(str(Q3DI_PATH), AA_ORDER)
    out = {"gene_a": np.array(P_ga), "gene_b": np.array(P_gb), "status": np.full(n_p, "ok"),
           "aa_true": np.array(P_daa), "tdi_true": np.array(P_d3)}
    for half in ("hA", "hB"):
        ca = simulate_counts(lg, np.array(P_daa), L // 2, alpha=1.0, rng=rng)
        ct = simulate_counts(q3, np.maximum(np.array(P_d3), 0), L // 2, alpha=1.0, rng=rng)
        out[f"aa_ml_{half}"] = ml_distance(ca, lg, alpha=1.0, d_max=20).d
        out[f"tdi_ml_{half}"] = ml_distance(ct, q3, alpha=1.0, d_max=20).d
    T = {k: (np.array(v, dtype=object) if k in ("o1", "o2") else np.asarray(v)) for k, v in T.items()}
    return T, out
