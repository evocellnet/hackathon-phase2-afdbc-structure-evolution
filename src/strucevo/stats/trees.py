"""Panel-tree utilities: pruning to S_h with unary-node suppression, HOG-edge -> panel mapping,
path-incidence design matrices, neighbour joining, node-age draws (ANALYSIS_PLAN §2.2, §2.5, §4.4, R3-M4).

The panel tree is held as flat arrays (preorder) so that pruning thousands of HOGs is cheap; dendropy is
used only for I/O.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class PanelTree:
    parent: np.ndarray                 # (N,) int, -1 for root; preorder (parent index < child index)
    label: list                        # node ids (str); leaves = taxid as str
    is_leaf: np.ndarray                # (N,) bool
    taxid: np.ndarray                  # (N,) int (-1 for internal)
    age: np.ndarray                    # (N,) float median age (Mya); leaves 0; NaN if unknown
    lo: np.ndarray | None = None       # 95% CI bounds (Mya) for age draws
    hi: np.ndarray | None = None
    children: list = field(default_factory=list)
    depth: np.ndarray | None = None
    leaf_of_taxid: dict = field(default_factory=dict)
    clade: list = field(default_factory=list)   # python int bitmask over leaf order

    def __post_init__(self):
        N = len(self.parent)
        self.children = [[] for _ in range(N)]
        for i, p in enumerate(self.parent):
            if p >= 0:
                self.children[p].append(i)
        self.depth = np.zeros(N, int)
        for i in range(N):
            if self.parent[i] >= 0:
                self.depth[i] = self.depth[self.parent[i]] + 1
        leaves = np.flatnonzero(self.is_leaf)
        self.leaf_order = leaves
        self.leaf_pos = {int(n): k for k, n in enumerate(leaves)}
        self.leaf_of_taxid = {int(self.taxid[n]): int(n) for n in leaves}
        clade = [0] * N
        for i in range(N - 1, -1, -1):
            if self.is_leaf[i]:
                clade[i] = 1 << self.leaf_pos[i]
            else:
                m = 0
                for c in self.children[i]:
                    m |= clade[c]
                clade[i] = m
        self.clade = clade
        self.index = {l: i for i, l in enumerate(self.label)}

    @property
    def n_nodes(self):
        return len(self.parent)

    # ------------------------------------------------------------------ construction
    @classmethod
    def from_dendropy(cls, tree, ages: dict | None = None):
        """tree: dendropy.Tree; leaf taxon labels must be taxids. Internal node labels are used as ids
        (generated 'n{k}' if absent). ages: {node_id: (median, lo95, hi95)} (leaves default to 0)."""
        nodes = list(tree.preorder_node_iter())
        idx = {id(n): k for k, n in enumerate(nodes)}
        parent = np.array([idx[id(n.parent_node)] if n.parent_node is not None else -1 for n in nodes])
        is_leaf = np.array([n.is_leaf() for n in nodes])
        taxid = np.array([int(n.taxon.label) if n.is_leaf() else -1 for n in nodes])
        label = []
        for k, n in enumerate(nodes):
            if n.is_leaf():
                label.append(str(int(n.taxon.label)))
            else:
                label.append(str(n.label) if n.label not in (None, "") else f"n{k}")
        N = len(nodes)
        age = np.full(N, np.nan); lo = np.full(N, np.nan); hi = np.full(N, np.nan)
        age[is_leaf] = lo[is_leaf] = hi[is_leaf] = 0.0
        if ages:
            for k, l in enumerate(label):
                if l in ages and not is_leaf[k]:
                    age[k], lo[k], hi[k] = ages[l]
        return cls(parent=parent, label=label, is_leaf=is_leaf, taxid=taxid, age=age, lo=lo, hi=hi)

    @classmethod
    def from_newick(cls, path_or_str, ages: dict | None = None, is_path=True):
        import dendropy
        kw = dict(schema="newick", preserve_underscores=True, suppress_internal_node_taxa=True,
                  suppress_leaf_node_taxa=False)
        t = dendropy.Tree.get(path=path_or_str, **kw) if is_path else dendropy.Tree.get(data=path_or_str, **kw)
        return cls.from_dendropy(t, ages)

    def with_ages(self, age, lo=None, hi=None):
        import copy
        t = copy.copy(self)
        t.age = np.asarray(age, float)
        if lo is not None:
            t.lo = np.asarray(lo, float)
        if hi is not None:
            t.hi = np.asarray(hi, float)
        return t

    def panel_edge_ids(self):
        """Panel edges are identified by their child node index (all non-root nodes)."""
        return np.flatnonzero(self.parent >= 0)


@dataclass
class HogTree:
    """T_h: panel tree pruned to S_h with unary nodes suppressed."""
    species: np.ndarray            # (n,) taxids in the order used for the membership matrix
    edge_child: np.ndarray         # (E,) panel node index of each HOG-edge child (= panel MRCA of child clade)
    edge_parent: np.ndarray        # (E,) panel node index of the HOG-edge parent
    n_panel_edges: np.ndarray      # (E,)
    is_terminal: np.ndarray        # (E,) bool
    member: np.ndarray             # (n, E) bool: species k in child clade of edge e
    root: int                      # panel index of MRCA(S_h)
    param_of_edge: np.ndarray      # (E,) index of the WLS parameter that carries edge e
    share: np.ndarray              # (E,) fraction of the parameter apportioned to edge e (1 except merged root edges)
    n_params: int
    root_merged: bool

    def edge_ids(self, panel: PanelTree):
        out = []
        for e in range(len(self.edge_child)):
            sp = np.sort(self.species[self.member[:, e]])
            out.append(",".join(map(str, sp)))
        return out

    def design(self, ia, ib):
        return design_matrix(self, ia, ib)


def prune(panel: PanelTree, species) -> HogTree:
    """Induced subtree of the panel on `species` (taxids) with unary nodes suppressed.

    If the induced root is bifurcating, its two child edges are not separately identifiable from pairwise
    distances: they are carried by one WLS parameter and apportioned by panel time span
    (t_root - t_child), or 50/50 when ages are unknown (DEVIATION D1 in the stats report).
    """
    species = np.asarray(species, dtype=np.int64)
    leaves = [panel.leaf_of_taxid[int(s)] for s in species]
    mask = 0
    for l in leaves:
        mask |= panel.clade[l]
    N = panel.n_nodes
    restricted = [panel.clade[i] & mask for i in range(N)]
    keep = np.zeros(N, bool)
    for i in range(N):
        if restricted[i] == 0:
            continue
        if panel.is_leaf[i]:
            keep[i] = True
        else:
            nz = sum(1 for c in panel.children[i] if restricted[c] != 0)
            keep[i] = nz >= 2
    kept = np.flatnonzero(keep)
    # nearest kept ancestor
    kpar = {}
    for i in kept:
        p = panel.parent[i]
        while p >= 0 and not keep[p]:
            p = panel.parent[p]
        kpar[int(i)] = int(p)
    # the induced root is the kept node whose nearest kept ancestor does not exist
    root = int(min(kept, key=lambda i: panel.depth[i]))
    edges = [int(i) for i in kept if i != root]
    E = len(edges)
    edge_child = np.array(edges, dtype=np.int64)
    edge_parent = np.array([kpar[i] for i in edges], dtype=np.int64)
    n_pe = panel.depth[edge_child] - panel.depth[edge_parent]
    is_term = panel.is_leaf[edge_child]
    sp_pos = {panel.leaf_pos[l]: k for k, l in enumerate(leaves)}
    member = np.zeros((len(species), E), bool)
    for e, c in enumerate(edges):
        m = restricted[c]
        while m:
            low = m & -m
            member[sp_pos[low.bit_length() - 1], e] = True
            m ^= low
    root_children = np.flatnonzero(edge_parent == root)
    param = np.arange(E)
    share = np.ones(E)
    merged = False
    if len(root_children) == 2:
        merged = True
        e1, e2 = root_children
        t_r = panel.age[root]
        s1 = t_r - panel.age[edge_child[e1]]
        s2 = t_r - panel.age[edge_child[e2]]
        if np.isfinite(s1) and np.isfinite(s2) and (s1 + s2) > 0:
            f1 = s1 / (s1 + s2)
        else:
            f1 = 0.5
        share[e1], share[e2] = f1, 1 - f1
        # renumber parameters: e2 shares e1's parameter
        param = np.empty(E, int)
        k = 0
        for e in range(E):
            if e == e2:
                continue
            param[e] = k
            k += 1
        param[e2] = param[e1]
        n_params = E - 1
    else:
        n_params = E
    return HogTree(species=species, edge_child=edge_child, edge_parent=edge_parent, n_panel_edges=n_pe,
                   is_terminal=is_term, member=member, root=root, param_of_edge=param, share=share,
                   n_params=n_params, root_merged=merged)


def design_matrix(ht: HogTree, ia, ib) -> np.ndarray:
    """Pairs x params path-incidence matrix (merged root edges share one column)."""
    Xe = (ht.member[ia] ^ ht.member[ib]).astype(float)
    X = np.zeros((len(ia), ht.n_params))
    np.add.at(X.T, ht.param_of_edge, Xe.T)
    return np.minimum(X, 1.0)  # merged root columns are identical -> OR


# ------------------------------------------------------------------------------------------ neighbour joining
def nj_splits(D: np.ndarray) -> np.ndarray:
    """Neighbour joining (Saitou & Nei 1987; same criterion as dendropy's PhylogeneticDistanceMatrix.nj_tree).
    Returns the (n, E) membership matrix of the unrooted NJ tree's edges (one side of every split, 2n-3 edges).
    """
    D = np.array(D, float)
    n = D.shape[0]
    clusters = [1 << i for i in range(n)]
    active = list(range(n))
    splits = [1 << i for i in range(n)]          # terminal edges
    while len(active) > 3:
        m = len(active)
        Dm = D[np.ix_(active, active)]
        r = Dm.sum(1)
        Qm = (m - 2) * Dm - r[:, None] - r[None, :]
        np.fill_diagonal(Qm, np.inf)
        a, b = np.unravel_index(np.argmin(Qm), Qm.shape)
        i, j = active[a], active[b]
        newd = 0.5 * (D[i] + D[j] - D[i, j])
        D[i, :] = newd
        D[:, i] = newd
        D[i, i] = 0.0
        clusters[i] = clusters[i] | clusters[j]
        active.remove(j)
        splits.append(clusters[i])
    # remaining 3 clusters joined at a trifurcation; the internal splits recorded so far are complete
    full = (1 << n) - 1
    out = []
    seen = set()
    for s in splits:
        key = min(s, full ^ s)
        if key in seen or s == full or s == 0:
            continue
        seen.add(key)
        out.append(s)
    M = np.zeros((n, len(out)), bool)
    for e, s in enumerate(out):
        for k in range(n):
            M[k, e] = (s >> k) & 1
    return M


# ------------------------------------------------------------------------------------------ ages
def edge_ages(panel: PanelTree, ht: HogTree, age=None):
    """t_child, t_parent, tau (arithmetic midpoint; R3-M5) for each HOG edge from panel node ages."""
    age = panel.age if age is None else age
    tc = age[ht.edge_child]
    tp = age[ht.edge_parent]
    return tc, tp, 0.5 * (tc + tp)


def draw_ages(panel: PanelTree, rng: np.random.Generator, repair=True) -> np.ndarray:
    """One fresh set of node ages: lognormal around the median with sigma=(ln hi - ln lo)/3.92, then a
    top-down repair child = min(child, parent - 0.1) (ANALYSIS_PLAN §4.4). Leaves stay 0."""
    med = panel.age.copy()
    lo, hi = panel.lo, panel.hi
    out = med.copy()
    internal = ~panel.is_leaf & np.isfinite(med) & (med > 0)
    sig = np.zeros_like(med)
    ok = internal & np.isfinite(lo) & np.isfinite(hi) & (lo > 0) & (hi > 0)
    sig[ok] = (np.log(hi[ok]) - np.log(lo[ok])) / 3.92
    z = rng.standard_normal(len(med))
    out[internal] = med[internal] * np.exp(sig[internal] * z[internal])
    if repair:
        for i in range(len(out)):           # preorder: parents first
            p = panel.parent[i]
            if p >= 0 and internal[i]:
                out[i] = min(out[i], out[p] - 0.1)
    return out
