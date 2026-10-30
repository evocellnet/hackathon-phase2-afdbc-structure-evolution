"""Pigeonhole bootstrap (Owen 2007; ANALYSIS_PLAN §2.4, R2-M6, R3-M9).

Three independent multinomial weight vectors per replicate:
    W_h over HOGs, M_s over species, Lambda_e over panel edges  (each ~ Multinomial(n, 1/n)).
Pair weight  = W_h * M_i * M_j.
Edge weight  = W_h * mean(Lambda over the panel edges the HOG edge spans).
No refit of tree topologies. Mixed models: each HOG is replicated W_h times under new group ids.
"""
from __future__ import annotations

import numpy as np


def multinomial_weights(rng: np.random.Generator, n: int, B: int | None = None) -> np.ndarray:
    if B is None:
        return rng.multinomial(n, np.full(n, 1.0 / n)).astype(float)
    return rng.multinomial(n, np.full(n, 1.0 / n), size=B).astype(float)


class Pigeonhole:
    """Factory of pigeonhole replicate weights for fixed index spaces.

    n_hogs, n_species, n_panel_edges: sizes of the three index spaces (pair/edge arrays carry integer codes).
    """

    def __init__(self, n_hogs: int, n_species: int = 0, n_panel_edges: int = 0):
        self.n_hogs, self.n_species, self.n_panel_edges = n_hogs, n_species, n_panel_edges

    def draw(self, rng: np.random.Generator, hogs=True, species=True, edges=True):
        W = multinomial_weights(rng, self.n_hogs) if hogs else np.ones(self.n_hogs)
        M = multinomial_weights(rng, self.n_species) if (species and self.n_species) else np.ones(max(self.n_species, 1))
        L = multinomial_weights(rng, self.n_panel_edges) if (edges and self.n_panel_edges) else np.ones(max(self.n_panel_edges, 1))
        return W, M, L


def pair_weights(W, M, hog_idx, sp_a_idx, sp_b_idx):
    """W_h * M_i * M_j for every pair (index arrays)."""
    return W[hog_idx] * M[sp_a_idx] * M[sp_b_idx]


def edge_weights(W, L, hog_idx, span_ptr, span_edges):
    """W_h * mean Lambda over spanned panel edges.

    span_ptr/span_edges: CSR layout -- the panel edges spanned by HOG edge k are span_edges[span_ptr[k]:span_ptr[k+1]].
    """
    lam = L[span_edges]
    sums = np.add.reduceat(lam, span_ptr[:-1]) if len(lam) else np.zeros(len(span_ptr) - 1)
    cnt = np.diff(span_ptr)
    return W[hog_idx] * sums / np.maximum(cnt, 1)


def replicate_groups(W, hog_codes):
    """Mixed-model replication: returns (row_index, new_group_id) for HOG rows repeated W_h times.

    hog_codes: (n_rows,) integer HOG code of each data row. A HOG drawn k times appears k times with k distinct ids.
    """
    W = np.asarray(W).astype(int)
    order = np.argsort(hog_codes, kind="stable")
    codes_sorted = hog_codes[order]
    starts = np.searchsorted(codes_sorted, np.arange(len(W)), "left")
    ends = np.searchsorted(codes_sorted, np.arange(len(W)), "right")
    rows, gids = [], []
    gid = 0
    for h in np.flatnonzero(W > 0):
        r = order[starts[h]:ends[h]]
        for _ in range(W[h]):
            rows.append(r)
            gids.append(np.full(len(r), gid))
            gid += 1
    if not rows:
        return np.zeros(0, int), np.zeros(0, int)
    return np.concatenate(rows), np.concatenate(gids)
