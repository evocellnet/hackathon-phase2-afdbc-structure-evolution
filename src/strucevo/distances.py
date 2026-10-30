"""Pairwise divergence measures.

Sequence-level (amino acid or 3Di, any 20-letter alphabet)
    encode / pair_counts            aligned pairs -> (N, 20, 20) count matrices of aligned states
                                    (only columns where both are valid letters, optional per-column mask)
    simple_distances                p-distance, Poisson, Kimura (1983) protein correction from counts
    p_distance / poisson_distance / kimura_distance   single-pair convenience wrappers
    RateModel, load_paml, load_model                    reversible 20-state models (PAML .dat)
    ml_distance                     batched ML distance under model + discrete Gamma (+I), PyTorch/GPU,
                                    with observed-information SE and saturation flags
    estimate_gtr                    estimate a GTR-20 model (e.g. for 3Di) from count matrices of close pairs
    simulate_counts                 simulate count matrices under a model (tests / benchmarks)
Structure-level
    one_minus, neg_log              TM-score / LDDT transforms
    kabsch_rmsd                     batched, masked Kabsch superposition RMSD (torch)

Conventions
    Distances d are expected substitutions per site; Q is normalised to mean rate 1 (sum_i pi_i q_i = 1) and
    the discrete-Gamma category rates (Yang 1994, category means, as IQ-TREE/AliSim default) have mean 1.
    With +I, variable-site rates are scaled by 1/(1-pinv) so the overall mean rate is 1 (IQ-TREE convention).
"""
from __future__ import annotations

import math
import os
from dataclasses import dataclass, field

import numpy as np
import torch
from scipy import special, stats

AA_ORDER = "ARNDCQEGHILKMFPSTWYV"        # PAML / IQ-TREE amino-acid order
THREEDI_ORDER = "ACDEFGHIKLMNPQRSTVWY"   # Foldseek 3Di letters (alphabetical, as Foldseek writes them)
_DATA_DIR = os.path.join(os.path.dirname(__file__), "data")


# ----------------------------------------------------------------------------------------------------------
# Encoding and count matrices
# ----------------------------------------------------------------------------------------------------------
def _lut(alphabet: str) -> np.ndarray:
    lut = np.full(256, -1, dtype=np.int16)
    for i, c in enumerate(alphabet):
        lut[ord(c.upper())] = i
        lut[ord(c.lower())] = i
    return lut


def encode(seq: str, alphabet: str = AA_ORDER) -> np.ndarray:
    """Encode a string to state indices; gaps / unknown letters (-, ., X, B, Z, *, ...) -> -1."""
    b = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
    return _lut(alphabet)[b]


def pair_counts(seqs_a, seqs_b, alphabet: str = AA_ORDER, masks=None) -> np.ndarray:
    """Count matrices for aligned pairs.

    seqs_a, seqs_b : sequences of equal-length aligned strings (gapped), or a single pair of strings.
    masks          : optional per-pair boolean arrays (True = use column), e.g. pLDDT>=70 in both.
    Returns int64 array (N, K, K) with N[i, j] = #columns where A has state i and B has state j.
    """
    if isinstance(seqs_a, str):
        seqs_a, seqs_b = [seqs_a], [seqs_b]
        masks = None if masks is None else [masks]
    K = len(alphabet)
    lut = _lut(alphabet)
    idx_all = []
    for n, (a, b) in enumerate(zip(seqs_a, seqs_b)):
        if len(a) != len(b):
            raise ValueError(f"pair {n}: aligned sequences differ in length ({len(a)} vs {len(b)})")
        ea = lut[np.frombuffer(a.encode("ascii"), dtype=np.uint8)].astype(np.int64)
        eb = lut[np.frombuffer(b.encode("ascii"), dtype=np.uint8)].astype(np.int64)
        ok = (ea >= 0) & (eb >= 0)
        if masks is not None and masks[n] is not None:
            m = np.asarray(masks[n], dtype=bool)
            if m.shape[0] != len(a):
                raise ValueError(f"pair {n}: mask length {m.shape[0]} != alignment length {len(a)}")
            ok &= m
        idx_all.append(n * K * K + ea[ok] * K + eb[ok])
    N = len(idx_all)
    flat = np.concatenate(idx_all) if idx_all else np.zeros(0, dtype=np.int64)
    return np.bincount(flat, minlength=N * K * K).reshape(N, K, K)


def plddt_mask(plddt_a, plddt_b, threshold: float = 70.0) -> np.ndarray:
    """Column mask for an aligned pair: both residues have pLDDT >= threshold (NaN at gaps -> False)."""
    a = np.asarray(plddt_a, dtype=float)
    b = np.asarray(plddt_b, dtype=float)
    with np.errstate(invalid="ignore"):
        return (a >= threshold) & (b >= threshold)


# ----------------------------------------------------------------------------------------------------------
# Simple distances
# ----------------------------------------------------------------------------------------------------------
def simple_distances(counts) -> dict:
    """p, Poisson and Kimura (1983) protein distances from (N,K,K) counts.

    Kimura: d = -ln(1 - p - 0.2 p^2); undefined (NaN, flagged) when the argument is <= 0 (p >~ 0.854).
    Poisson: d = -ln(1 - p); infinite when p = 1 (NaN, flagged).
    """
    c = np.asarray(counts, dtype=float)
    if c.ndim == 2:
        c = c[None]
    n = c.sum(axis=(1, 2))
    ident = np.trace(c, axis1=1, axis2=2)
    with np.errstate(divide="ignore", invalid="ignore"):
        p = 1.0 - ident / n
        arg_p = 1.0 - p
        arg_k = 1.0 - p - 0.2 * p * p
        poisson = np.where(arg_p > 0, -np.log(np.where(arg_p > 0, arg_p, 1.0)), np.nan)
        kimura = np.where(arg_k > 0, -np.log(np.where(arg_k > 0, arg_k, 1.0)), np.nan)
    p = np.where(n > 0, p, np.nan)
    return dict(n_sites=n.astype(np.int64), p=p, poisson=poisson, kimura=kimura,
                poisson_saturated=~np.isfinite(poisson), kimura_saturated=~np.isfinite(kimura))


def _single(fn_key, a, b, mask=None, alphabet=AA_ORDER):
    return float(simple_distances(pair_counts(a, b, alphabet, masks=mask))[fn_key][0])


def p_distance(a: str, b: str, mask=None, alphabet: str = AA_ORDER) -> float:
    return _single("p", a, b, mask, alphabet)


def poisson_distance(a: str, b: str, mask=None, alphabet: str = AA_ORDER) -> float:
    return _single("poisson", a, b, mask, alphabet)


def kimura_distance(a: str, b: str, mask=None, alphabet: str = AA_ORDER) -> float:
    return _single("kimura", a, b, mask, alphabet)


# ----------------------------------------------------------------------------------------------------------
# Rate models
# ----------------------------------------------------------------------------------------------------------
@dataclass
class RateModel:
    """Reversible K-state model: symmetric exchangeabilities S (K,K) and stationary frequencies pi (K,)."""
    S: np.ndarray
    pi: np.ndarray
    alphabet: str
    name: str = "custom"
    source: str = ""
    meta: dict = field(default_factory=dict)

    def __post_init__(self):
        S = np.asarray(self.S, dtype=float)
        S = 0.5 * (S + S.T)
        np.fill_diagonal(S, 0.0)
        pi = np.asarray(self.pi, dtype=float)
        if np.any(pi <= 0):
            raise ValueError("all frequencies must be > 0")
        self.pi = pi / pi.sum()
        if np.any(S < 0):
            raise ValueError("exchangeabilities must be >= 0")
        self.S = S
        if len(self.alphabet) != len(self.pi) or S.shape != (len(self.pi),) * 2:
            raise ValueError("alphabet / S / pi dimension mismatch")

    @property
    def K(self) -> int:
        return len(self.pi)

    def Q(self) -> np.ndarray:
        """Rate matrix normalised to one expected substitution per unit time."""
        Q = self.S * self.pi[None, :]
        np.fill_diagonal(Q, -Q.sum(axis=1))
        return Q / (-(self.pi * np.diag(Q)).sum())

    def normalized_S(self) -> np.ndarray:
        """S scaled consistently with the normalised Q."""
        Q = self.S * self.pi[None, :]
        return self.S / (Q.sum(axis=1) * self.pi).sum()

    def P(self, d: float) -> np.ndarray:
        from scipy.linalg import expm
        return expm(self.Q() * d)

    def reorder(self, alphabet: str) -> "RateModel":
        """Return the same model with states permuted to a new letter order."""
        perm = [self.alphabet.index(c) for c in alphabet]
        return RateModel(self.S[np.ix_(perm, perm)], self.pi[perm], alphabet, self.name, self.source, dict(self.meta))

    def to_paml(self, path: str, comment: str = ""):
        with open(path, "w") as fh:
            S = self.normalized_S()
            for i in range(1, self.K):
                fh.write(" ".join(f"{S[i, j]:.6f}" for j in range(i)) + "\n")
            fh.write("\n" + " ".join(f"{p:.6f}" for p in self.pi) + "\n\n")
            fh.write(" ".join(self.alphabet) + "\n")
            fh.write(f"\n{self.name} {comment}\n{self.source}\n")


def load_paml(path: str, alphabet: str = AA_ORDER, name: str | None = None) -> RateModel:
    """Parse a PAML-format .dat file: K(K-1)/2 lower-triangle exchangeabilities then K frequencies.

    Numbers are read in order until K(K-1)/2 + K have been collected; any trailing text is ignored.
    The alphabet must be given in the order used by the file (PAML amino-acid files: AA_ORDER).
    """
    K = len(alphabet)
    need = K * (K - 1) // 2 + K
    vals = []
    with open(path) as fh:
        for line in fh:
            for tok in line.replace(",", " ").split():
                try:
                    vals.append(float(tok))
                except ValueError:
                    if len(vals) > 0 and len(vals) < need:
                        raise ValueError(f"{path}: non-numeric token {tok!r} after {len(vals)} numbers")
                    continue
                if len(vals) == need:
                    break
            if len(vals) == need:
                break
    if len(vals) < need:
        raise ValueError(f"{path}: found {len(vals)} numbers, need {need}")
    S = np.zeros((K, K))
    it = iter(vals[: K * (K - 1) // 2])
    for i in range(1, K):
        for j in range(i):
            S[i, j] = S[j, i] = next(it)
    pi = np.array(vals[K * (K - 1) // 2:])
    if abs(pi.sum() - 1) > 1e-3:
        raise ValueError(f"{path}: frequencies sum to {pi.sum():.4f}")
    return RateModel(S, pi, alphabet, name or os.path.splitext(os.path.basename(path))[0].upper(), path)


def load_model(name: str) -> RateModel:
    """Bundled models: 'LG', 'WAG' (see data/PROVENANCE.txt)."""
    path = os.path.join(_DATA_DIR, f"{name.lower()}.dat")
    if not os.path.exists(path):
        raise KeyError(f"no bundled model {name!r}; use load_paml(path, alphabet)")
    return load_paml(path, AA_ORDER, name.upper())


def gamma_rates(alpha: float | None, ncat: int = 4, pinv: float = 0.0):
    """Discrete-Gamma category rates (means; Yang 1994) and weights, with optional invariant-site proportion.

    Returns (rates, weights) for the variable-site categories; rates are scaled by 1/(1-pinv).
    """
    if alpha is None or ncat <= 1:
        r = np.ones(1)
    else:
        q = stats.gamma.ppf(np.arange(1, ncat) / ncat, a=alpha, scale=1.0 / alpha)
        cdf1 = np.concatenate([[0.0], special.gammainc(alpha + 1.0, alpha * q), [1.0]])
        r = ncat * np.diff(cdf1)
    w = np.full(len(r), (1.0 - pinv) / len(r))
    return r / (1.0 - pinv), w


# ----------------------------------------------------------------------------------------------------------
# Batched ML distance
# ----------------------------------------------------------------------------------------------------------
@dataclass
class MLResult:
    d: np.ndarray            # ML distance (capped at d_max); NaN if no sites
    se: np.ndarray           # sqrt(1 / observed information); NaN at boundaries (d=0 or d=d_max)
    saturated: np.ndarray    # sat_cap | sat_identity | (lr_unrelated < lr_crit): distance not informative
    sat_cap: np.ndarray      # d hit d_max
    sat_identity: np.ndarray # identity <= expected identity of a random composition-matched pair
    lr_unrelated: np.ndarray # 2[L(d_hat) - L(d=inf)], L(inf) = sum N_ij log(pi_i pi_j): LR vs 'unrelated'
    at_zero: np.ndarray      # no differences -> d = 0
    n_sites: np.ndarray
    identity: np.ndarray     # fraction identical among counted sites
    random_identity: np.ndarray  # sum_i fA_i fB_i (marginal compositions of the pair)
    loglik: np.ndarray       # sum_ij N_ij log(pi_i P_ij(d))
    converged: np.ndarray

    def as_dict(self):
        return {k: getattr(self, k) for k in self.__dataclass_fields__}


class _Engine:
    """Eigendecomposition of Q once; evaluates P, dP/dd, d2P/dd2 for batches of distances."""

    def __init__(self, model: RateModel, alpha, ncat, pinv, device, dtype):
        Q = model.Q()
        pi = model.pi
        sq = np.sqrt(pi)
        A = (sq[:, None] * Q) / sq[None, :]
        A = 0.5 * (A + A.T)
        lam, U = np.linalg.eigh(A)
        V = U / sq[:, None]               # Pi^{-1/2} U
        W = U.T * sq[None, :]             # U^T Pi^{1/2}
        r, w = gamma_rates(alpha, ncat, pinv)
        t = lambda x: torch.as_tensor(x, device=device, dtype=dtype)
        self.V, self.W, self.lam = t(V), t(W), t(lam)
        self.r, self.w = t(r), t(w)
        self.pinv = float(pinv)
        self.K = model.K
        self.logpi = t(np.log(pi))
        self.eye = torch.eye(model.K, device=device, dtype=dtype)

    def _e(self, d, order=2):
        # d: (n,) -> e^{(o)}_m(d) = sum_k w_k (lam_m r_k)^o exp(lam_m r_k d), shapes (n, K)
        lr = self.lam[None, :] * self.r[:, None]                    # (C, K)
        ex = torch.exp(d[:, None, None] * lr[None])                   # (n, C, K)
        out = [torch.einsum("c,nck->nk", self.w, ex)]
        if order >= 1:
            out.append(torch.einsum("c,nck->nk", self.w, ex * lr[None]))
        if order >= 2:
            out.append(torch.einsum("c,nck->nk", self.w, ex * lr[None] ** 2))
        return out

    def mats(self, d, order=2):
        es = self._e(d, order)
        Ms = [torch.einsum("im,nm,mj->nij", self.V, e, self.W) for e in es]
        Ms[0] = Ms[0] + self.pinv * self.eye
        return Ms

    def logP_grid(self, grid):
        P = self.mats(grid, order=0)[0]
        return torch.log(P.clamp_min(1e-300))


def ml_distance(counts, model: RateModel, alpha: float | None = 1.0, ncat: int = 4, pinv: float = 0.0,
                d_max: float = 10.0, device: str | None = None, chunk: int = 200_000, n_grid: int = 400,
                max_iter: int = 40, tol: float = 1e-7, lr_crit: float = 2.706, dtype=torch.float64) -> MLResult:
    """Batched ML pairwise distance.

    Maximises L(d) = sum_ij N_ij log(pi_i P_ij(d)), P(d) = pinv*I + sum_k w_k exp(Q r_k d), for every pair.
    Algorithm: log-spaced grid (single matmul) -> bracket -> safeguarded Newton (bisection fallback).

    counts : (N,K,K) array/tensor of aligned-state counts (rows = sequence A, columns = sequence B).
    lr_crit: pairs whose likelihood-ratio statistic against 'unrelated' (d = infinity, P_ij = pi_j) is below
             lr_crit are flagged saturated (default 2.706 = 5% level for the boundary 50:50 chi2_0/chi2_1 mixture).
    Returns MLResult of numpy arrays.
    """
    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"
    C_all = torch.as_tensor(np.asarray(counts) if not torch.is_tensor(counts) else counts)
    if C_all.ndim == 2:
        C_all = C_all[None]
    K = model.K
    if C_all.shape[1:] != (K, K):
        raise ValueError(f"counts must be (N,{K},{K})")
    eng = _Engine(model, alpha, ncat, pinv, device, dtype)
    grid = torch.cat([torch.zeros(1, dtype=dtype, device=device),
                      torch.logspace(-5, math.log10(d_max), n_grid - 1, dtype=dtype, device=device)])
    grid[-1] = d_max
    logPg = eng.logP_grid(grid).reshape(len(grid), K * K)                # (G, K*K)
    G = len(grid)
    outs = []
    for s in range(0, C_all.shape[0], chunk):
        C = C_all[s:s + chunk].to(device=device, dtype=dtype)
        n = C.shape[0]
        Cf = C.reshape(n, K * K)
        nsites = Cf.sum(1)
        ident = torch.diagonal(C, dim1=1, dim2=2).sum(1)
        LLg = Cf @ logPg.T                                                 # (n, G)
        gi = LLg.argmax(1)
        d = grid[gi].clone()
        lo = grid[(gi - 1).clamp_min(0)].clone()
        hi = grid[(gi + 1).clamp_max(G - 1)].clone()
        offdiag = nsites - ident
        zero = offdiag <= 0
        top = gi == G - 1
        active = ~zero & ~top & (nsites > 0)
        conv = ~active.clone()
        for _ in range(max_iter):
            idx = torch.nonzero(active & ~conv).squeeze(1)
            if idx.numel() == 0:
                break
            di, Ci = d[idx], C[idx]
            P, P1, P2 = eng.mats(di)
            P = P.clamp_min(1e-300)
            R1 = P1 / P
            g = (Ci * R1).sum((1, 2))
            h = (Ci * (P2 / P - R1 * R1)).sum((1, 2))
            loi = torch.where(g > 0, di, lo[idx])
            hii = torch.where(g > 0, hi[idx], di)
            step = torch.where(h < 0, -g / h, torch.zeros_like(g))
            dn = di + step
            bad = (h >= 0) | (dn <= loi) | (dn >= hii)
            dn = torch.where(bad, 0.5 * (loi + hii), dn)
            lo[idx], hi[idx], d[idx] = loi, hii, dn
            done = ((dn - di).abs() <= tol * (1 + di)) | ((hii - loi) <= tol * (1 + di))
            conv[idx] = done
        d = torch.where(zero & (nsites > 0), torch.zeros_like(d), d)
        d = torch.where(top, torch.full_like(d, d_max), d)
        # final quantities
        P, P1, P2 = eng.mats(d)
        P = P.clamp_min(1e-300)
        R1 = P1 / P
        h = (C * (P2 / P - R1 * R1)).sum((1, 2))
        ll = (C * torch.log(P)).sum((1, 2)) + (C.sum(2) * eng.logpi[None]).sum(1)
        interior = active & (d < d_max * (1 - 1e-9)) & (d > 0)
        se = torch.where(interior & (h < 0), 1.0 / torch.sqrt((-h).clamp_min(1e-300)),
                         torch.full_like(d, float("nan")))
        fa = C.sum(2) / nsites.clamp_min(1)[:, None]
        fb = C.sum(1) / nsites.clamp_min(1)[:, None]
        rid = (fa * fb).sum(1)
        idf = ident / nsites.clamp_min(1)
        ll_inf = (C.sum(2) * eng.logpi[None]).sum(1) + (C.sum(1) * eng.logpi[None]).sum(1)
        lr = 2.0 * (ll - ll_inf)
        sat_cap = d >= d_max * (1 - 1e-9)
        sat_id = (idf <= rid) & ~zero
        nos = nsites <= 0
        sat = sat_cap | sat_id | (lr < lr_crit) | nos
        d = torch.where(nos, torch.full_like(d, float("nan")), d)
        outs.append([d, se, sat, sat_cap, sat_id, lr, zero & ~nos, nsites, idf, rid, ll, conv])
    cat = [torch.cat([o[i] for o in outs]).cpu().numpy() for i in range(12)]
    return MLResult(d=cat[0], se=cat[1], saturated=cat[2].astype(bool), sat_cap=cat[3].astype(bool),
                    sat_identity=cat[4].astype(bool), lr_unrelated=cat[5], at_zero=cat[6].astype(bool),
                    n_sites=cat[7].astype(np.int64), identity=cat[8], random_identity=cat[9],
                    loglik=cat[10], converged=cat[11].astype(bool))


def profile_alpha(counts, model: RateModel, alphas=(0.25, 0.35, 0.5, 0.7, 1.0, 1.4, 2.0, 3.0, 5.0, 1e3),
                  ncat: int = 4, **kw):
    """Total log-likelihood over pairs (d profiled per pair) for a grid of Gamma shapes.

    Returns (alphas, total_loglik). NB: with only two sequences, alpha is weakly identified; use many pairs.
    """
    lls = []
    for a in alphas:
        res = ml_distance(counts, model, alpha=a, ncat=ncat, **kw)
        ok = np.isfinite(res.loglik)
        lls.append(res.loglik[ok].sum())
    return np.asarray(alphas, float), np.asarray(lls)


# ----------------------------------------------------------------------------------------------------------
# Simulation (for tests / benchmarks)
# ----------------------------------------------------------------------------------------------------------
def simulate_counts(model: RateModel, d, length, alpha: float | None = 1.0, ncat: int = 4, pinv: float = 0.0,
                    rng=None) -> np.ndarray:
    """Simulate (N,K,K) count matrices for pairs at distances d (array) with `length` sites each."""
    from scipy.linalg import expm
    rng = np.random.default_rng(rng)
    d = np.atleast_1d(np.asarray(d, float))
    length = np.broadcast_to(np.asarray(length, dtype=np.int64), d.shape)
    K = model.K
    r, w = gamma_rates(alpha, ncat, pinv)
    cats_r = np.concatenate([[0.0], r]) if pinv > 0 else r
    cats_w = np.concatenate([[pinv], w]) if pinv > 0 else w
    Q = model.Q()
    lam, U = np.linalg.eigh((np.sqrt(model.pi)[:, None] * Q) / np.sqrt(model.pi)[None, :])
    V = U / np.sqrt(model.pi)[:, None]
    Wm = U.T * np.sqrt(model.pi)[None, :]
    nc = rng.multinomial(length, cats_w / cats_w.sum())          # (N, C)
    out = np.zeros((len(d), K * K), dtype=np.int64)
    for c, rc in enumerate(cats_r):
        P = np.einsum("im,nm,mj->nij", V, np.exp(lam[None] * rc * d[:, None]), Wm)
        J = np.clip(model.pi[None, :, None] * P, 0, None).reshape(len(d), K * K)
        J /= J.sum(1, keepdims=True)
        out += rng.multinomial(nc[:, c], J)
    return out.reshape(len(d), K, K)


# ----------------------------------------------------------------------------------------------------------
# GTR estimation (e.g. 3Di matrix from close ortholog pairs)
# ----------------------------------------------------------------------------------------------------------
def counting_gtr(counts, alphabet: str, pseudocount: float = 0.5, name: str = "GTR-count") -> RateModel:
    """Pooled counting estimator (Dayhoff / Kosiol & Goldman 2005 style, matrix-log variant).

    F = sum_pairs (N + N^T)/2 + pseudocount (symmetric joint counts); pi = row margins; P = F / pi;
    Q = logm(P) via the symmetric eigendecomposition (eigenvalues clipped to (0,1]); S_ij = Q_ij / pi_j.
    Biased if pooled pairs differ much in divergence -- use close pairs, then refine with estimate_gtr().
    """
    c = np.asarray(counts, float)
    if c.ndim == 3:
        c = c.sum(0)
    F = 0.5 * (c + c.T) + pseudocount
    pi = F.sum(1) / F.sum()
    J = F / F.sum()
    sq = np.sqrt(pi)
    B = J / sq[:, None] / sq[None, :]                    # Pi^{-1/2} J Pi^{-1/2} = Pi^{1/2} P Pi^{-1/2}
    lam, U = np.linalg.eigh(0.5 * (B + B.T))
    lam = np.clip(lam, 1e-6, 1.0)
    L = (U * np.log(lam)[None]) @ U.T                     # Pi^{1/2} logm(P) Pi^{-1/2}
    Q = L * sq[None, :] / sq[:, None]
    S = Q / pi[None, :]
    S = 0.5 * (S + S.T)
    np.fill_diagonal(S, 0)
    S = np.clip(S, 1e-4 * np.abs(S).max(), None)
    np.fill_diagonal(S, 0)
    m = RateModel(S, pi, alphabet, name, "counting_gtr", dict(pseudocount=pseudocount))
    return RateModel(m.normalized_S(), pi, alphabet, name, "counting_gtr", dict(pseudocount=pseudocount))


def estimate_gtr(counts, alphabet: str, alpha: float | None = 1.0, ncat: int = 4, pseudocount: float = 0.5,
                 n_rounds: int = 6, lbfgs_iter: int = 60, device: str | None = None, init: RateModel | None = None,
                 name: str = "GTR-ML", verbose: bool = False) -> RateModel:
    """ML estimate of a GTR-K model from count matrices of (preferably close) pairs.

    pi is fixed at the pooled empirical frequencies (+pseudocount). Alternates (i) per-pair ML distances under
    the current model (ml_distance) and (ii) L-BFGS on log-exchangeabilities given distances, maximising
    sum_pairs sum_ij N_ij log P_ij(d_p) with fixed discrete Gamma(alpha). Starts from counting_gtr().
    """
    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"
    c = np.asarray(counts, float)
    K = len(alphabet)
    m0 = init if init is not None else counting_gtr(c, alphabet, pseudocount)
    pooled = 0.5 * (c.sum(0) + c.sum(0).T) + pseudocount
    pi = pooled.sum(1) / pooled.sum()
    iu = np.triu_indices(K, 1)
    theta = torch.tensor(np.log(np.maximum(m0.S[iu], 1e-8)), dtype=torch.float64, device=device,
                         requires_grad=True)
    pi_t = torch.tensor(pi, dtype=torch.float64, device=device)
    C = torch.tensor(c, dtype=torch.float64, device=device)
    r, w = gamma_rates(alpha, ncat, 0.0)
    r_t = torch.tensor(r, dtype=torch.float64, device=device)
    w_t = torch.tensor(w, dtype=torch.float64, device=device)
    iu_t = (torch.tensor(iu[0], device=device), torch.tensor(iu[1], device=device))
    sq = torch.sqrt(pi_t)

    def build_S(th):
        S = torch.zeros(K, K, dtype=torch.float64, device=device)
        S = S.index_put(iu_t, torch.exp(th))
        return S + S.T

    def current_model():
        S = build_S(theta.detach()).cpu().numpy()
        m = RateModel(S, pi, alphabet, name, "estimate_gtr")
        return RateModel(m.normalized_S(), pi, alphabet, name, "estimate_gtr",
                         dict(alpha=alpha, ncat=ncat, pseudocount=pseudocount, n_pairs=int(c.shape[0])))

    for rnd in range(n_rounds):
        res = ml_distance(c, current_model(), alpha=alpha, ncat=ncat, device=device)
        ok = np.isfinite(res.d) & ~res.saturated
        dd = torch.tensor(np.where(ok, res.d, 0.0), dtype=torch.float64, device=device)
        Cw = C[torch.tensor(ok, device=device)]
        dd = dd[torch.tensor(ok, device=device)]
        opt = torch.optim.LBFGS([theta], max_iter=lbfgs_iter, line_search_fn="strong_wolfe")

        def closure():
            opt.zero_grad()
            S = build_S(theta)
            Q = S * pi_t[None, :]
            Q = Q - torch.diag(Q.sum(1))
            Q = Q / (-(pi_t * torch.diagonal(Q)).sum())
            A = (sq[:, None] * Q) / sq[None, :]
            lam, U = torch.linalg.eigh(0.5 * (A + A.T))
            V = U / sq[:, None]
            Wm = U.T * sq[None, :]
            e = (w_t[None, :, None] * torch.exp(dd[:, None, None] * lam[None, None, :] * r_t[None, :, None])).sum(1)
            P = torch.einsum("im,nm,mj->nij", V, e, Wm).clamp_min(1e-300)
            loss = -(Cw * torch.log(P)).sum() / Cw.sum()
            loss.backward()
            return loss

        loss = opt.step(closure)
        if verbose:
            print(f"round {rnd}: mean -loglik/site {float(loss):.6f}, pairs used {int(ok.sum())}")
    return current_model()


# ----------------------------------------------------------------------------------------------------------
# Structure-level
# ----------------------------------------------------------------------------------------------------------
def one_minus(x):
    """1 - x (TM-score / LDDT -> distance)."""
    return 1.0 - np.asarray(x, dtype=float)


def neg_log(x, floor: float = 1e-6):
    """-ln(x) (FoldTree-style saturation-correctable transform); x clipped at `floor`."""
    return -np.log(np.clip(np.asarray(x, dtype=float), floor, None))


def kabsch_rmsd(X, Y, mask=None, return_transform: bool = False):
    """Batched Kabsch RMSD between paired coordinate sets.

    X, Y : (B, L, 3) or (L, 3) arrays/tensors of *aligned* (equivalenced) coordinates.
    mask : optional (B, L) bool / weights; masked-out positions (e.g. gaps, low pLDDT) are ignored.
    Returns rmsd (B,) [and R (B,3,3), t (B,3) such that X @ R^T + t ~ Y]. Proper rotations only (det R = +1).
    """
    is_np = not torch.is_tensor(X)
    X = torch.as_tensor(X, dtype=torch.float64)
    Y = torch.as_tensor(Y, dtype=torch.float64, device=X.device)
    single = X.ndim == 2
    if single:
        X, Y = X[None], Y[None]
        if mask is not None:
            mask = torch.as_tensor(mask)[None]
    if mask is None:
        w = torch.ones(X.shape[:2], dtype=X.dtype, device=X.device)
    else:
        w = torch.as_tensor(mask, device=X.device).to(X.dtype)
    W = w.sum(1).clamp_min(1e-12)
    cx = (w[..., None] * X).sum(1) / W[:, None]
    cy = (w[..., None] * Y).sum(1) / W[:, None]
    Xc, Yc = X - cx[:, None], Y - cy[:, None]
    H = torch.einsum("bl,bli,blj->bij", w, Xc, Yc)
    U, S, Vh = torch.linalg.svd(H)
    dsign = torch.sign(torch.linalg.det(Vh.transpose(1, 2) @ U.transpose(1, 2)))
    dsign = torch.where(dsign == 0, torch.ones_like(dsign), dsign)
    D = torch.diag_embed(torch.stack([torch.ones_like(dsign), torch.ones_like(dsign), dsign], -1))
    R = Vh.transpose(1, 2) @ D @ U.transpose(1, 2)
    Xr = Xc @ R.transpose(1, 2)
    msd = (w * ((Xr - Yc) ** 2).sum(-1)).sum(1) / W
    rmsd = torch.sqrt(msd.clamp_min(0))
    rmsd = torch.where(w.sum(1) > 0, rmsd, torch.full_like(rmsd, float("nan")))
    t = cy - torch.einsum("bij,bj->bi", R, cx)
    if single:
        rmsd, R, t = rmsd[0], R[0], t[0]
    if is_np:
        rmsd, R, t = rmsd.cpu().numpy(), R.cpu().numpy(), t.cpu().numpy()
    return (rmsd, R, t) if return_transform else rmsd
