"""X7 extra summaries (descriptive):
 (1) TK origin segment: the two adjacent supported branches at the base of the holozoan TK clade (MRCA clade of 34
     with one non-TK intruder, and the pure 32-TK clade) pooled: rho_seg = sum b_3Di / sum b_AA / scale, compared with
     pooled rho of all 2-branch parent-child paths of supported internal within-euk-ePK branches (percentile).
 (2) Question (ii): pooled and median rho, prokaryotic vs eukaryotic ePK internal (supported, b_AA >= 0.02) and
     terminal branches, with a branch bootstrap (B = 2000, seed 20261030) 95 % interval of the ratio of pooled rho
     (prok/euk); and the same restricted to the overlapping b_AA range (b_AA-matched: 0.05-0.5).
Writes results/kinase/tk_segment_by_source.tsv, results/kinase/prok_vs_euk_by_source.tsv"""
import numpy as np, pandas as pd
from ete3 import Tree
R = "/scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution"
OUT = f"{R}/results/kinase"
rng = np.random.default_rng(20261030)
rows, rows2 = [], []
for src in ("linsi", "famsa", "foldmason"):
    b = pd.read_csv(f"{OUT}/branches_{src}.tsv", sep="\t")
    st = pd.read_csv(f"{OUT}/group_stems_{src}.tsv", sep="\t").set_index("stem")
    scale = b.b_3di.sum() / b.b_aa.sum()
    seg = st.loc[["TK holozoan MRCA", "TK holozoan (largest pure clade)"]]
    r_seg = seg.b_3di.sum() / seg.b_aa.sum() / scale
    bg = b[(b.bclass == "within_euk_ePK") & b.supported & ~b.terminal & (b.b_aa >= 0.02)]
    # background: single branches and branch pairs (random pairs of background branches, approximate)
    i = rng.integers(0, len(bg), (5000, 2))
    pair = (bg.b_3di.values[i].sum(1) / bg.b_aa.values[i].sum(1)) / scale
    rows.append(dict(source=src, b_aa_seg=seg.b_aa.sum(), b_3di_seg=seg.b_3di.sum(), rho_seg=r_seg,
                     pct_vs_single_bg=(bg.rho < r_seg).mean() * 100, pct_vs_pair_bg=(pair < r_seg).mean() * 100,
                     bg_median=bg.rho.median(), bg_q95=bg.rho.quantile(.95)))
    for kind, cp, ce in (("internal", "within_prok_ePK", "within_euk_ePK"),
                         ("terminal", "terminal_prok_ePK", "terminal_euk_ePK")):
        for rng_lab, lo, hi in (("all b_AA>=0.02", 0.02, 99), ("b_AA 0.05-0.5", 0.05, 0.5)):
            s = b[b.supported & (b.b_aa >= lo) & (b.b_aa <= hi)]
            P, E = s[s.bclass == cp], s[s.bclass == ce]
            f = lambda x: x.b_3di.sum() / x.b_aa.sum() / scale
            bs = []
            for _ in range(2000):
                p_ = P.iloc[rng.integers(0, len(P), len(P))]; e_ = E.iloc[rng.integers(0, len(E), len(E))]
                bs.append(f(p_) / f(e_))
            rows2.append(dict(source=src, kind=kind, b_aa_range=rng_lab, n_prok=len(P), n_euk=len(E),
                              pooled_prok=f(P), pooled_euk=f(E), ratio_pooled=f(P) / f(E),
                              ci_lo=np.quantile(bs, .025), ci_hi=np.quantile(bs, .975),
                              median_prok=P.rho.median(), median_euk=E.rho.median(),
                              median_baa_prok=P.b_aa.median(), median_baa_euk=E.b_aa.median()))
a = pd.DataFrame(rows).round(3); c = pd.DataFrame(rows2).round(3)
a.to_csv(f"{OUT}/tk_segment_by_source.tsv", sep="\t", index=False)
c.to_csv(f"{OUT}/prok_vs_euk_by_source.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
print(a.to_string()); print(c.to_string())
