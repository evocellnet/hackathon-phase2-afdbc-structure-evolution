"""X6 node 2, focal question (ii): lineage-specific structural acceleration within eukaryotes.

Uses results/casestudy/branch_ratios_{sf}_ml.tsv (17_ratio_analysis.py). Eukaryote-only comparable branches
(supported, 0.02 <= b_AA <= 2) whose clade lies within one supergroup (terminal branches and lineage-pure clades) are
grouped by supergroup and by selected genera (divergent protists). Reports median rho, the share flagged high/low
(comparable-length window) and robust flags; lists every robust within-eukaryote flag.
Outputs results/casestudy/lineage_rho.tsv, results/casestudy/flagged_branches.tsv
"""
import json
import numpy as np, pandas as pd

R = "/scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution"
OUT = f"{R}/results/casestudy"
GENERA = ["Giardia", "Trichomonas", "Spironucleus", "Plasmodium", "Toxoplasma", "Cryptosporidium", "Theileria", "Babesia",
          "Trypanosoma", "Leishmania", "Naegleria", "Tetrahymena", "Paramecium", "Arabidopsis", "Saccharomyces",
          "Schizosaccharomyces", "Entamoeba", "Dictyostelium", "Homo"]
sel = pd.read_csv(f"{OUT}/sample_members.tsv", sep="\t").set_index("seqid")
rows, flags = [], []
for sf in ("actin", "tubulin"):
    B = pd.read_csv(f"{OUT}/branch_ratios_{sf}_ml.tsv", sep="\t")
    s = sel[sel.superfamily == sf]
    leaves = sorted(s.index)
    def mem(h):
        m = int(h, 16)
        return [l for i, l in enumerate(leaves) if (m >> i) & 1]
    B["members"] = [mem(h) for h in B.clade_hex]
    B["euk_only"] = [all(s.at[l, "domain"] == "Eukaryota" and s.at[l, "subfamily"] != "euk_FtsZ" for l in ms) for ms in B.members]
    B["group1"] = [s.loc[ms, "group"].iloc[0] if s.loc[ms, "group"].nunique() == 1 else "" for ms in B.members]
    B["genus1"] = [s.loc[ms, "name"].str.split().str[0].iloc[0] if s.loc[ms, "name"].str.split().str[0].nunique() == 1 else ""
                   for ms in B.members]
    comp = B.supported & B.iq_linsi_b_aa.between(0.02, 2) & B.euk_only
    E = B[comp].copy()
    E["lr"] = np.log(E.iq_linsi_rho)
    allmed = float(np.exp(E.lr.median()))
    for key, col in (("supergroup", "group1"), ("genus", "genus1")):
        for g, d in E[E[col] != ""].groupby(col):
            if key == "genus" and g not in GENERA:
                continue
            rows.append(dict(sf=sf, level=key, lineage=g, n_branches=len(d), n_terminal=int(d.terminal.sum()),
                             median_rho=float(np.exp(d.lr.median())), q25=float(np.exp(d.lr.quantile(.25))),
                             q75=float(np.exp(d.lr.quantile(.75))), median_rho_all_euk=allmed,
                             median_b_aa=float(d.iq_linsi_b_aa.median()),
                             n_flag_high=int((d.iq_linsi_flag == "high").sum()), n_flag_low=int((d.iq_linsi_flag == "low").sum()),
                             n_robust_high=int((d.robust & (d.iq_linsi_flag == "high")).sum()),
                             n_robust_low=int((d.robust & (d.iq_linsi_flag == "low")).sum())))
    # AF-model quality check for terminal branches (3Di tip lengths include prediction noise)
    from scipy.stats import spearmanr
    Tm = B[B.terminal & B.supported & B.iq_linsi_b_aa.between(0.02, 2)].copy()
    Tm["tip_plddt"] = [s.loc[ms, "mean_plddt"].iloc[0] if len(ms) == 1 else s.loc[[l for l in leaves if l not in ms], "mean_plddt"].iloc[0]
                       for ms in Tm.members]
    rs = spearmanr(np.log(Tm.iq_linsi_rho), Tm.tip_plddt)
    rows.append(dict(sf=sf, level="QC", lineage="spearman(log rho_terminal, tip pLDDT)", n_branches=len(Tm),
                     spearman_r=float(rs[0]), spearman_p=float(rs[1])))
    B["tip_plddt"] = [s.loc[ms, "mean_plddt"].min() for ms in B.members]
    F = B[B.robust].copy()
    F["sf"] = sf
    F["where"] = np.where(F.euk_only, "within_eukaryote", "other")
    F["taxa_short"] = [", ".join(sorted({s.at[l, "name"].split(" (")[0] for l in ms}))[:150] for ms in F.members]
    flags.append(F)
L = pd.DataFrame(rows)
L.to_csv(f"{OUT}/lineage_rho.tsv", sep="\t", index=False)
F = pd.concat(flags)
keep = ["sf", "branch_id", "where", "clade_name", "n_clade", "terminal", "ufboot", "subfamilies", "groups", "taxa_short",
        "iq_linsi_b_aa", "iq_linsi_b_3di", "iq_linsi_rho", "iq_linsi_rho_lo", "iq_linsi_rho_hi", "iq_linsi_win_lo",
        "iq_linsi_win_hi", "iq_linsi_flag", "iq_linsi_p_beyond", "n_aln_agree", "iq_famsa_rho", "iq_fmg_rho",
        "iq_foldmason_rho", "tip_plddt", "eng_tdiaa_rho", "eng_tdiaa_logrho_se", "eng_tdiaa_rho_lddt", "b_aa_LG_C20"]
F = F[[c for c in keep if c in F]]
F.to_csv(f"{OUT}/flagged_branches.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 70)
print(L.round(3).to_string())
print(F.drop(columns=["subfamilies", "groups"]).round(3).to_string())
