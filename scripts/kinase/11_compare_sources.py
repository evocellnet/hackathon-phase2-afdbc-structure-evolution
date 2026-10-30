"""X7: alignment sensitivity summary (same fixed topology, three column sets).
Writes results/kinase/stems_by_source.tsv (group stems x source) and results/kinase/class_rho_by_source.tsv
(pooled and median rho per branch class x source; supported branches with b_AA >= 0.02), plus
results/kinase/sensitivity_summary.json (Spearman of per-branch rho between sources; flagged-branch overlap)."""
import json
import numpy as np, pandas as pd
R = "/scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution"
OUT = f"{R}/results/kinase"
SRC = ["linsi", "famsa", "foldmason"]
st = pd.concat([pd.read_csv(f"{OUT}/group_stems_{s}.tsv", sep="\t").assign(source=s) for s in SRC])
w = st.pivot_table(index="stem", columns="source", values=["rho", "rho_pct_vs_within_euk_ePK", "b_aa"], aggfunc="first")
w.columns = [f"{a}_{b}" for a, b in w.columns]
base = st[st.source == "linsi"].set_index("stem")[["n_members", "clade_size", "purity", "ufboot"]]
w = base.join(w).round(3)
w.to_csv(f"{OUT}/stems_by_source.tsv", sep="\t")
cl = pd.concat([pd.read_csv(f"{OUT}/class_rho_{s}.tsv", sep="\t").assign(source=s) for s in SRC])
cl = cl[cl.subset == "supported_baa>=0.02"].pivot_table(index="bclass", columns="source",
                                                           values=["n", "median_rho", "pooled_rho"], aggfunc="first")
cl.columns = [f"{a}_{b}" for a, b in cl.columns]
cl.round(3).to_csv(f"{OUT}/class_rho_by_source.tsv", sep="\t")
B = {s: pd.read_csv(f"{OUT}/branches_{s}.tsv", sep="\t") for s in SRC}
m = B["linsi"].supported & ~B["linsi"].terminal & (B["linsi"].b_aa >= 0.02)
summ = {}
for s in ("famsa", "foldmason"):
    summ[f"spearman_rho_internal_linsi_vs_{s}"] = float(B["linsi"].rho[m].corr(B[s].rho[m], method="spearman"))
    summ[f"spearman_baa_internal_linsi_vs_{s}"] = float(B["linsi"].b_aa[m].corr(B[s].b_aa[m], method="spearman"))
    summ[f"spearman_rho_terminal_linsi_vs_{s}"] = float(B["linsi"].rho[B["linsi"].terminal].corr(
        B[s].rho[B["linsi"].terminal], method="spearman"))
    f0 = set(B["linsi"].index[B["linsi"].flag.fillna("").str.contains("rho_high")])
    f1 = set(B[s].index[B[s].flag.fillna("").str.contains("rho_high")])
    summ[f"rho_high_flags_linsi_{s}_overlap"] = [len(f0), len(f1), len(f0 & f1)]
f_all = set.intersection(*[set(B[s].index[B[s].flag.fillna("").str.contains("rho_high")]) for s in SRC])
summ["rho_high_in_all_three"] = len(f_all)
keep = B["linsi"].loc[sorted(f_all), ["clade_size", "ufboot", "b_aa", "b_3di", "rho", "bclass", "cats", "egroups"]]
for s in ("famsa", "foldmason"):
    keep[f"rho_{s}"] = B[s].loc[sorted(f_all), "rho"]
keep.round(3).to_csv(f"{OUT}/flagged_rho_high_all_sources.tsv", sep="\t")
json.dump(summ, open(f"{OUT}/sensitivity_summary.json", "w"), indent=1)
pd.set_option("display.width", 250)
print(w.to_string()); print(cl.round(3).to_string()); print(json.dumps(summ, indent=1)); print(keep.round(3).to_string())
