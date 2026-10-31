"""V11 input sanity: split membership, primary pair completeness, E6 pair/triplet counts, branch-table coverage."""
import numpy as np
from v11_common import *

res = {}
split = pl.read_csv(ROOT / "data/interim/split.tsv", separator="\t")
conf = set(split.filter(pl.col("split") == "confirm")["hog"].to_list())
tune = set(split.filter(pl.col("split") == "tune")["hog"].to_list())
h = hogs()
hc = h.filter(pl.col("split") == "confirm")
res["split_tsv_confirm"] = len(conf); res["hogs_parquet_confirm"] = hc.height
res["hogs_parquet_confirm_eq_split_tsv"] = set(hc["hog"].to_list()) == conf
expected = int(sum(len(s) * (len(s) - 1) // 2 for s in hc["selected"].to_list()))
for aln in ("tdiaa", "seqonly", "tmalign"):
    df = load_pairs(aln, SPLIT, columns=["pair_id", "pair_set", "hog", "split", "status", "species_a", "species_b"])
    pr = df.filter(pl.col("pair_set") == "primary")
    r = {"rows": df.height, "primary": pr.height, "primary_expected_sum_C(n,2)": expected,
         "primary_dup_ids": pr.height - pr["pair_id"].n_unique(),
         "hogs_not_confirm_in_split_tsv": int((~pr["hog"].is_in(list(conf))).sum()),
         "hogs_in_tune": int(pr["hog"].is_in(list(tune)).sum()),
         "confirm_hogs_without_pairs": len(conf - set(pr["hog"].unique().to_list())),
         "split_values": df["split"].unique().to_list(),
         "e6": df.filter(pl.col("pair_set").str.starts_with("e6_")).group_by("pair_set").len().to_dicts()}
    # per-HOG completeness: pairs per hog == C(n,2)
    n_sel = dict(zip(hc["hog"].to_list(), [len(s) for s in hc["selected"].to_list()]))
    g = pr.group_by("hog").len()
    r["hogs_pair_count_mismatch"] = int(sum(1 for hh, n in g.iter_rows() if n != n_sel[hh] * (n_sel[hh] - 1) // 2))
    res[aln] = r
T = pl.read_csv(ROOT / "data/interim/pairlists/e6_triplets.tsv", separator="\t", infer_schema_length=0)
res["e6_triplets_total"] = T.height; res["e6_triplets_by_split"] = T.group_by("split").len().to_dicts()
res["e6_triplet_hogs_not_confirm"] = int((~T.filter(pl.col("split") == "confirm")["hog"].is_in(list(conf))).sum())
log = pl.read_csv(ROOT / "logs/filters/hogs_e6.tsv", separator="\t")
res["e6_filter_log_triplets_all_splits"] = int(log.filter(pl.col("step") == "triplets")["n"][0])
b = load_branches(PRIMARY, SPLIT)
res["branch_table_hogs"] = b["hog"].n_unique(); res["branch_table_hogs_all_confirm"] = set(b["hog"].unique().to_list()) <= conf
rm = json.loads((ROOT / "results/endpoints/run_meta.json").read_text())
res["run_meta_n_pairs"] = rm["n_pairs"]; res["run_meta_n_hogs"] = rm["n_hogs"]; res["run_meta_split"] = rm["split"]
print(json.dumps(res, indent=1, default=str))
wjson(res, OUT / "V11_inputs.json")
