#!/usr/bin/env python3
"""Validate hog_pair_labels.py on Ancestra's toy simulation (known truth).
Usage: python scripts/ancestra/validate_labels_toy.py PRED_LABELS.tsv TRUTH_LABELS.tsv PRED_ORTHOXML"""
import sys
import pandas as pd
from ancestra.orthoxml import hogs_to_pairs
pred, truth, ox = sys.argv[1:4]
p = pd.read_csv(pred, sep="\t"); t = pd.read_csv(truth, sep="\t")
o, pa = hogs_to_pairs(ox, with_paralogs=True)
po = set(map(tuple, p[p.relation == "ortholog"][["gene_a", "gene_b"]].values))
pp = set(map(tuple, p[p.relation.str.endswith("paralog")][["gene_a", "gene_b"]].values))
print("orthologs == ancestra.hogs_to_pairs:", po == o, "| paralogs == hogs_to_pairs:", pp == pa,
      "| ortholog TaxRange == species MRCA for all:", bool(p[p.relation == "ortholog"].lca_ok.all()))
m = p.merge(t, on=["gene_a", "gene_b"], how="outer", suffixes=("_pred", "_truth"))
m = m.fillna({"relation_pred": "none", "relation_truth": "none"})
print(pd.crosstab(m.relation_truth, m.relation_pred, margins=True))
both = m[(m.relation_pred == "ortholog") & (m.relation_truth == "ortholog")]
print("ortholog pairs in both: %d, speciation level agrees: %.4f" % (len(both), (both.event_level_pred == both.event_level_truth).mean()))
bp = m[m.relation_pred.str.endswith("paralog") & m.relation_truth.str.endswith("paralog")]
print("paralog pairs in both: %d, duplication level agrees: %.4f" % (len(bp), (bp.event_level_pred == bp.event_level_truth).mean()))
# precision/recall of paralog call among pairs that are homologous in truth
for lab in ("ortholog", "paralog"):
    pr = m.relation_pred.str.endswith(lab) if lab == "paralog" else m.relation_pred == lab
    tr = m.relation_truth.str.endswith(lab) if lab == "paralog" else m.relation_truth == lab
    tp = (pr & tr).sum()
    print(f"{lab}: P={tp/pr.sum():.4f} R={tp/tr.sum():.4f}")
