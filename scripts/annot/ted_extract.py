#!/usr/bin/env python3
"""Step 2a: extract TED domain rows and per-chain segmentation rows for panel proteins.
Outputs data/interim/annotations/ted_domains_panel.parquet and ted_chains_panel.parquet.
TED domain summary columns (21, no header; verified on data and scout report):
 1 id (AF-<acc>-F1-model_v4_TEDnn) 2 md5(domain seq) 3 consensus level (high/medium) 4 chopping (e.g. 2-143_287-314)
 5 nres (domain residues) 6 n_segments 7 domain pLDDT 8-12 packing/SSE metrics 13 proteome 14 CATH label (C.A.T.H, or C.A.T at level T, '-' = none)
 15 CATH level (H/T/-) 16 assignment method (foldseek/foldclass/-) 17 norm_rg 18 globularity 19 '-' 20 species 21 lineage
Per-chain file: id, md5, nres, n_high, n_medium, n_low, chop_high, chop_medium, chop_low, proteome.
"""
import duckdb, sys, time
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]; RAW = ROOT/"data/raw/ted"
OUT = ROOT/"data/interim/annotations"; OUT.mkdir(parents=True, exist_ok=True)
which = sys.argv[1]
con = duckdb.connect(); con.execute("set threads=4; set memory_limit='30GB'")
con.execute(f"""create table panel as select acc from read_parquet('{ROOT}/data/processed/panel/protein_table.parquet')
 union select acc from read_parquet('{ROOT}/data/processed/panel/protein_table_prok.parquet')""")
t = time.time()
if which == "domains":
    cols = {f"c{i}": "VARCHAR" for i in range(1, 22)}
    con.execute(f"""copy (select regexp_extract(c1,'^AF-(.+)-F1-model_v4_TED',1) acc, c1 ted_id, c2 dom_md5, c3 consensus, c4 chopping,
      try_cast(c5 as integer) dom_nres, try_cast(c6 as integer) n_seg, try_cast(c7 as double) dom_plddt, c14 cath_label, c15 cath_level, c16 cath_method,
      try_cast(c17 as double) norm_rg, try_cast(c18 as double) globularity
      from read_csv('{RAW}/ted_365m.domain_summary.cath.globularity.taxid.tsv.gz',delim='\t',header=false,quote='',escape='',columns={cols})
      where regexp_extract(c1,'^AF-(.+)-F1-model_v4_TED',1) in (select acc from panel)) to '{OUT}/ted_domains_panel.parquet' (format parquet)""")
else:
    cols = {f"c{i}": "VARCHAR" for i in range(1, 11)}
    con.execute(f"""copy (select regexp_extract(c1,'^AF-(.+)-F1-model_v4$',1) acc, c2 chain_md5, try_cast(c3 as integer) nres, try_cast(c4 as integer) n_high,
      try_cast(c5 as integer) n_medium, try_cast(c6 as integer) n_low, c7 chop_high, c8 chop_medium, c9 chop_low
      from read_csv('{RAW}/ted_214m_per_chain_segmentation.tsv.gz',delim='\t',header=false,quote='',escape='',columns={cols})
      where regexp_extract(c1,'^AF-(.+)-F1-model_v4$',1) in (select acc from panel)) to '{OUT}/ted_chains_panel.parquet' (format parquet)""")
print(which, time.time()-t)
