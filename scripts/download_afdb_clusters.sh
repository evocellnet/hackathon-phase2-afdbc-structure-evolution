#!/usr/bin/env bash
# Download AFDB Clusters (Barrio-Hernandez et al. 2023) metadata tables, v3 and v6 releases.
# Source: https://afdb-cluster.steineggerlab.workers.dev/ -> https://steineggerlab.s3.amazonaws.com/afdb-cluster/
# License: CC-BY-4.0
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$ROOT/data/raw/afdb_clusters"
BASE="https://steineggerlab.s3.amazonaws.com/afdb-cluster"
mkdir -p "$OUT/v3" "$OUT/v6"
curl -s "https://steineggerlab.s3.amazonaws.com/?prefix=afdb-cluster/" > "$OUT/s3_listing.xml"
FILES=(
v3/1-AFDBClusters-entryId_repId_taxId.tsv.gz
v3/2-repId_isDark_nMem_repLen_avgLen_repPlddt_avgPlddt_LCAtaxId.tsv.gz
v3/3-sapId_sapGO_repId_cluFlag_LCAtaxId.tsv.gz
v3/4-domain-clustering.zip
v3/5-allmembers-repId-entryId-cluFlag-taxId.tsv.gz
v3/6-all-vs-all-similarity-queryId_targetId_eValue-Dec16_2024.tsv.gz
v3/7-AFDB50-repId_memId.tsv.gz
v3/8-cluster_alignment_tmscore_lddt.tar.gz
v3/commands.gz
v6/1-AFDBClusters-repId_entryId_cluFlag_taxId.tsv.gz
v6/2-repId_isDark_nMem_repLen_avgLen_repPlddt_avgPlddt_LCAtaxId.tsv.gz
v6/5-allmembers-repId-entryId-cluFlag-taxId.tsv.gz
v6/6-all-vs-all-similarity-queryId_targetId_eValue.tsv.gz
v6/7-AFDB50-repId_memId.tsv.gz
)
dl() { f="$1"; curl -sS -L --retry 5 -C - -o "$OUT/$f" "$BASE/$f" && echo "done $f"; }
export -f dl; export OUT BASE
printf '%s\n' "${FILES[@]}" | xargs -P 6 -I{} bash -c 'dl {}'
cd "$OUT" && md5sum v3/* v6/* > md5sums.txt
echo ALL_DONE
