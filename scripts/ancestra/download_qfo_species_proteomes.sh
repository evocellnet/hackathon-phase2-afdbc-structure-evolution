#!/usr/bin/env bash
# Download current UniProt reference proteomes (canonical, one protein per gene) for the
# 78 leaves of Ancestra's QfO 2022 guide tree. Output: data/raw/ancestra_test/proteomes/<UPID>_<taxid>.fasta.gz
# Usage: bash scripts/ancestra/download_qfo_species_proteomes.sh [OUTDIR]
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
OUT="${1:-$ROOT/data/raw/ancestra_test/proteomes}"
mkdir -p "$OUT"
TREE="$ROOT/tools/ancestra-mvp-ai/inputs/qfo2022/guide_tree.nwk"
BASE=https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes
grep -oE 'UP[0-9]+_[0-9]+' "$TREE" | sort -u | while read leaf; do
  up=${leaf%_*}
  [ -s "$OUT/$leaf.fasta.gz" ] && continue
  ok=0
  for k in Eukaryota Bacteria Archaea; do
    if wget -q -O "$OUT/$leaf.fasta.gz.tmp" "$BASE/$k/$up/$leaf.fasta.gz"; then
      mv "$OUT/$leaf.fasta.gz.tmp" "$OUT/$leaf.fasta.gz"; ok=1; break
    fi
  done
  rm -f "$OUT/$leaf.fasta.gz.tmp"
  [ $ok = 1 ] && echo "ok $leaf" || echo "MISSING $leaf"
done
wget -q -O "$OUT/../RELEASE.metalink" "$BASE/RELEASE.metalink" || true
