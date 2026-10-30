#!/usr/bin/env bash
# Download UniProt reference proteomes (canonical, one protein per gene; current release) for every upid in
# one or more panel TSVs (columns upid, taxid, superkingdom). Resumable; 8 parallel downloads.
# Output: data/raw/proteomes/<UPID>_<taxid>.fasta.gz ; md5 list data/raw/proteomes/md5sums.txt
# Usage: bash scripts/panel/download_proteomes.sh data/processed/panel/species_panel.tsv [more.tsv ...]
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
OUT="$ROOT/data/raw/proteomes"; mkdir -p "$OUT"
BASE=https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes
fetch() {  # upid taxid kingdom
  local up=$1 tx=$2 k=$3 f="$OUT/${1}_${2}.fasta.gz"
  [ -s "$f" ] && gzip -t "$f" 2>/dev/null && return 0
  for i in 1 2 3; do
    if wget -q -O "$f.tmp" "$BASE/$k/$up/${up}_${tx}.fasta.gz" && gzip -t "$f.tmp"; then mv "$f.tmp" "$f"; echo "ok $up"; return 0; fi
    sleep 3
  done
  rm -f "$f.tmp"; echo "MISSING $up $tx $k"
}
export -f fetch; export OUT BASE
for tsv in "$@"; do
  awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}{print $h["upid"], $h["taxid"], $h["superkingdom"]}' "$tsv"
done | sort -u | xargs -P 8 -n 3 bash -c 'fetch "$@"' _
wget -q -O "$OUT/RELEASE.metalink" "$BASE/RELEASE.metalink" || true
(cd "$OUT" && md5sum *.fasta.gz > md5sums.txt)
