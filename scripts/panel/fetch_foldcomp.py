#!/usr/bin/env python3
"""Fetch AFDB v4 structures from the remote FoldComp DB afdb_uniprot_v4 by HTTP range requests and build a
local FoldComp DB (data/interim/structures/<name> + .index/.lookup/.dbtype), resumable.

Input: data/interim/panel/proteins_2026.parquet (fc_key of every panel protein with an AFDB v4 model),
restricted to the leaves of the given panel TSV(s).
Method (data scout report §3): byte offsets from the local afdb_uniprot_v4.index; S3 returns 206 for Range
requests.  --procs processes x --threads threads; each process appends to its own shard (data + done-list),
so a rerun skips entries already fetched.  `--merge` concatenates shards into one FoldComp DB sorted by key.
"""
from __future__ import annotations

import argparse
import multiprocessing as mp
import os
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor

import polars as pl
import requests

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import ROOT  # noqa: E402

URL = "https://steineggerlab.s3.amazonaws.com/foldcomp/afdb_uniprot_v4"
INDEX = os.path.join(ROOT, "data/raw/foldcomp/afdb_uniprot_v4.index")
INTERIM = os.path.join(ROOT, "data/interim/panel")


def targets(panels):
    leaves = set()
    for p in panels:
        leaves |= set(pl.read_csv(p, separator="\t")["leaf"].to_list())
    prot = pl.read_parquet(os.path.join(INTERIM, "proteins_2026.parquet"), columns=["leaf", "fc_key", "entry"])
    prot = prot.filter(pl.col("leaf").is_in(list(leaves)) & pl.col("fc_key").is_not_null()).unique("fc_key")
    idx = (pl.scan_csv(INDEX, separator="\t", has_header=False, new_columns=["fc_key", "off", "len"],
                       schema_overrides={"fc_key": pl.Int64, "off": pl.Int64, "len": pl.Int64})
           .join(prot.lazy(), on="fc_key", how="inner").collect().sort("off"))
    return idx


def worker(args):
    wid, rows, shard_dir, threads = args
    data_p = os.path.join(shard_dir, f"shard{wid:02d}.data")
    done_p = os.path.join(shard_dir, f"shard{wid:02d}.done")
    done = set()
    if os.path.exists(done_p):
        # a done line is written only after its bytes were flushed to the data file
        with open(done_p) as fh:
            for line in fh:
                p = line.split("\t")
                if len(p) == 4:
                    done.add(int(p[0]))
        # truncate data to last recorded end (drop partial writes)
        end = 0
        with open(done_p) as fh:
            for line in fh:
                p = line.split("\t")
                if len(p) == 4:
                    end = max(end, int(p[1]) + int(p[2]))
        if os.path.exists(data_p):
            with open(data_p, "r+b") as fh:
                fh.truncate(end)
    todo = [r for r in rows if r[0] not in done]
    if not todo:
        return wid, 0, 0
    lock = threading.Lock()
    local = threading.local()
    fdata = open(data_p, "ab")
    fdone = open(done_p, "a")
    pos = fdata.tell()
    n_ok, n_fail = 0, 0

    def get(r):
        key, off, ln, entry = r
        s = getattr(local, "s", None)
        if s is None:
            s = local.s = requests.Session()
        for i in range(6):
            try:
                resp = s.get(URL, headers={"Range": f"bytes={off}-{off + ln - 1}"}, timeout=60)
                if resp.status_code == 206 and len(resp.content) == ln:
                    return r, resp.content
            except Exception:
                pass
            time.sleep(1 + 2 * i)
        return r, None
    t0 = time.time()
    with ThreadPoolExecutor(threads) as ex:
        for k, (r, b) in enumerate(ex.map(get, todo)):
            if b is None:
                n_fail += 1
                continue
            with lock:
                fdata.write(b)
                fdata.flush()
                fdone.write(f"{r[0]}\t{pos}\t{len(b)}\t{r[3]}\n")
                pos += len(b)
                n_ok += 1
            if (k + 1) % 20000 == 0:
                fdone.flush()
                print(f"[w{wid}] {k+1}/{len(todo)} {(k+1)/(time.time()-t0):.0f}/s", flush=True)
    fdata.close()
    fdone.close()
    return wid, n_ok, n_fail


def merge(shard_dir, out):
    os.makedirs(os.path.dirname(out), exist_ok=True)
    recs = []
    for f in sorted(os.listdir(shard_dir)):
        if f.endswith(".done"):
            wid = f.replace(".done", "")
            with open(os.path.join(shard_dir, f)) as fh:
                for line in fh:
                    p = line.rstrip("\n").split("\t")
                    if len(p) == 4:
                        recs.append((int(p[0]), wid, int(p[1]), int(p[2]), p[3]))
    recs.sort()
    seen = set()
    handles = {}
    off = 0
    with open(out + ".tmp", "wb") as fo, open(out + ".index.tmp", "w") as fi, open(out + ".lookup.tmp", "w") as fl:
        for key, wid, o, ln, entry in recs:
            if key in seen:
                continue
            seen.add(key)
            h = handles.get(wid) or handles.setdefault(wid, open(os.path.join(shard_dir, wid + ".data"), "rb"))
            h.seek(o)
            b = h.read(ln)
            fo.write(b)
            fi.write(f"{key}\t{off}\t{ln}\n")
            fl.write(f"{key}\t{entry}\t0\n")
            off += ln
    for h in handles.values():
        h.close()
    for ext in ("", ".index", ".lookup"):
        os.replace(out + ext + ".tmp", out + ext)
    with open(os.path.join(ROOT, "data/raw/foldcomp/afdb_uniprot_v4.dbtype"), "rb") as fh, open(out + ".dbtype", "wb") as fo:
        fo.write(fh.read())
    print(f"merged {len(seen)} entries into {out} ({off/1e9:.2f} GB)")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--panel", nargs="+", required=True)
    ap.add_argument("--name", default="panel_v4")
    ap.add_argument("--procs", type=int, default=8)
    ap.add_argument("--threads", type=int, default=96)
    ap.add_argument("--merge", action="store_true")
    a = ap.parse_args()
    sdir = os.path.join(ROOT, "data/interim/structures", a.name + "_shards")
    os.makedirs(sdir, exist_ok=True)
    out = os.path.join(ROOT, "data/interim/structures", a.name)
    if not a.merge:
        idx = targets(a.panel)
        rows = list(idx.select("fc_key", "off", "len", "entry").iter_rows())
        print(f"{len(rows)} entries to fetch ({idx['len'].sum()/1e9:.2f} GB)", flush=True)
        # deterministic sharding by key so that resumption reuses the same shard
        shards = [[] for _ in range(a.procs)]
        for r in rows:
            shards[r[0] % a.procs].append(r)
        t = time.time()
        with mp.Pool(a.procs) as pool:
            for wid, ok, fail in pool.imap_unordered(worker, [(i, s, sdir, a.threads) for i, s in enumerate(shards)]):
                print(f"[w{wid}] done: fetched {ok}, failed {fail}", flush=True)
        print(f"fetch finished in {time.time()-t:.0f}s", flush=True)
    merge(sdir, out)


if __name__ == "__main__":
    main()
