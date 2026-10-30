"""Shared helpers for the species-panel pipeline: NCBI taxonomy and a cached TimeTree API client.

TimeTree data may not be redistributed: the cache lives under data/interim/timetree/ (never committed).
"""
from __future__ import annotations

import concurrent.futures as cf
import csv
import io
import json
import os
import threading
import time
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
TAXDUMP = os.path.join(ROOT, "data/raw/ncbi_taxonomy/new_taxdump")
TT_CACHE_DIR = os.path.join(ROOT, "data/interim/timetree/api_cache")
TT_BASE = "https://timetree.org/api"
UA = {"User-Agent": "Mozilla/5.0 (phase2-afdbc-evolution research; ETH Zurich)"}


# ----------------------------------------------------------------------------- taxonomy
class Taxonomy:
    def __init__(self, taxdump: str = TAXDUMP):
        self.parent: dict[int, int] = {}
        self.rank: dict[int, str] = {}
        with open(os.path.join(taxdump, "nodes.dmp")) as fh:
            for line in fh:
                f = line.split("\t|\t", 3)
                t = int(f[0])
                self.parent[t] = int(f[1])
                self.rank[t] = f[2]
        self.name: dict[int, str] = {}
        with open(os.path.join(taxdump, "names.dmp")) as fh:
            for line in fh:
                f = line.rstrip("\t|\n").split("\t|\t")
                if f[3] == "scientific name":
                    self.name[int(f[0])] = f[1]
        self.merged: dict[int, int] = {}
        with open(os.path.join(taxdump, "merged.dmp")) as fh:
            for line in fh:
                f = line.rstrip("\t|\n").split("\t|\t")
                self.merged[int(f[0])] = int(f[1])

    def resolve(self, t: int) -> int:
        return self.merged.get(t, t)

    def lineage(self, t: int) -> list[int]:
        """root-first list of taxids from root (1) to t."""
        t = self.resolve(t)
        out = [t]
        while t != 1 and t in self.parent:
            t = self.parent[t]
            out.append(t)
        return out[::-1]

    def species_of(self, t: int) -> int:
        """taxid of rank 'species' at or above t (t itself if none found)."""
        for x in reversed(self.lineage(t)):
            if self.rank.get(x) == "species":
                return x
        return self.resolve(t)

    def rank_ancestor(self, t: int, rank: str) -> int | None:
        for x in self.lineage(t):
            if self.rank.get(x) == rank:
                return x
        return None

    def mrca(self, a: int, b: int) -> int:
        la, lb = self.lineage(a), self.lineage(b)
        m = 1
        for x, y in zip(la, lb):
            if x != y:
                break
            m = x
        return m


# ----------------------------------------------------------------------------- TimeTree
class TimeTree:
    """Thread-safe, disk-cached client for the (undocumented) TimeTree 6 REST API.

    pairwise(a, b) -> dict(all_total, precomputed_age, ci_low, ci_high, adjusted_age, name_a, name_b, version)
                      or {} when TimeTree returns no row.
    taxon(t)       -> dict (JSON of /api/taxon) or {} if not found.
    """

    def __init__(self, cache_dir: str = TT_CACHE_DIR, max_workers: int = 8):
        os.makedirs(cache_dir, exist_ok=True)
        self.pw_path = os.path.join(cache_dir, "pairwise.jsonl")
        self.tx_path = os.path.join(cache_dir, "taxon.jsonl")
        self.max_workers = max_workers
        self.lock = threading.Lock()
        self.pw: dict[tuple[int, int], dict] = {}
        self.tx: dict[int, dict] = {}
        for p, d, key in ((self.pw_path, self.pw, lambda r: tuple(r["key"])), (self.tx_path, self.tx, lambda r: r["key"])):
            if os.path.exists(p):
                with open(p) as fh:
                    for line in fh:
                        try:
                            r = json.loads(line)
                        except json.JSONDecodeError:
                            continue
                        d[key(r)] = r["val"]

    @staticmethod
    def _get(url: str, tries: int = 4) -> str | None:
        for i in range(tries):
            try:
                req = urllib.request.Request(url, headers=UA)
                with urllib.request.urlopen(req, timeout=60) as r:
                    return r.read().decode()
            except urllib.error.HTTPError as e:
                # TimeTree answers 500 (deterministically) for pairs it cannot date (e.g. human vs E. coli);
                # treat 400/404/500 as "no answer" after one retry.
                if e.code in (400, 404) or (e.code == 500 and i >= 1):
                    return e.read().decode(errors="replace")
                time.sleep(1)
            except Exception:
                time.sleep(2 * (i + 1))
        return None

    def _pairwise_fetch(self, a: int, b: int) -> dict | None:
        txt = self._get(f"{TT_BASE}/pairwise/{a}/{b}")
        if txt is None:
            return None  # network failure: do not cache
        rows = list(csv.DictReader(io.StringIO(txt)))
        if not rows or "precomputed_age" not in rows[0]:
            return {}

        def fl(x):
            try:
                return float(x)
            except (TypeError, ValueError):
                return None
        r = rows[0]
        return {"name_a": r.get("scientific_name_a"), "name_b": r.get("scientific_name_b"),
                "all_total": fl(r.get("all_total")), "precomputed_age": fl(r.get("precomputed_age")),
                "ci_low": fl(r.get("precomputed_ci_low")), "ci_high": fl(r.get("precomputed_ci_high")),
                "adjusted_age": fl(r.get("adjusted_age"))}

    def _taxon_fetch(self, t: int) -> dict | None:
        txt = self._get(f"{TT_BASE}/taxon/{t}")
        if txt is None:
            return None
        try:
            j = json.loads(txt)
        except json.JSONDecodeError:
            return None
        if "taxon_id" not in j:
            return {}
        return {k: j.get(k) for k in ("taxon_id", "scientific_name", "taxonomic_rank", "ncbi_id", "isLeaf", "version")} | {
            "n_study_times": len(j.get("study_times") or [])}

    def _run(self, keys, fetch, store, path):
        todo = [k for k in dict.fromkeys(keys) if k not in store]
        if todo:
            with cf.ThreadPoolExecutor(self.max_workers) as ex, open(path, "a", buffering=1) as fh:
                futs = {ex.submit(fetch, *(k if isinstance(k, tuple) else (k,))): k for k in todo}
                for i, f in enumerate(cf.as_completed(futs)):
                    k = futs[f]
                    v = f.result()
                    if v is None:
                        continue
                    with self.lock:
                        store[k] = v
                        fh.write(json.dumps({"key": list(k) if isinstance(k, tuple) else k, "val": v}) + "\n")
                    if (i + 1) % 1000 == 0:
                        fh.flush()
                        print(f"  timetree {os.path.basename(path)}: {i+1}/{len(todo)}", flush=True)
        return {k: store.get(k) for k in keys}

    def pairwise_many(self, pairs):
        keys = [(min(a, b), max(a, b)) for a, b in pairs]
        return self._run(keys, self._pairwise_fetch, self.pw, self.pw_path)

    def taxon_many(self, taxids):
        return self._run(list(taxids), self._taxon_fetch, self.tx, self.tx_path)


def parse_busco(s: str | None) -> float | None:
    """'C:99.6%[S:97.4%,D:2.2%],F:0.1%,M:0.3%,n:2137' -> 99.6"""
    if not s or not s.startswith("C:"):
        return None
    try:
        return float(s[2:s.index("%")])
    except ValueError:
        return None
