#!/usr/bin/env python3
"""Run Ancestra (tools/ancestra-mvp-ai, frozen freeze-v1 method) on a protein set and label all pairs.

Three input modes (exactly one):
  --proteomes DIR   one FASTA per species; file stem = species name. Stems like UP000005640_9606
                    (UniProt reference proteome) give the taxid automatically; otherwise pass --taxmap
                    with lines "<file stem>\t<taxid>".
  --fasta FILE      any (multi-species) FASTA, e.g. the members of one AFDB cluster / Pfam family.
                    Taxid per sequence from the UniProt header (OX=<taxid>) or from --taxmap
                    ("<protein id>\t<taxid>"). Split into one proteome per taxid (leaf name t<taxid>).
  --structures DIR  PDB/mmCIF files (e.g. AFDB models). Sequences are extracted with foldseek
                    (createdb + convert2fasta); the 3Di strings are written too (usable with
                    --family-rescue 3di, off by default since it lowered precision in the Ancestra
                    benchmarks). Needs --taxmap "<structure id>\t<taxid>"; AFDB names
                    AF-<ACC>-F1-model_v4 are reduced to <ACC>.

Then: NCBI-taxonomy guide tree (ancestra.speciestree.ncbi_tree) -> `ancestra run` -> pair labels
(scripts/ancestra/hog_pair_labels.py) -> timing.json (wall, peak RSS of the whole run).

Example:
  conda_envs/ancestra/bin/python scripts/ancestra/run_ancestra.py --fasta fam.fasta \
      --taxdump data/interim/ancestra/taxdump --out results/ancestra/runs/fam1 --threads 48 --preset fast
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import re
import resource
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
ENV_BIN = os.path.join(ROOT, "conda_envs", "ancestra", "bin")


def read_fasta(path):
    op = gzip.open if path.endswith(".gz") else open
    name, seq, hdr = None, [], None
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, hdr, "".join(seq)
                hdr = line[1:].strip()
                tok = hdr.split()[0]
                name = tok.split("|")[1] if tok.count("|") >= 2 else tok
                seq = []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, hdr, "".join(seq)


def load_taxmap(path):
    m = {}
    if path:
        with open(path) as fh:
            for line in fh:
                if line.strip() and not line.startswith("#"):
                    k, v = line.rstrip("\n").split("\t")[:2]
                    m[k] = int(v)
    return m


def split_by_taxid(records, taxmap, outdir, min_len=30):
    os.makedirs(outdir, exist_ok=True)
    by = {}
    n_skip = 0
    for name, hdr, seq in records:
        t = taxmap.get(name)
        if t is None:
            mm = re.search(r"OX=(\d+)", hdr or "")
            t = int(mm.group(1)) if mm else None
        if t is None or len(seq) < min_len:
            n_skip += 1
            continue
        by.setdefault(t, []).append((name, seq.replace("*", "")))
    leaf_taxid = {}
    for t, recs in by.items():
        leaf = f"t{t}"
        with open(os.path.join(outdir, leaf + ".fa"), "w") as fh:
            for n, s in recs:
                fh.write(f">{n}\n{s}\n")
        leaf_taxid[leaf] = t
    return leaf_taxid, n_skip


def structures_to_fasta(sdir, work, threads):
    os.makedirs(work, exist_ok=True)
    fs = shutil.which("foldseek") or os.path.join(ENV_BIN, "foldseek")
    db = os.path.join(work, "sdb")
    subprocess.run([fs, "createdb", sdir, db, "--threads", str(threads)], check=True,
                   stdout=subprocess.DEVNULL)
    aa, di = os.path.join(work, "aa.fasta"), os.path.join(work, "3di.fasta")
    subprocess.run([fs, "convert2fasta", db, aa], check=True, stdout=subprocess.DEVNULL)
    subprocess.run([fs, "lndb", db + "_h", db + "_ss_h"], check=True, stdout=subprocess.DEVNULL)
    subprocess.run([fs, "convert2fasta", db + "_ss", di], check=True, stdout=subprocess.DEVNULL)

    def clean(p):   # AF-P12345-F1-model_v4[.pdb|.cif][ chain] -> P12345
        out = p + ".clean"
        with open(p) as fi, open(out, "w") as fo:
            for line in fi:
                if line.startswith(">"):
                    tok = line[1:].split()[0]
                    tok = re.sub(r"\.(pdb|cif|ent)(\.gz)?$", "", tok)
                    mm = re.match(r"AF-([A-Z0-9]+)-F\d+", tok)
                    fo.write(">" + (mm.group(1) if mm else tok) + "\n")
                else:
                    fo.write(line)
        os.replace(out, p)
    clean(aa)
    clean(di)
    return aa, di


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--proteomes")
    g.add_argument("--fasta")
    g.add_argument("--structures")
    ap.add_argument("--taxmap")
    ap.add_argument("--taxdump", default=os.path.join(ROOT, "data", "interim", "ancestra", "taxdump"))
    ap.add_argument("--tree", help="use this guide tree instead of building the NCBI one")
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=48)
    ap.add_argument("--preset", default="accurate", choices=["accurate", "fast"])
    ap.add_argument("--families", help="precomputed gene->family TSV (e.g. AFDB cluster membership)")
    ap.add_argument("--extra", default="", help="extra args passed to `ancestra run` (quoted)")
    ap.add_argument("--no-labels", action="store_true", help="skip the all-pairs label table")
    ap.add_argument("--root-split", action="store_true",
                    help="also list same-family pairs in different root HOGs (O(sum family_size^2) rows: "
                         "fine for one family, huge for whole proteomes)")
    a = ap.parse_args(argv)
    os.environ["PATH"] = ENV_BIN + os.pathsep + os.environ["PATH"]
    py = os.path.join(ENV_BIN, "python")
    out = os.path.abspath(a.out)
    os.makedirs(out, exist_ok=True)
    taxmap = load_taxmap(a.taxmap)
    t0 = time.time()
    threedi = None
    if a.proteomes:
        pdir = os.path.abspath(a.proteomes)
        stems = {}
        for f in os.listdir(pdir):
            m = re.match(r"(.+?)\.(fa|fasta|faa|fas|pep)(\.gz)?$", f)
            if m:
                st = m.group(1)
                stems[st] = taxmap.get(st) or int(st.rsplit("_", 1)[1])
        leaf_taxid = stems
    else:
        if a.structures:
            aa, threedi = structures_to_fasta(a.structures, os.path.join(out, "structures"), a.threads)
            recs = read_fasta(aa)
        else:
            recs = read_fasta(a.fasta)
        pdir = os.path.join(out, "proteomes")
        leaf_taxid, nskip = split_by_taxid(recs, taxmap, pdir)
        print(f"[run_ancestra] {len(leaf_taxid)} species; {nskip} sequences skipped (no taxid/short)",
              file=sys.stderr)
    tree = a.tree
    if not tree:
        sys.path.insert(0, os.path.join(ROOT, "tools", "ancestra-mvp-ai"))
        from ancestra.speciestree import ncbi_tree, write_tree
        t = ncbi_tree(leaf_taxid, a.taxdump)
        tree = os.path.join(out, "guide_tree.nwk")
        write_tree(t, tree, lengths=False)
    t_prep = time.time() - t0
    cmd = [py, "-m", "ancestra", "run", "--proteomes", pdir, "--tree", tree, "--out",
           os.path.join(out, "ancestra"), "--threads", str(a.threads), "--preset", a.preset]
    if a.families:
        cmd += ["--families", a.families]
    cmd += a.extra.split()
    env = dict(os.environ, OMP_NUM_THREADS=str(a.threads), NUMBA_NUM_THREADS=str(a.threads))
    t1 = time.time()
    with open(os.path.join(out, "ancestra.log"), "w") as log:
        r = subprocess.run(["/usr/bin/time", "-v"] + cmd, stderr=log, stdout=log, env=env)
    t_run = time.time() - t1
    log_txt = open(os.path.join(out, "ancestra.log")).read()
    m = re.search(r"Maximum resident set size \(kbytes\): (\d+)", log_txt)
    cpu = sum(float(x) for x in re.findall(r"(?:User|System) time \(seconds\): ([\d.]+)", log_txt))
    timing = {"returncode": r.returncode, "prep_s": round(t_prep, 1), "ancestra_wall_s": round(t_run, 1),
              "ancestra_cpu_s": round(cpu, 1), "peak_rss_gb_max_single_process": round(int(m.group(1)) / 1e6, 2) if m else None,
              "n_species": len(leaf_taxid), "threads": a.threads, "preset": a.preset, "cmd": cmd}
    man = os.path.join(out, "ancestra", "manifest.json")
    if os.path.exists(man):
        mf = json.load(open(man))
        timing["n_genes"] = mf.get("n_genes")
        timing["stages"] = {k: {"wall_s": v["wall_s"], "peak_rss_mb": v["peak_rss_mb"],
                                "children_max_rss_mb": v.get("children_max_rss_mb")}
                            for k, v in mf["stages"].items()}
        timing["merge_stats"] = mf["stages"].get("merge", {}).get("stats")
        timing["families_stats"] = {k: v for k, v in mf["stages"].get("families", {}).get("stats", {}).items()
                                    if k != "s0_timings"}
    if r.returncode == 0 and not a.no_labels:
        t2 = time.time()
        lcmd = [py, os.path.join(HERE, "hog_pair_labels.py"), os.path.join(out, "ancestra", "hogs.orthoxml"),
                os.path.join(out, "pair_labels.tsv.gz")]
        if a.root_split:
            lcmd += ["--families-from-run", os.path.join(out, "ancestra")]
        subprocess.run(lcmd, check=True, stderr=open(os.path.join(out, "labels.log"), "w"))
        timing["labels_wall_s"] = round(time.time() - t2, 1)
    if threedi:
        timing["threedi_fasta"] = threedi
    json.dump(timing, open(os.path.join(out, "timing.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in timing.items() if k != "cmd"}, indent=1))
    return r.returncode


if __name__ == "__main__":
    sys.exit(main())
