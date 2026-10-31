"""V11 independent re-derivation: shared input loading (only project DATA helpers are used:
strucevo.stats.common.load_pairs / load_branches; no endpoint/driver code)."""
import os, json
from pathlib import Path
import numpy as np
import polars as pl
import yaml

os.environ.setdefault("CONFIRM_RUN", "1")
from strucevo.stats.common import load_pairs, load_branches  # data loading only

ROOT = Path(__file__).resolve().parents[2]
TUN = yaml.safe_load((ROOT / "config/tuning.yaml").read_text())
PRIMARY = TUN["primary_alignment"]
SPLIT = "confirm"
OUT = ROOT / "results" / "verify"


def binlook(edges, vals, x):
    edges = np.asarray(edges, float); vals = np.asarray(vals, float)
    x = np.nan_to_num(np.asarray(x, float), nan=0.0)
    b = np.searchsorted(edges, x, side="right") - 1
    return vals[np.clip(b, 0, len(vals) - 1)]


def pairs(aln=PRIMARY):
    """All confirm rows of the alignment table, with kappa(bin of aa_ml) x Fisher SE applied (§2.3)."""
    df = load_pairs(aln, SPLIT)
    assert set(df["split"].unique().to_list()) == {"confirm"}
    kap = TUN["kappa"][aln]
    x = df["aa_ml"].to_numpy().astype(float)
    df = df.with_columns(
        pl.Series("aa_ml_se_k", df["aa_ml_se"].to_numpy() * binlook(kap["bin_edges"], kap["aa_ml"], x)),
        pl.Series("tdi_ml_se_k", df["tdi_ml_se"].to_numpy() * binlook(kap["bin_edges"], kap["tdi_ml"], x)))
    return df


def hogcols():
    return pl.read_parquet(ROOT / "data/interim" / f"hogs_endpoint_cols_{SPLIT}.parquet")


def hogs():
    h = pl.read_parquet(ROOT / "data/interim/hogs.parquet")
    return h


def endpoint(name):
    return json.loads((ROOT / "results/endpoints" / f"{name}.json").read_text())


def wjson(obj, path):
    def j(o):
        if isinstance(o, dict): return {str(k): j(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)): return [j(v) for v in o]
        if isinstance(o, (np.floating, float)): return None if not np.isfinite(o) else float(o)
        if isinstance(o, np.integer): return int(o)
        if isinstance(o, np.bool_): return bool(o)
        return o
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(j(obj), indent=2))
