#!/usr/bin/env python3
"""Interpretation (ANALYSIS_PLAN §5; DESIGN §11.3): match endpoint decisions to §11.3 row ids (ids only, no
free text) and write <dir>/interpretation.json.

Rules: the E4 column is read on beta_k = -beta; the joint interpretation requires spike-in power >= 0.80
(otherwise an E5 acceptance becomes 'inconclusive', E1-E4 decisions and an E5 rejection stand);
E6 > 0 (detected after Holm) with N0 accepted adds the 'duplication' modifier.

  python scripts/endpoints/interpret.py [--dir ROOT/results/endpoints]
"""
import argparse, json, sys
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT / "src"))
from strucevo.stats.endpoints_misc import interpret  # noqa: E402
from strucevo.stats.common import write_json  # noqa: E402


def interpret_dir(d):
    d = Path(d)
    r = {e: json.loads((d / f"{e}.json").read_text()) for e in ("E1", "E2", "E4", "E5", "E6") if (d / f"{e}.json").exists()}
    e5 = r.get("E5", {})
    pw = ((e5.get("power") or {}).get("f=1.5") or {}).get("detect_rate")
    e6 = r.get("E6", {})
    dec = {"E1": r.get("E1", {}).get("decision"), "E2": r.get("E2", {}).get("decision"),
           "E4": r.get("E4", {}).get("decision"),
           "E4_robust_alpha_timescale": bool(r.get("E4", {}).get("robust_alpha_timescale", False)),
           "E5": e5.get("decision"), "E6": e6.get("decision"),
           "E6_positive": (e6.get("E6a") or 0) > 0 or (e6.get("E6b") or 0) > 0, "power": pw}
    out = interpret(dec)
    out["inputs"] = dec
    out["power_rule"] = "power >= 0.80 required for joint interpretation" + ("" if (pw is not None and pw >= 0.80) else
                                                                           " -> NOT met: E5 acceptance treated as inconclusive")
    out["missing_endpoints"] = [e for e in ("E1", "E2", "E4", "E5", "E6") if e not in r]
    out["not_estimable"] = [e for e, v in r.items() if v.get("decision") == "not_estimable"]
    labels = sorted({f for v in r.values() for f in (v.get("flags") or []) if f.startswith("exploratory") or f == "provisional_tuning"})
    out["labels"] = labels
    write_json(out, d / "interpretation.json")
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", default=str(PROJECT / "results" / "endpoints"))
    a = ap.parse_args()
    print(json.dumps(interpret_dir(a.dir), indent=1))
