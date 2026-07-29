#!/usr/bin/env python3
"""
Stage 4 - Apply inhibitor-binder gates and rank the manifest.

Robust to missing columns: each gate is applied ONLY if its column exists in the
manifest, so this runs at any stage of completion (the report says which gates were
active). Ranking is a transparent weighted sum over min-max-normalised signals that
are present. AMP-database similarity is deliberately NOT a criterion.

Usage:
  python apply_filters.py --manifest results/manifest_NDM5.csv [--filters filters.json]
Outputs: results/filtered_<T>.csv (passing, ranked) and results/filtered_<T>_all.csv (all + reasons)
"""
from __future__ import annotations

import argparse, csv, json, os

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))


def resolve(row, candidates):
    """Find the first present, non-empty value whose column matches a candidate
    (exact or as a suffix after build_manifest's 'tag.' prefixing)."""
    if isinstance(candidates, str):
        candidates = [candidates]
    for cand in candidates:
        for col, val in row.items():
            if (col == cand or col.endswith("." + cand)) and str(val).strip() != "":
                return val
    return None


def as_float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def gate_pass(row, spec):
    val = resolve(row, spec["column"])
    if val is None:
        return None  # gate inactive (no data)
    if "equals" in spec:
        return str(val) == str(spec["equals"])
    if "not_equals" in spec:
        return str(val) != str(spec["not_equals"])
    f = as_float(val)
    if f is None:
        return None
    if "min" in spec and f < spec["min"]:
        return False
    if "max" in spec and f > spec["max"]:
        return False
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--filters", default=os.path.join(HERE, "filters.json"))
    args = ap.parse_args()

    cfg = json.load(open(args.filters))
    gates = {k: v for k, v in cfg["gates"].items()}
    weights = cfg["ranking_weights"]
    rows = list(csv.DictReader(open(args.manifest)))
    tkey = os.path.splitext(os.path.basename(args.manifest))[0].replace("manifest_", "")

    # which gates have data?
    active = [g for g, spec in gates.items()
              if any(gate_pass(r, spec) is not None for r in rows)]

    # evaluate gates
    for r in rows:
        reasons = []
        for g in active:
            res = gate_pass(r, gates[g])
            if res is False:
                reasons.append(gates[g]["reason"])
        r["_fail_reasons"] = "; ".join(reasons)
        r["_pass"] = (reasons == [])

    # ranking: min-max normalise each present signal, weighted sum
    present = {}
    for sig in weights:
        vals = [as_float(resolve(r, sig)) for r in rows]
        vals = [v for v in vals if v is not None]
        if len(vals) >= 2 and max(vals) > min(vals):
            present[sig] = (min(vals), max(vals))
    for r in rows:
        score, used = 0.0, 0
        for sig, (lo, hi) in present.items():
            v = as_float(resolve(r, sig))
            if v is None:
                continue
            norm = (v - lo) / (hi - lo)
            score += weights[sig] * norm
            used += 1
        r["rank_score"] = round(score, 4) if used else ""

    rows.sort(key=lambda r: (r["_pass"], r["rank_score"] if r["rank_score"] != "" else -1e9), reverse=True)

    base_cols = [c for c in rows[0].keys() if not c.startswith("_") and c != "rank_score"] if rows else []
    cols = ["rank_score"] + base_cols + ["_pass", "_fail_reasons"]
    out_all = os.path.join(REPO, "results", f"filtered_{tkey}_all.csv")
    out_pass = os.path.join(REPO, "results", f"filtered_{tkey}.csv")
    os.makedirs(os.path.dirname(out_all), exist_ok=True)
    for path, subset in ((out_all, rows), (out_pass, [r for r in rows if r["_pass"]])):
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
            w.writeheader()
            w.writerows(subset)

    n_pass = sum(r["_pass"] for r in rows)
    print(f"[{tkey}] {len(rows)} designs; gates active: {active or 'none (run Stages 2-3 first)'}")
    print(f"  ranking signals used: {list(present) or 'none'}")
    print(f"  PASS {n_pass}/{len(rows)} -> {os.path.relpath(out_pass, REPO)} ; all+reasons -> {os.path.relpath(out_all, REPO)}")


if __name__ == "__main__":
    main()
