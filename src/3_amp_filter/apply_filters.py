#!/usr/bin/env python3
"""Keep designs that pass filters.json hard gates (no ranking).

Gates mix AMP classifier scores with the developability / inhibitor floors from
select_permeable.py (catalytic_ok, net charge, length, liabilities, aggregation).
A missing value (e.g. HydrAMP not scored, or catalytic_ok not yet computed)
skips that gate for that row.

Usage:
  python apply_filters.py --manifest results/manifest_NDM5.csv
  python apply_filters.py --manifest results/manifest_KPC3.csv --out results/filtered_KPC3.csv

Writes the passing set as CSV + FASTA for the structure-prediction step, and
an _all table with pass/fail reasons.
"""
from __future__ import annotations

import argparse
import csv
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))


def resolve(row, candidates):
    """First present, non-empty value whose column matches a candidate name."""
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


def _norm(v):
    return str(v).strip()


def _boolish_equal(got, expected):
    """Match True/False the way select_permeable.py does (case-insensitive)."""
    exp = expected.lower()
    if exp in ("true", "false"):
        return got.lower() == exp
    return got == expected


def gate_pass(row, spec):
    val = resolve(row, spec["column"])
    if val is None:
        return None  # no data on this row -> skip this gate
    got = _norm(val)
    if "equals" in spec:
        return _boolish_equal(got, _norm(spec["equals"]))
    if "not_equals" in spec:
        return not _boolish_equal(got, _norm(spec["not_equals"]))
    f = as_float(val)
    if f is None:
        return None
    if "min" in spec and f < spec["min"]:
        return False
    if "max" in spec and f > spec["max"]:
        return False
    return True


def target_key(path):
    stem = os.path.splitext(os.path.basename(path))[0]
    for prefix in ("manifest_", "filtered_"):
        if stem.startswith(prefix):
            return stem[len(prefix):]
    return stem


def write_fasta(rows, path):
    n = 0
    with open(path, "w") as fh:
        for r in rows:
            did = (r.get("design_id") or r.get("seq_id") or "").strip()
            seq = (r.get("sequence") or "").strip()
            if not did or not seq:
                continue
            fh.write(f">{did}\n{seq}\n")
            n += 1
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--filters", default=os.path.join(HERE, "filters.json"))
    ap.add_argument(
        "--out",
        default=None,
        help="passing CSV (default: results/filtered_<TARGET>.csv)",
    )
    args = ap.parse_args()

    cfg = json.load(open(args.filters))
    gates = cfg["gates"]
    rows = list(csv.DictReader(open(args.manifest)))
    tkey = target_key(args.out or args.manifest)

    for r in rows:
        reasons = []
        if not (r.get("sequence") or "").strip():
            reasons.append("no sequence")
        for spec in gates.values():
            res = gate_pass(r, spec)
            if res is False:
                reasons.append(spec["reason"])
        r["_fail_reasons"] = "; ".join(reasons)
        r["_pass"] = reasons == []

    passing = [r for r in rows if r["_pass"]]

    out_pass = args.out or os.path.join(REPO, "results", f"filtered_{tkey}.csv")
    out_all = os.path.splitext(out_pass)[0] + "_all.csv"
    if out_all == out_pass:
        out_all = out_pass.replace(".csv", "_all.csv")
    out_fa = os.path.splitext(out_pass)[0] + ".fasta"
    os.makedirs(os.path.dirname(os.path.abspath(out_pass)) or ".", exist_ok=True)

    base_cols = [c for c in (rows[0].keys() if rows else []) if not c.startswith("_")]
    cols = base_cols + ["_pass", "_fail_reasons"]
    for path, subset in ((out_all, rows), (out_pass, passing)):
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
            w.writeheader()
            w.writerows(subset)
    n_fa = write_fasta(passing, out_fa)

    print(f"[{tkey}] {len(rows)} designs; gates: {list(gates)}")
    print(
        f"  PASS {len(passing)}/{len(rows)} -> {os.path.relpath(out_pass, REPO)} "
        f"+ {n_fa} sequences in {os.path.relpath(out_fa, REPO)}"
    )
    print(f"  all+reasons -> {os.path.relpath(out_all, REPO)}")


if __name__ == "__main__":
    main()
