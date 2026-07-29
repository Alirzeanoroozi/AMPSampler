#!/usr/bin/env python3
"""
Stage 5 variant - pick a panel biased toward outer-membrane PERMEABILITY.

Same inhibitor philosophy as select_candidates.py (catalytic_ok required, never dropped)
but re-weights the ranking so designs with a chance of crossing the outer membrane to
reach the periplasmic target are promoted. Useful when you want candidates that may
work in cellular MIC-rescue WITHOUT a co-administered permeabilizer.

What this does, in order:
  1. Load results/manifest_<T>.csv.
  2. Apply hard gates:
       - catalytic_ok == True               (inhibitor gate)
       - net_charge_pH7.4 >= --min-charge   (default -10, i.e. essentially off)
       - length in [--min-length, --max-length]
       - n_liabilities <= --max-liabilities
       - aggregation_proxy <= --max-aggregation
  3. Re-rank survivors by a permeability-heavy weighted sum (delivery_proxy x2,
     net_charge x1.5, hydrophobic_moment x1, plus the standard inhibitor/quality terms).
  4. Greedy diversity selection (<= --max-identity pairwise) -> top --n.

Outputs (separate from the standard panel so nothing gets clobbered):
  results/selected_<T>_permeable.csv
  results/selected_<T>_permeable.fasta   (record headers include charge for at-a-glance)

Usage:
  python pipeline/stage5_selection/select_permeable.py --target NDM5
  python pipeline/stage5_selection/select_permeable.py --target KPC3 --min-charge -3 --n 16
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from difflib import SequenceMatcher

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

# Permeability-heavy weights. Inhibitor terms still present so we don't pick
# designs that score great on charge but barely touch the active site.
WEIGHTS = {
    "delivery_proxy":         2.0,   # charge + amphipathicity composite (Stage 4)
    "net_charge_pH7.4":       1.5,
    "hydrophobic_moment":     1.0,
    "epitope_recall":         1.0,
    "n_catalytic_contacts":   0.7,
    "n_liabilities":         -0.5,
    "aggregation_proxy":     -0.5,
}


def find_col(keys, name):
    for k in keys:
        if k == name or k.endswith("." + name):
            return k
    return None


def as_float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def identity(a, b):
    return SequenceMatcher(None, a, b).ratio()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True)
    ap.add_argument("--manifest", default=None,
                    help="default: results/manifest_<T>.csv")
    ap.add_argument("--min-charge", type=float, default=-10.0,
                    help="hard floor on net charge at pH 7.4 (default -10 = off; "
                         "use 0 to require non-negative, +1 for definitely cationic)")
    ap.add_argument("--min-length", type=int, default=12)
    ap.add_argument("--max-length", type=int, default=45)
    ap.add_argument("--max-liabilities", type=int, default=2)
    ap.add_argument("--max-aggregation", type=float, default=2.5)
    ap.add_argument("--n", type=int, default=32)
    ap.add_argument("--max-identity", type=float, default=0.8)
    args = ap.parse_args()

    manifest = args.manifest or os.path.join(REPO, "results", f"manifest_{args.target}.csv")
    rows = list(csv.DictReader(open(manifest)))
    if not rows:
        print(f"No rows in {manifest}", file=sys.stderr)
        return 1
    keys = list(rows[0].keys())

    cols = {
        "catalytic_ok":       find_col(keys, "catalytic_ok"),
        "net_charge_pH7.4":   find_col(keys, "net_charge_pH7.4"),
        "length":             find_col(keys, "length"),
        "n_liabilities":      find_col(keys, "n_liabilities"),
        "aggregation_proxy":  find_col(keys, "aggregation_proxy"),
    }

    counters = {"start": len(rows), "no_sequence": 0, "catalytic_ok": 0,
                "min_charge": 0, "length": 0, "n_liabilities": 0,
                "aggregation_proxy": 0, "survived": 0}

    survivors = []
    for r in rows:
        seq = (r.get("sequence") or "").upper()
        if not seq:
            counters["no_sequence"] += 1
            continue
        if cols["catalytic_ok"]:
            if str(r.get(cols["catalytic_ok"], "")).strip().lower() != "true":
                counters["catalytic_ok"] += 1
                continue
        ch = as_float(r.get(cols["net_charge_pH7.4"])) if cols["net_charge_pH7.4"] else None
        if ch is not None and ch < args.min_charge:
            counters["min_charge"] += 1
            continue
        L = as_float(r.get(cols["length"])) if cols["length"] else None
        if L is not None and not (args.min_length <= L <= args.max_length):
            counters["length"] += 1
            continue
        nl = as_float(r.get(cols["n_liabilities"])) if cols["n_liabilities"] else None
        if nl is not None and nl > args.max_liabilities:
            counters["n_liabilities"] += 1
            continue
        ag = as_float(r.get(cols["aggregation_proxy"])) if cols["aggregation_proxy"] else None
        if ag is not None and ag > args.max_aggregation:
            counters["aggregation_proxy"] += 1
            continue
        survivors.append(r)
    counters["survived"] = len(survivors)

    print(f"[{args.target}] gates: {counters['start']} input -> {counters['survived']} pass")
    for k in ("no_sequence", "catalytic_ok", "min_charge", "length",
              "n_liabilities", "aggregation_proxy"):
        if counters[k]:
            print(f"  -{counters[k]:>4} dropped by {k}")

    if not survivors:
        print("\nNothing passed. Suggestions:")
        print("  - lower --min-charge (e.g. -3 for KPC3, whose designs are mostly -3 to -7)")
        print("  - check that results/manifest_<T>.csv has catalytic_ok populated")
        return 1

    # Min-max normalize each present signal and weighted-sum
    present = {}
    for sig, w in WEIGHTS.items():
        col = find_col(keys, sig)
        if not col:
            continue
        vals = [as_float(r.get(col)) for r in survivors]
        vals = [v for v in vals if v is not None]
        if len(vals) >= 2 and max(vals) > min(vals):
            present[sig] = (col, min(vals), max(vals), w)

    for r in survivors:
        score = 0.0
        for sig, (col, lo, hi, w) in present.items():
            v = as_float(r.get(col))
            if v is None:
                continue
            score += w * (v - lo) / (hi - lo)
        r["_permeability_score"] = round(score, 4)
    survivors.sort(key=lambda r: -r["_permeability_score"])

    # Greedy diversity selection
    selected = []
    for r in survivors:
        seq = r["sequence"].upper()
        if all(identity(seq, s["sequence"].upper()) < args.max_identity for s in selected):
            selected.append(r)
        if len(selected) >= args.n:
            break

    out_csv = os.path.join(REPO, "results", f"selected_{args.target}_permeable.csv")
    out_fa = os.path.join(REPO, "results", f"selected_{args.target}_permeable.fasta")
    out_cols = ["_permeability_score"] + [c for c in survivors[0].keys() if c != "_permeability_score"]
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=out_cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(selected)
    with open(out_fa, "w") as fh:
        for r in selected:
            ch = r.get(cols["net_charge_pH7.4"], "") if cols["net_charge_pH7.4"] else ""
            fh.write(f">{r['design_id']} charge={ch} permscore={r['_permeability_score']}\n"
                     f"{r['sequence']}\n")
    print(f"\n[{args.target}] selected {len(selected)} diverse permeability-biased designs")
    print(f"  -> {os.path.relpath(out_csv, REPO)}")
    print(f"  -> {os.path.relpath(out_fa, REPO)}")

    if cols["net_charge_pH7.4"]:
        charges = [as_float(r.get(cols["net_charge_pH7.4"])) for r in selected]
        charges = [c for c in charges if c is not None]
        if charges:
            print(f"  Net charge in panel: min {min(charges):+.1f}, "
                  f"median {sorted(charges)[len(charges)//2]:+.1f}, max {max(charges):+.1f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
