#!/usr/bin/env python3
"""Parse Boltz-2 predictions into a design_id-keyed score table.

Boltz-2 writes, per prediction <name>:
  predictions/<name>/confidence_<name>_model_0.json   (iptm, ptm, complex_plddt, ...)
  predictions/<name>/affinity_<name>.json             (affinity_pred_value, affinity_probability_binary)

Usage (from AMPBinderDesign):
  python src/4_structure_prediction/parse_boltz2.py --target NDM5
  python src/4_structure_prediction/parse_boltz2.py --target KPC3
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))


def find(d, *keys):
    for k in keys:
        if k in d:
            return d[k]
    return ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--boltz_out",
        default=os.path.join(REPO, "boltz_results"),
        help="Boltz-2 output dir (searched recursively for confidence_*.json)",
    )
    ap.add_argument("--target", choices=["NDM5", "KPC3"], default=None)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    conf_files = glob.glob(os.path.join(args.boltz_out, "**", "confidence_*.json"), recursive=True)
    rows = []
    for cf in sorted(conf_files):
        name = os.path.basename(os.path.dirname(cf))
        if args.target and args.target not in name:
            continue
        try:
            conf = json.load(open(cf))
        except Exception:
            continue
        aff_path = glob.glob(os.path.join(os.path.dirname(cf), "affinity_*.json"))
        aff = {}
        if aff_path:
            try:
                aff = json.load(open(aff_path[0]))
            except Exception:
                aff = {}
        rows.append({
            "design_id": name,
            "boltz2_iptm": find(conf, "iptm", "complex_iptm"),
            "boltz2_ptm": find(conf, "ptm"),
            "boltz2_plddt": find(conf, "complex_plddt", "plddt"),
            "boltz2_affinity_pred_value": find(aff, "affinity_pred_value"),
            "boltz2_affinity_prob_binary": find(aff, "affinity_probability_binary"),
        })

    if args.out:
        out = args.out
    elif args.target:
        out = os.path.join(REPO, "results", f"boltz2_{args.target}.csv")
    else:
        out = os.path.join(REPO, "results", "boltz2_all.csv")
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)

    cols = ["design_id", "boltz2_iptm", "boltz2_ptm", "boltz2_plddt",
            "boltz2_affinity_pred_value", "boltz2_affinity_prob_binary"]
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    tag = args.target or "all"
    print(f"[{tag}] Parsed {len(rows)} Boltz-2 predictions -> {os.path.relpath(out, REPO)}")


if __name__ == "__main__":
    main()
