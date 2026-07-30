#!/usr/bin/env python3
"""
Stage 5 - Select AMP binders for the wet-lab panel.

Primary path (default): keep only designs classified as AMP by AMPScanner
and/or Macrel from results/classifiers/, ranked by AMPScanner probability then
BoltzGen iPTM, with greedy sequence-diversity filtering.

Fallback: if no classifier rows match the target, fall back to filtered_/manifest_
ranking (inhibitor panel) — same as the old behaviour.

Usage:
  python select_candidates.py --target NDM5 [--n 32] [--max-identity 0.8]
  python select_candidates.py --target KPC3 --require both   # AMPScanner AND Macrel

Outputs: results/selected_<T>.csv, results/selected_<T>.fasta
"""
from __future__ import annotations

import argparse
import csv
import os
from difflib import SequenceMatcher

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
CLS_DIR = os.path.join(REPO, "results", "classifiers")


def identity(a, b):
    return SequenceMatcher(None, a, b).ratio()


def as_float(v, default=-1e9):
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def is_true(v):
    return str(v).strip().lower() in ("true", "1", "yes")


def is_amp(row, mode: str) -> bool:
    scanner = str(row.get("ampscanner_class", "")).strip() == "AMP"
    macrel = is_true(row.get("macrel_is_amp")) or is_true(row.get("macrel_is_amp_bool"))
    if mode == "ampscanner":
        return scanner
    if mode == "macrel":
        return macrel
    if mode == "both":
        return scanner and macrel
    return scanner or macrel  # either


def load_amp_rows(tkey: str, mode: str):
    path = os.path.join(CLS_DIR, "merged_binder_selection.csv")
    if not os.path.exists(path):
        return [], None
    rows = [
        r for r in csv.DictReader(open(path))
        if r.get("target") == tkey and is_amp(r, mode)
    ]
    # prefer high AMP confidence, then binder–target iPTM
    rows.sort(
        key=lambda r: (
            as_float(r.get("ampscanner_prob")),
            as_float(r.get("bg_design_to_target_iptm")),
            as_float(r.get("boltz2_iptm")),
        ),
        reverse=True,
    )
    # normalise id field used downstream
    for r in rows:
        r.setdefault("design_id", r.get("bg_id") or r.get("seq_id") or "")
        r.setdefault("method", "boltzgen")
        r.setdefault("rank_score", r.get("ampscanner_prob", ""))
    return rows, "classifiers/merged_binder_selection.csv (AMP only)"


def load_ranked_fallback(tkey):
    for name in (f"filtered_{tkey}.csv", f"filtered_{tkey}_all.csv", f"manifest_{tkey}.csv"):
        p = os.path.join(REPO, "results", name)
        if not os.path.exists(p):
            continue
        rows = list(csv.DictReader(open(p)))
        if rows and "_pass" in rows[0]:
            rows = [r for r in rows if is_true(r.get("_pass"))]
        if rows and rows[0].get("rank_score", "") not in ("", None):
            rows.sort(
                key=lambda r: as_float(r.get("rank_score")),
                reverse=True,
            )
        return rows, name
    return [], None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True)
    ap.add_argument("--n", type=int, default=50)
    ap.add_argument("--max-identity", type=float, default=0.8)
    ap.add_argument(
        "--require",
        choices=("either", "ampscanner", "macrel", "both"),
        default="either",
        help="AMP definition: AMPScanner and/or Macrel (default: either)",
    )
    ap.add_argument(
        "--from-filtered",
        action="store_true",
        help="ignore classifiers; select from filtered/manifest (old inhibitor path)",
    )
    args = ap.parse_args()

    if args.from_filtered:
        rows, src = load_ranked_fallback(args.target)
    else:
        rows, src = load_amp_rows(args.target, args.require)
        if not rows:
            print(
                f"[{args.target}] no AMP rows in classifiers (require={args.require}); "
                "falling back to filtered/manifest"
            )
            rows, src = load_ranked_fallback(args.target)

    if not rows:
        print(f"No designs found for {args.target}.")
        return 1

    selected = []
    for r in rows:
        seq = (r.get("sequence") or "").upper()
        if not seq:
            continue
        if all(identity(seq, s.get("sequence", "").upper()) < args.max_identity for s in selected):
            selected.append(r)
        if len(selected) >= args.n:
            break

    out_csv = os.path.join(REPO, "results", f"selected_{args.target}.csv")
    out_fa = os.path.join(REPO, "results", f"selected_{args.target}.fasta")
    cols = list(dict.fromkeys(
        ["design_id", "seq_id", "target", "sequence", "length",
         "ampscanner_class", "ampscanner_prob", "macrel_is_amp", "macrel_amp_prob",
         "macrel_hemolytic", "macrel_hemo_prob", "toxinpred_class",
         "bg_design_to_target_iptm", "boltz2_iptm", "rank_score", "method"]
        + [c for c in rows[0].keys() if c not in (
            "design_id", "seq_id", "target", "sequence", "length",
            "ampscanner_class", "ampscanner_prob", "macrel_is_amp", "macrel_amp_prob",
            "macrel_hemolytic", "macrel_hemo_prob", "toxinpred_class",
            "bg_design_to_target_iptm", "boltz2_iptm", "rank_score", "method"
        )]
    ))
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(selected)
    with open(out_fa, "w") as fh:
        for r in selected:
            did = r.get("design_id") or r.get("seq_id")
            fh.write(
                f">{did} amp={r.get('ampscanner_class', '')}:"
                f"{r.get('ampscanner_prob', '')} "
                f"macrel={r.get('macrel_is_amp', '')} "
                f"iptm={r.get('bg_design_to_target_iptm', '')}\n"
                f"{r['sequence']}\n"
            )
    print(
        f"[{args.target}] selected {len(selected)} AMP binders "
        f"(require={args.require}, <{args.max_identity} identity) "
        f"from {len(rows)} candidates ({src}) -> results/selected_{args.target}.csv/.fasta"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
