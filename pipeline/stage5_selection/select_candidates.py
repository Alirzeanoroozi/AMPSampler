#!/usr/bin/env python3
"""
Stage 5 - Select a diverse wet-lab panel from the filtered, ranked designs.

Greedy diversity selection: walk designs best-rank-first, keep one only if it is below
`max_identity` to everything already kept (drops near-duplicate motifs). Targets a panel
of ~24-48 per target.

Usage:
  python select_candidates.py --target NDM5 [--n 32] [--max-identity 0.8]
Inputs : results/filtered_<T>.csv (passing+ranked) | falls back to filtered_<T>_all.csv | manifest_<T>.csv
Outputs: results/selected_<T>.csv, results/selected_<T>.fasta
"""
from __future__ import annotations

import argparse, csv, os
from difflib import SequenceMatcher

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def load_ranked(tkey):
    for name in (f"filtered_{tkey}.csv", f"filtered_{tkey}_all.csv", f"manifest_{tkey}.csv"):
        p = os.path.join(REPO, "results", name)
        if os.path.exists(p):
            rows = list(csv.DictReader(open(p)))
            if "_pass" in (rows[0].keys() if rows else []):
                rows = [r for r in rows if r.get("_pass") in ("True", "true", True)]
            # keep ranking order if rank_score present
            if rows and rows[0].get("rank_score", "") not in ("", None):
                rows.sort(key=lambda r: float(r["rank_score"]) if r.get("rank_score") not in ("", None) else -1e9,
                          reverse=True)
            return rows, name
    return [], None


def identity(a, b):
    return SequenceMatcher(None, a, b).ratio()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True)
    ap.add_argument("--n", type=int, default=32)
    ap.add_argument("--max-identity", type=float, default=0.8)
    args = ap.parse_args()

    rows, src = load_ranked(args.target)
    if not rows:
        print(f"No ranked designs found for {args.target} (run Stages 2-4 first).")
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
    cols = list(rows[0].keys())
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader(); w.writerows(selected)
    with open(out_fa, "w") as fh:
        for r in selected:
            fh.write(f">{r['design_id']} method={r.get('method','')} rank={r.get('rank_score','')}\n{r['sequence']}\n")
    print(f"[{args.target}] selected {len(selected)} diverse designs (<{args.max_identity} identity) "
          f"from {len(rows)} ranked ({src}) -> results/selected_{args.target}.csv/.fasta")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
