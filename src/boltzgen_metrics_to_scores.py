#!/usr/bin/env python3
"""
Map BoltzGen's `final_designs_metrics_<B>.csv` into a design_id-keyed score CSV
that build_manifest.py can join. The only transformation is renaming the BoltzGen
`id` column to `design_id` (which equals the complex .cif filename stem and is the
primary key across all downstream scores).

Usage:
  python boltzgen_metrics_to_scores.py --metrics <final_designs_metrics_<B>.csv> --out <scores_dir>/boltzgen_metrics.csv
"""
import argparse, csv, os, sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metrics", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.metrics)))
    if not rows:
        print(f"empty metrics CSV: {args.metrics}", file=sys.stderr)
        return 1
    src = "id" if "id" in rows[0] else "design_id" if "design_id" in rows[0] else None
    if src is None:
        print(f"no 'id' column in {args.metrics}; have {list(rows[0])}", file=sys.stderr)
        return 1

    cols = ["design_id"] + [c for c in rows[0].keys() if c not in ("id", "design_id")]
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=cols)
        w.writeheader()
        for r in rows:
            out = {c: r.get(c, "") for c in cols}
            out["design_id"] = r[src]
            w.writerow(out)
    print(f"Wrote {len(rows)} BoltzGen-metrics rows -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
