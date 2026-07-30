#!/usr/bin/env python3
"""
Stage 3 - Build the canonical design manifest (the bookkeeping the old pipeline lacked).

The old repo stored Boltz2 scores keyed by an integer `peptide_num` with no link back
to sequences (FASTAs used `seq_rank_N`), so a score could not be traced to a design.
This tool makes `design_id` the single primary key and left-joins every score table on
it, producing one row per design with sequence + method + all metrics.

Inputs:
  --designs : a FASTA, or a directory of FASTAs named <method>_<target>.fasta
              (record id -> design_id; sequence -> sequence; filename -> method/target)
  --scores  : directory of CSVs that each contain a `design_id` column
              (e.g. <T>_active_site_overlap.csv, boltz2_<T>.csv, af_<T>.csv, ddg_<T>.csv,
               toxinpred, hemolysis, developability ...). All are merged on design_id.

Output: results/manifest_<target>.csv  (one row per design, all scores joined)

Runs in the base conda env. This is the table Stage 4 filters and Stage 5 selects from.
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import sys

from Bio import SeqIO

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def read_designs(designs_path):
    rows = {}
    files = ([designs_path] if os.path.isfile(designs_path)
             else sorted(glob.glob(os.path.join(designs_path, "*.fasta"))))
    for f in files:
        base = os.path.splitext(os.path.basename(f))[0]
        parts = base.split("_")
        method = parts[0] if parts else "unknown"
        for rec in SeqIO.parse(f, "fasta"):
            did = rec.id
            rows[did] = {"design_id": did, "method": method,
                         "sequence": str(rec.seq).upper(), "length": len(rec.seq)}
    return rows


def merge_scores(rows, scores_dir):
    joined_cols = []
    for csv_path in sorted(glob.glob(os.path.join(scores_dir, "*.csv"))):
        try:
            reader = list(csv.DictReader(open(csv_path)))
        except Exception:
            continue
        if not reader or "design_id" not in reader[0]:
            continue
        tag = os.path.splitext(os.path.basename(csv_path))[0]
        for r in reader:
            did = r["design_id"]
            if did not in rows:
                rows[did] = {"design_id": did, "method": "", "sequence": "", "length": ""}
            for k, v in r.items():
                if k == "design_id":
                    continue
                col = k if k not in ("method", "sequence", "length") else f"{tag}.{k}"
                rows[did][col] = v
                if col not in joined_cols:
                    joined_cols.append(col)
    return joined_cols


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True)
    ap.add_argument("--designs", required=True)
    ap.add_argument("--scores", default=None)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    rows = read_designs(args.designs)
    extra = merge_scores(rows, args.scores) if args.scores else []

    base_cols = ["design_id", "method", "length", "sequence"]
    cols = base_cols + [c for c in extra if c not in base_cols]
    out = args.out or os.path.join(REPO, "results", f"manifest_{args.target}.csv")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for did in sorted(rows):
            w.writerow({c: rows[did].get(c, "") for c in cols})
    print(f"[{args.target}] manifest: {len(rows)} designs x {len(cols)} columns -> {os.path.relpath(out, REPO)}")
    if extra:
        print(f"  joined score columns: {', '.join(extra[:12])}{' ...' if len(extra) > 12 else ''}")


if __name__ == "__main__":
    main()
