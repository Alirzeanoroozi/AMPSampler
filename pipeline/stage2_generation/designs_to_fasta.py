#!/usr/bin/env python3
"""
Bridge Stage 2 -> Stage 3: collect designed binder sequences into one FASTA whose record
ids ARE the design_id (so every downstream score joins back). Replaces the old
src/utils/csv_to_fasta.py (which ran an example call on import and used a stale path).

BoltzGen writes final_ranked_designs/final_<budget>_designs/ + final_designs_metrics_<budget>.csv
with a designed-sequence column. This reads that CSV (or any CSV with id+sequence columns)
and writes results/stage2_designs/<method>_<target>.fasta.

Usage:
  python designs_to_fasta.py --csv out_NDM5/final_ranked_designs/final_designs_metrics_200.csv \
      --method boltzgen --target NDM5 [--id-col id --seq-col designed_sequence]
"""
import argparse, csv, os

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def pick(cols, candidates):
    for c in candidates:
        if c in cols:
            return c
    raise SystemExit(f"none of {candidates} in columns {cols}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True)
    ap.add_argument("--method", required=True)
    ap.add_argument("--target", required=True)
    ap.add_argument("--id-col", default=None)
    ap.add_argument("--seq-col", default=None)
    ap.add_argument("--prefix", default="", help="optional id prefix for multi-method merges "
                    "(default: none, so design_id == the generator's id == the complex .cif stem, "
                    "which is what Stage 3 active_site_overlap keys on)")
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.csv)))
    if not rows:
        raise SystemExit("empty CSV")
    cols = rows[0].keys()
    id_col = args.id_col or pick(cols, ["id", "design_id", "name", "final_rank"])
    seq_col = args.seq_col or pick(cols, ["designed_sequence", "sequence", "seq", "binder_sequence"])

    out_dir = os.path.join(REPO, "results", "stage2_designs")
    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, f"{args.method}_{args.target}.fasta")
    n = 0
    with open(out, "w") as fh:
        for r in rows:
            seq = (r.get(seq_col) or "").strip().upper()
            if not seq:
                continue
            # design_id MUST equal the complex .cif filename stem so Stage 3 scores join.
            did = f"{args.prefix}{r.get(id_col, n)}"
            fh.write(f">{did}\n{seq}\n")
            n += 1
    print(f"Wrote {n} designs -> results/stage2_designs/{args.method}_{args.target}.fasta")


if __name__ == "__main__":
    main()
