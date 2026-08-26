#!/usr/bin/env python3
"""Build the structure-prediction manifest from AMP-filtered FASTAs.

One row per design_id from results/filtered_<TARGET>.fasta. Left-joins any
score CSVs that have a design_id column (filtered AMP table, Boltz-2, ipSAE,
active-site overlap, ...).

Usage (from AMPBinderDesign):
  python src/4_structure_prediction/build_manifest.py
  python src/4_structure_prediction/build_manifest.py --scores results
"""
from __future__ import annotations

import argparse
import csv
import glob
import os

from Bio import SeqIO

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))

TARGETS = ("NDM5", "KPC3")
SKIP_SCORE_STEMS = {
    "merged_classifier_results",
    "developability_NDM5",
    "developability_KPC3",
    "manifest_NDM5",
    "manifest_KPC3",
    "filtered_NDM5_all",
    "filtered_KPC3_all",
    "structure_manifest_NDM5",
    "structure_manifest_KPC3",
}


def read_designs(fasta_path: str, target: str) -> dict:
    rows = {}
    for rec in SeqIO.parse(fasta_path, "fasta"):
        seq = "".join(str(rec.seq).split()).upper()
        did = rec.id
        rows[did] = {
            "design_id": did,
            "target": target,
            "sequence": seq,
            "length": len(seq),
            "filtered_fasta": os.path.relpath(fasta_path, REPO),
        }
    return rows


def merge_csv(rows: dict, csv_path: str, joined_cols: list) -> None:
    try:
        reader = list(csv.DictReader(open(csv_path)))
    except Exception:
        return
    if not reader or "design_id" not in reader[0]:
        return
    tag = os.path.splitext(os.path.basename(csv_path))[0]
    reserved = {"design_id", "target", "sequence", "length", "filtered_fasta"}
    for r in reader:
        did = r.get("design_id", "").strip()
        if not did or did not in rows:
            continue
        for k, v in r.items():
            if k in reserved:
                continue
            col = k if k not in rows[did] else f"{tag}.{k}"
            rows[did][col] = v
            if col not in joined_cols:
                joined_cols.append(col)


def score_files(target: str, scores_dir: str | None) -> list[str]:
    named = [
        os.path.join(REPO, "results", f"filtered_{target}.csv"),
        os.path.join(REPO, "results", "boltz_config_manifest.csv"),
        os.path.join(REPO, "results", f"boltz2_{target}.csv"),
        os.path.join(REPO, "results", f"active_site_overlap_{target}.csv"),
        os.path.join(REPO, "results", f"ipsae_{target}.csv"),
    ]
    extra = []
    if scores_dir and os.path.isdir(scores_dir):
        extra = sorted(glob.glob(os.path.join(scores_dir, "*.csv")))
    seen = set()
    out = []
    for path in named + extra:
        if not os.path.isfile(path):
            continue
        stem = os.path.splitext(os.path.basename(path))[0]
        if stem in SKIP_SCORE_STEMS:
            continue
        if path in seen:
            continue
        seen.add(path)
        out.append(path)
    return out


def write_manifest(rows: dict, extra_cols: list, out: str, target: str) -> None:
    base = ["design_id", "target", "length", "sequence", "filtered_fasta"]
    cols = base + [c for c in extra_cols if c not in base]
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for did in sorted(rows):
            w.writerow({c: rows[did].get(c, "") for c in cols})
    print(
        f"[{target}] structure manifest: {len(rows)} designs x {len(cols)} columns "
        f"-> {os.path.relpath(out, REPO)}"
    )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--targets", nargs="+", default=list(TARGETS))
    ap.add_argument(
        "--filtered-dir",
        default=os.path.join(REPO, "results"),
        help="directory with filtered_<TARGET>.fasta",
    )
    ap.add_argument(
        "--scores",
        default=None,
        help="optional extra directory of design_id CSVs to join",
    )
    ap.add_argument(
        "--out-dir",
        default=os.path.join(REPO, "results"),
        help="writes structure_manifest_<TARGET>.csv here",
    )
    args = ap.parse_args()

    for target in args.targets:
        fasta = os.path.join(args.filtered_dir, f"filtered_{target}.fasta")
        if not os.path.isfile(fasta):
            raise FileNotFoundError(fasta)
        rows = read_designs(fasta, target)
        if not rows:
            raise ValueError(f"{fasta}: no sequences")
        extra_cols = []
        for path in score_files(target, args.scores):
            before = len(extra_cols)
            merge_csv(rows, path, extra_cols)
            if len(extra_cols) > before:
                print(f"  [{target}] joined {os.path.relpath(path, REPO)}")
        out = os.path.join(args.out_dir, f"structure_manifest_{target}.csv")
        write_manifest(rows, extra_cols, out, target)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
