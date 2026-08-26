#!/usr/bin/env python3
"""Write Boltz-2 YAMLs for AMP-filtered binders (does not run prediction).

Each filtered FASTA record becomes one config:
  chain A = target design-domain sequence
  chain B = binder
  job name = design_id  (so Boltz-2 / parse_boltz2.py key on the same id)

Usage (from AMPBinderDesign):
  python src/4_structure_prediction/prepare_configs.py
"""
from __future__ import annotations

import argparse
import csv
import os

from Bio import SeqIO

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))

TARGETS = ("NDM5", "KPC3")
TARGET_CHAIN = "A"
BINDER_CHAIN = "B"


def read_seq(path: str) -> str:
    rec = next(SeqIO.parse(path, "fasta"))
    return "".join(str(rec.seq).split()).upper()


def read_binders(path: str) -> list[tuple[str, str]]:
    out = []
    for rec in SeqIO.parse(path, "fasta"):
        seq = "".join(str(rec.seq).split()).upper()
        if not seq:
            continue
        out.append((rec.id, seq))
    return out


def render_yaml(target_seq: str, binder_seq: str) -> str:
    return (
        "version: 1\n"
        "sequences:\n"
        "  - protein:\n"
        f"      id: {TARGET_CHAIN}\n"
        f"      sequence: {target_seq}\n"
        "      msa: empty\n"
        "  - protein:\n"
        f"      id: {BINDER_CHAIN}\n"
        f"      sequence: {binder_seq}\n"
        "      msa: empty\n"
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
        "--out-dir",
        default=os.path.join(REPO, "boltz_inputs"),
        help="Boltz-2 YAML directory (slurm: boltz predict boltz_inputs)",
    )
    ap.add_argument(
        "--manifest",
        default=os.path.join(REPO, "results", "boltz_config_manifest.csv"),
    )
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    os.makedirs(os.path.dirname(args.manifest) or ".", exist_ok=True)

    rows = []
    for target in args.targets:
        fasta = os.path.join(args.filtered_dir, f"filtered_{target}.fasta")
        if not os.path.isfile(fasta):
            raise FileNotFoundError(fasta)
        target_fa = os.path.join(REPO, "targets", target, "design_domain.fasta")
        if not os.path.isfile(target_fa):
            raise FileNotFoundError(target_fa)
        target_seq = read_seq(target_fa)
        binders = read_binders(fasta)
        if not binders:
            raise ValueError(f"{fasta}: no sequences")

        for design_id, binder_seq in binders:
            yaml_name = f"{design_id}.yaml"
            yaml_path = os.path.join(args.out_dir, yaml_name)
            with open(yaml_path, "w") as fh:
                fh.write(render_yaml(target_seq, binder_seq))
            rows.append(
                {
                    "design_id": design_id,
                    "stem": design_id,
                    "target": target,
                    "binder_length": len(binder_seq),
                    "binder_sequence": binder_seq,
                    "receptor_sequence": target_seq,
                    "yaml_name": yaml_name,
                    "yaml_path": os.path.relpath(yaml_path, REPO),
                    "target_chain": TARGET_CHAIN,
                    "binder_chain": BINDER_CHAIN,
                    "filtered_fasta": os.path.relpath(fasta, REPO),
                }
            )
        print(f"[{target}] {len(binders)} configs from {os.path.relpath(fasta, REPO)}")

    with open(args.manifest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    print(f"Wrote {len(rows)} YAMLs -> {os.path.relpath(args.out_dir, REPO)}")
    print(f"Manifest -> {os.path.relpath(args.manifest, REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
