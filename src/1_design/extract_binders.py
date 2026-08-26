#!/usr/bin/env python3
"""
Walk every CIF under boltzgen_outputs/ (all subdirectories) and write one FASTA
per target. A file is assigned to NDM5.fasta or KPC3.fasta if its path/name
contains that target string.

Binder = shortest protein chain (design outputs: chain B; inverse-fold: chain A).
Native CIFs (*_native.cif) are skipped.

Usage (from AMPBinderDesign):
  conda activate ampbinder
  python src/1_design/boltzgen_extract_binders.py
"""
from __future__ import annotations

import argparse
import os
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO
from Bio.PDB import MMCIFParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))

THREE_TO_ONE = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
    "MSE": "M", "SEP": "S", "TPO": "T", "PTR": "Y", "CSO": "C", "KCX": "K",
    "UNK": "X",
}


def chain_sequence(chain) -> str:
    letters = []
    for res in chain:
        if res.id[0] != " ":
            continue
        letters.append(THREE_TO_ONE.get(res.get_resname().strip(), "X"))
    return "".join(letters)


def binder_sequence(cif_path: str, parser: MMCIFParser) -> Tuple[str, str]:
    """Return (binder_seq, chain_id) using the shortest amino-acid chain."""
    structure = parser.get_structure(os.path.basename(cif_path), cif_path)
    model = next(structure.get_models())
    chains: List[Tuple[str, str]] = []
    for chain in model:
        seq = chain_sequence(chain)
        if seq:
            chains.append((chain.id, seq))
    if not chains:
        raise ValueError(f"no protein chains in {cif_path}")
    chain_id, seq = min(chains, key=lambda x: len(x[1]))
    return seq, chain_id


def target_from_path(path: str, targets: List[str]) -> Optional[str]:
    """Assign a CIF to a target if NDM5 or KPC3 appears in the file path."""
    blob = path.replace("\\", "/").upper()
    hits = [t for t in targets if t.upper() in blob]
    if not hits:
        return None
    return max(hits, key=len)


def list_cifs(directory: str, skip_native: bool = True) -> List[str]:
    """All .cif files under directory, recursively."""
    out = []
    if not os.path.isdir(directory):
        return out
    for root, _dirs, files in os.walk(directory):
        for name in files:
            if not name.endswith(".cif"):
                continue
            if skip_native and "_native" in name:
                continue
            out.append(os.path.join(root, name))
    out.sort()
    return out


def extract_all(
    directory: str,
    targets: List[str],
) -> Dict[str, List[SeqRecord]]:
    parser = MMCIFParser(QUIET=True)
    records: Dict[str, List[SeqRecord]] = {t: [] for t in targets}
    seen_ids: Dict[str, Dict[str, int]] = {t: {} for t in targets}
    paths = list_cifs(directory)
    skipped = 0
    unmatched = 0
    print(f"Found {len(paths)} CIFs under {directory}")

    for i, path in enumerate(paths, 1):
        target = target_from_path(path, targets)
        if target is None:
            unmatched += 1
            continue
        stem = os.path.splitext(os.path.basename(path))[0]
        rec_id = stem
        if rec_id in seen_ids[target]:
            seen_ids[target][rec_id] += 1
            rec_id = f"{stem}_{seen_ids[target][stem]}"
        else:
            seen_ids[target][stem] = 1
        try:
            seq, chain_id = binder_sequence(path, parser)
        except Exception as exc:
            print(f"  skip {stem}: {exc}")
            skipped += 1
            continue
        if not seq:
            print(f"  skip {stem}: empty binder sequence")
            skipped += 1
            continue
        rel = os.path.relpath(path, directory)
        records[target].append(
            SeqRecord(
                Seq(seq),
                id=rec_id,
                description=f"target={target} chain={chain_id} length={len(seq)} file={rel}",
            )
        )
        if i % 200 == 0 or i == len(paths):
            print(f"  parsed {i}/{len(paths)}")

    print(f"skipped={skipped} unmatched={unmatched}")
    return records


def write_fasta(records: List[SeqRecord], path: str) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    n = SeqIO.write(records, path, "fasta")
    print(f"  wrote {n} sequences -> {path}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--outputs-dir",
        default=os.path.join(REPO, "boltzgen_outputs"),
        help="root to recurse for CIF files",
    )
    ap.add_argument(
        "--out-dir",
        default=os.path.join(REPO, "fastas"),
        help="directory for NDM5.fasta and KPC3.fasta",
    )
    ap.add_argument("--targets", nargs="+", default=["NDM5", "KPC3"])
    args = ap.parse_args()

    by_target = extract_all(args.outputs_dir, args.targets)
    for target, recs in by_target.items():
        write_fasta(recs, os.path.join(args.out_dir, f"{target}.fasta"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
