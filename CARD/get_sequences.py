#!/usr/bin/env python3
"""
Extract representative sequences from 50% identity clusters.

- Input:
    - card_clusters_50.txt   (output from cluster.py)
    - protein_fasta_protein_homolog_model.fasta  (CARD protein FASTA)
- Output:
    - card_clusters_50_reps.fasta  (one representative per cluster)
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, List
import argparse


def read_fasta(path: str) -> Dict[str, str]:
    """Read FASTA into {id: sequence} where id is the first token after '>'."""
    seqs: Dict[str, str] = {}
    current_id = None
    parts: List[str] = []

    with open(path, "r") as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    seqs[current_id] = "".join(parts)
                current_id = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.strip())
        if current_id is not None:
            seqs[current_id] = "".join(parts)

    return seqs


def read_cluster_representatives(path: str) -> List[str]:
    """
    Read cluster file produced by cluster.py and return the representative
    sequence ID for each cluster (the first ID after the '# Cluster ...' line).
    """
    reps: List[str] = []
    with open(path, "r") as f:
        current_rep = None
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("# Cluster "):
                # Starting a new cluster; reset
                current_rep = None
                continue
            # First non-empty, non-comment line after a cluster header is the rep
            if current_rep is None:
                current_rep = line
                reps.append(current_rep)

    return reps


def write_fasta(path: str, ids: List[str], seqs: Dict[str, str]) -> None:
    """Write given IDs and sequences to a FASTA file."""
    with open(path, "w") as f:
        for sid in ids:
            seq = seqs.get(sid)
            if seq is None:
                # skip if not found in FASTA
                continue
            f.write(f">{sid}\n")
            # wrap to 60 chars per line for readability
            for i in range(0, len(seq), 60):
                f.write(seq[i : i + 60] + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Extract representative sequences for 50% clusters."
    )
    parser.add_argument(
        "--clusters",
        default="card_clusters_50.txt",
        help="Cluster file (default: card_clusters_50.txt)",
    )
    parser.add_argument(
        "--fasta",
        default="card-data/protein_fasta_protein_homolog_model.fasta",
        help="Original CARD protein FASTA (default: card-data/protein_fasta_protein_homolog_model.fasta)",
    )
    parser.add_argument(
        "--out",
        default="card_clusters_50_reps.fasta",
        help="Output FASTA with representatives (default: card_clusters_50_reps.fasta)",
    )
    args = parser.parse_args()

    clusters_path = Path(args.clusters)
    fasta_path = Path(args.fasta)

    if not clusters_path.is_file():
        raise FileNotFoundError(f"Cluster file not found: {clusters_path}")
    if not fasta_path.is_file():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    print(f"Reading FASTA:   {fasta_path}")
    seqs = read_fasta(str(fasta_path))
    print(f"Loaded {len(seqs)} sequences.")

    print(f"Reading clusters: {clusters_path}")
    reps = read_cluster_representatives(str(clusters_path))
    print(f"Found {len(reps)} cluster representatives.")

    print(f"Writing representatives to: {args.out}")
    write_fasta(args.out, reps, seqs)
    print("Done.")


if __name__ == "__main__":
    main()