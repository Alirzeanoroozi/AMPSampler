#!/usr/bin/env python3
"""
Cluster protein sequences from a FASTA file at several sequence-identity thresholds.

- Input: CARD protein FASTA (e.g. protein_fasta_protein_homolog_model.fasta)
- Output:
    clusters_50.txt   # clusters at >=50% identity
    clusters_90.txt   # clusters at >=90% identity
    clusters_100.txt  # clusters at 100% identity (identical sequences)
"""

from __future__ import annotations
from typing import Dict, List, Tuple
import argparse
from pathlib import Path


def read_fasta(path: str) -> Dict[str, str]:
    """Read FASTA into {id: sequence}."""
    seqs: Dict[str, str] = {}
    current_id: str | None = None
    parts: List[str] = []

    with open(path, "r") as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                # save previous
                if current_id is not None:
                    seqs[current_id] = "".join(parts)
                # new
                current_id = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.strip())
        # save last
        if current_id is not None:
            seqs[current_id] = "".join(parts)

    return seqs


def seq_identity(a: str, b: str) -> float:
    """
    Simple ungapped sequence identity between two sequences.

    Uses the length of the shorter sequence for normalisation:
        identity = matches / min(len(a), len(b))
    """
    if not a or not b:
        return 0.0
    L = min(len(a), len(b))
    matches = sum(1 for x, y in zip(a[:L], b[:L]) if x == y)
    return matches / L


def greedy_cluster(
    ids: List[str],
    seqs: Dict[str, str],
    cutoff: float,
) -> List[List[str]]:
    """
    Simple greedy clustering:
    - Take the first unassigned sequence as a cluster representative.
    - Any unassigned sequence with identity >= cutoff to the representative
      joins that cluster.
    - Repeat until all sequences are assigned.
    """
    remaining = set(ids)
    clusters: List[List[str]] = []

    while remaining:
        rep = next(iter(remaining))
        rep_seq = seqs[rep]
        cluster = [rep]
        remaining.remove(rep)

        to_remove = []
        for sid in remaining:
            iden = seq_identity(rep_seq, seqs[sid])
            if iden >= cutoff:
                cluster.append(sid)
                to_remove.append(sid)
        for sid in to_remove:
            remaining.remove(sid)

        clusters.append(cluster)

    return clusters


def write_clusters(path: str, clusters: List[List[str]]) -> None:
    """
    Write clusters to a simple text file:
        # Cluster 1 (size N)
        id1
        id2
        ...

        # Cluster 2 (size M)
        ...
    """
    with open(path, "w") as f:
        for i, cl in enumerate(clusters, start=1):
            f.write(f"# Cluster {i} (size {len(cl)})\n")
            for sid in cl:
                f.write(f"{sid}\n")
            f.write("\n")


def main():
    parser = argparse.ArgumentParser(
        description="Cluster CARD protein FASTA by sequence identity."
    )
    parser.add_argument(
        "-i",
        "--input_fasta",
        required=True,
        help="Input FASTA file (e.g. protein_fasta_protein_homolog_model.fasta)",
    )
    parser.add_argument(
        "--out_prefix",
        default="clusters",
        help="Prefix for output files (default: clusters)",
    )
    args = parser.parse_args()

    fasta_path = Path(args.input_fasta)
    if not fasta_path.is_file():
        raise FileNotFoundError(f"Input FASTA not found: {fasta_path}")

    print(f"Reading FASTA: {fasta_path}")
    seqs = read_fasta(str(fasta_path))
    ids = list(seqs.keys())
    print(f"Loaded {len(ids)} sequences.")

    # Thresholds: 50%, 90%, 100%
    thresholds: List[Tuple[float, str]] = [
        (0.50, "50"),
        (0.90, "90"),
        (1.00, "100"),
    ]

    for thr, label in thresholds:
        print(f"\nClustering at >= {int(float(label))}% identity (cutoff={thr:.2f})")
        clusters = greedy_cluster(ids, seqs, cutoff=thr)
        out_path = f"{args.out_prefix}_{label}.txt"
        write_clusters(out_path, clusters)
        print(f"  -> wrote {len(clusters)} clusters to {out_path}")


if __name__ == "__main__":
    main()