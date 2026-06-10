#!/usr/bin/env python3
"""
Generate a FASTA file with random protein sequences.

- Default: 1000 sequences
- Length: random between 10 and 50 (inclusive)
"""

import random
import argparse
from pathlib import Path

# Standard 20 amino acids
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


def random_protein(min_len: int, max_len: int) -> str:
    L = random.randint(min_len, max_len)
    return "".join(random.choice(AMINO_ACIDS) for _ in range(L))


def write_random_fasta(out_path: str, n_seqs: int, min_len: int, max_len: int) -> None:
    with open(out_path, "w") as f:
        for i in range(1, n_seqs + 1):
            seq = random_protein(min_len, max_len)
            f.write(f">rand_{i}\n")
            # wrap at 60 chars per line (just style)
            for j in range(0, len(seq), 60):
                f.write(seq[j:j+60] + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Generate random protein sequences in FASTA format."
    )
    parser.add_argument(
        "-o", "--out",
        default="random_proteins.fasta",
        help="Output FASTA file (default: random_proteins.fasta)",
    )
    parser.add_argument(
        "-n", "--num",
        type=int,
        default=1000,
        help="Number of sequences to generate (default: 1000)",
    )
    parser.add_argument(
        "--min_len",
        type=int,
        default=10,
        help="Minimum sequence length (default: 10)",
    )
    parser.add_argument(
        "--max_len",
        type=int,
        default=50,
        help="Maximum sequence length (default: 50)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed (optional, for reproducibility)",
    )
    args = parser.parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    if args.min_len <= 0 or args.max_len < args.min_len:
        raise ValueError("Lengths must satisfy: 0 < min_len <= max_len")

    out_path = Path(args.out)
    write_random_fasta(str(out_path), args.num, args.min_len, args.max_len)
    print(f"Wrote {args.num} sequences to {out_path} "
          f"(length {args.min_len}–{args.max_len})")


if __name__ == "__main__":
    main()