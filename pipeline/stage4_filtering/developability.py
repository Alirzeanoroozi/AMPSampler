#!/usr/bin/env python3
"""
Stage 4 - Developability, synthesizability, and delivery-proxy descriptors.

Reframing vs the old AMP pipeline:
  - "Is it an AMP?" / "similarity to an AMP database" are DROPPED as objectives. The
    designs are folded enzyme-inhibitor binders, not membrane lytics; AMP-likeness is
    not the goal and the old similarity-to-DBAASP ranking actively penalised novelty.
  - The cationic/amphipathic ("AMP-like") signals are KEPT but RELABELLED as a
    *periplasmic-delivery proxy*: NDM-5 and KPC-3 are periplasmic, so a binder must
    cross the outer membrane to reach them. Charge + amphipathicity correlate with OM
    permeation, so they inform delivery, not efficacy.
  - Real developability gates are ADDED: synthesis liabilities, aggregation proxy.

Self-contained (no scipy / modlamp). Usage:
  python developability.py --fasta designs.fasta --out developability_<T>.csv
  python developability.py --manifest results/manifest_<T>.csv --out developability_<T>.csv
"""
from __future__ import annotations

import argparse, csv, math, os, re
from collections import Counter

# Kyte-Doolittle hydropathy
KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
      "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
      "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}
# Eisenberg consensus (for hydrophobic moment / amphipathicity)
EIS = {"A": 0.62, "R": -2.53, "N": -0.78, "D": -0.90, "C": 0.29, "Q": -0.85, "E": -0.74,
       "G": 0.48, "H": -0.40, "I": 1.38, "L": 1.06, "K": -1.50, "M": 0.64, "F": 1.19,
       "P": 0.12, "S": -0.18, "T": -0.05, "W": 0.81, "Y": 0.26, "V": 1.08}
PKA_POS = {"Nterm": 9.0, "K": 10.5, "R": 12.5, "H": 6.0}
PKA_NEG = {"Cterm": 2.0, "D": 3.9, "E": 4.1, "C": 8.3, "Y": 10.1}


def net_charge(seq, ph=7.4):
    c = Counter(seq)
    pos = 1.0 / (1 + 10 ** (ph - PKA_POS["Nterm"]))
    neg = 1.0 / (1 + 10 ** (PKA_NEG["Cterm"] - ph))
    for aa, pk in (("K", PKA_POS["K"]), ("R", PKA_POS["R"]), ("H", PKA_POS["H"])):
        pos += c.get(aa, 0) * 1.0 / (1 + 10 ** (ph - pk))
    for aa, pk in (("D", PKA_NEG["D"]), ("E", PKA_NEG["E"]), ("C", PKA_NEG["C"]), ("Y", PKA_NEG["Y"])):
        neg += c.get(aa, 0) * 1.0 / (1 + 10 ** (pk - ph))
    return pos - neg


def gravy(seq):
    vals = [KD[a] for a in seq if a in KD]
    return sum(vals) / len(vals) if vals else 0.0


def hydrophobic_moment(seq, window=11, angle=100):
    v = [EIS.get(a, 0.0) for a in seq]
    if len(v) < 2:
        return 0.0
    w = min(window, len(v))
    best = 0.0
    for i in range(len(v) - w + 1):
        seg = v[i:i + w]
        sc = sum(h * math.cos(math.radians(angle * k)) for k, h in enumerate(seg))
        ss = sum(h * math.sin(math.radians(angle * k)) for k, h in enumerate(seg))
        best = max(best, math.sqrt(sc * sc + ss * ss) / w)
    return best


def max_hydrophobic_window(seq, window=7):
    """Aggregation proxy: most hydrophobic KD window mean (higher = more aggregation-prone)."""
    if len(seq) < window:
        return gravy(seq)
    return max(gravy(seq[i:i + window]) for i in range(len(seq) - window + 1))


def liabilities(seq):
    flags = []
    cys = seq.count("C")
    if cys % 2 == 1:
        flags.append("odd_Cys(free_thiol)")
    if cys > 2:
        flags.append(f"{cys}_Cys")
    if seq[:1] in ("Q", "E"):
        flags.append("N-term_Q/E(pyroglutamate)")
    if seq[:1] == "M":
        flags.append("N-term_Met(ox)")
    for motif in ("NG", "DG", "DP", "DS", "NS"):  # deamidation / isomerisation / cleavage
        if motif in seq:
            flags.append(f"motif_{motif}")
    if re.search(r"(.)\1{3,}", seq):
        flags.append("homopolymer_run>=4")
    if seq.count("M") + seq.count("W") >= 3:
        flags.append("oxidation_prone(M/W)")
    return flags


def descriptors(seq):
    seq = "".join(c for c in seq.upper() if c.isalpha())
    q = net_charge(seq)
    mu = hydrophobic_moment(seq)
    flags = liabilities(seq)
    # delivery proxy: periplasmic access favoured by net + charge AND amphipathicity
    delivery = round(min(max(q, 0) / 6.0, 1.0) * 0.5 + min(mu / 0.6, 1.0) * 0.5, 3)
    return {
        "length": len(seq),
        "net_charge_pH7.4": round(q, 2),
        "gravy": round(gravy(seq), 3),
        "hydrophobic_moment": round(mu, 3),
        "aggregation_proxy": round(max_hydrophobic_window(seq), 3),
        "cys_count": seq.count("C"),
        "delivery_proxy": delivery,
        "n_liabilities": len(flags),
        "liabilities": ";".join(flags),
    }


def iter_seqs(args):
    if args.fasta:
        from Bio import SeqIO
        for rec in SeqIO.parse(args.fasta, "fasta"):
            yield rec.id, str(rec.seq)
    else:
        for r in csv.DictReader(open(args.manifest)):
            if r.get("sequence"):
                yield r["design_id"], r["sequence"]


def main():
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--fasta")
    g.add_argument("--manifest")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cols = ["design_id", "length", "net_charge_pH7.4", "gravy", "hydrophobic_moment",
            "aggregation_proxy", "cys_count", "delivery_proxy", "n_liabilities", "liabilities"]
    n = 0
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for did, seq in iter_seqs(args):
            d = descriptors(seq)
            d["design_id"] = did
            w.writerow(d)
            n += 1
    print(f"Computed developability for {n} designs -> {args.out}")


if __name__ == "__main__":
    main()
