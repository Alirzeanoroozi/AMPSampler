#!/usr/bin/env python3
"""
Stage 5 - Generate alanine-scan specificity controls (mirrors the BoltzGen GyrA control).

For each selected design we make negative controls by mutating the binder residues that
form the designed interface to alanine. If binding/inhibition drops for the mutant but not
the parent, the effect is via the designed epitope interface - the key specificity check.

Two modes per design:
  * complex mode (preferred): a binder-target complex PDB/CIF is available -> mutate the
    binder residues within `cutoff` of the target to Ala (one combined 'interface-Ala'
    control + per-interface-residue single-Ala mutants).
  * sequence-only fallback: full single-alanine scan (every non-Ala position -> Ala).

Usage:
  python alanine_scan.py --target NDM5 --selected results/selected_NDM5.fasta [--complexes <dir>] [--cutoff 4.5]
Outputs: results/alanine_controls_<T>.fasta  +  results/alanine_controls_<T>.csv
"""
from __future__ import annotations

import argparse, csv, glob, os, sys

import numpy as np
from Bio import SeqIO

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, REPO)
from src.common import structure_io as sio  # noqa: E402


def find_complex(design_id, comp_dir):
    if not comp_dir:
        return None
    for ext in ("pdb", "cif"):
        hits = glob.glob(os.path.join(comp_dir, f"*{design_id}*.{ext}"))
        if hits:
            return hits[0]
    return None


def binder_interface_positions(path, cutoff):
    """Return (binder_seq, sorted 0-based interface positions) using length-based chain ID."""
    structure = sio.load_structure(path)
    model = next(structure.get_models())
    chains = {}
    for ch in model:
        seq, residues = sio.chain_sequence(structure, ch.id)
        if seq:
            chains[ch.id] = (seq, residues)
    if len(chains) < 2:
        return None, []
    # target = longest chain; binder = shortest
    target = max(chains, key=lambda c: len(chains[c][0]))
    binder = min(chains, key=lambda c: len(chains[c][0]))
    tgt_coords = np.array([a.coord for r in chains[target][1] for a in r.get_atoms() if a.element != "H"])
    bseq, bres = chains[binder]
    iface = []
    for idx, res in enumerate(bres):
        atoms = np.array([a.coord for a in res.get_atoms() if a.element != "H"])
        d = np.sqrt(((atoms[:, None, :] - tgt_coords[None, :, :]) ** 2).sum(-1))
        if d.min() <= cutoff:
            iface.append(idx)
    return bseq, iface


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True)
    ap.add_argument("--selected", required=True)
    ap.add_argument("--complexes", default=None)
    ap.add_argument("--cutoff", type=float, default=4.5)
    args = ap.parse_args()

    out_fa = os.path.join(REPO, "results", f"alanine_controls_{args.target}.fasta")
    out_csv = os.path.join(REPO, "results", f"alanine_controls_{args.target}.csv")
    os.makedirs(os.path.dirname(out_fa), exist_ok=True)

    fa, rows = open(out_fa, "w"), []
    n_designs = 0
    for rec in SeqIO.parse(args.selected, "fasta"):
        did, parent = rec.id, str(rec.seq).upper()
        n_designs += 1
        comp = find_complex(did, args.complexes)
        mode = "complex" if comp else "sequence-only"
        positions = []
        if comp:
            bseq, iface = binder_interface_positions(comp, args.cutoff)
            if bseq and iface:
                positions = iface
                # one combined interface-Ala control
                mut = list(parent)
                for i in iface:
                    if i < len(mut):
                        mut[i] = "A"
                fa.write(f">{did}_ifaceAla mode=complex iface={','.join(str(p+1) for p in iface)}\n{''.join(mut)}\n")
                rows.append({"design_id": did, "control": f"{did}_ifaceAla", "mode": mode,
                             "mutated_positions": ",".join(str(p + 1) for p in iface)})
        if not positions:  # sequence-only fallback: full single-Ala scan
            positions = [i for i, a in enumerate(parent) if a != "A"]
        # per-position single-Ala mutants (interface residues, or all positions in fallback)
        for i in positions:
            if parent[i] == "A":
                continue
            mut = parent[:i] + "A" + parent[i + 1:]
            name = f"{did}_{parent[i]}{i+1}A"
            fa.write(f">{name} mode={mode} parent={did}\n{mut}\n")
            rows.append({"design_id": did, "control": name, "mode": mode,
                         "mutated_positions": str(i + 1)})
    fa.close()
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["design_id", "control", "mode", "mutated_positions"])
        w.writeheader(); w.writerows(rows)
    print(f"[{args.target}] {n_designs} designs -> {len(rows)} alanine controls "
          f"({'complex+single' if args.complexes else 'single-Ala scan (no complexes given)'}) "
          f"-> results/alanine_controls_{args.target}.fasta")


if __name__ == "__main__":
    main()
