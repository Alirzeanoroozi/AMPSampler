#!/usr/bin/env python3
"""
Stage 3 - Active-site overlap scorer  (repurposed from the old hotspot.py).

The old hotspot.py just *reported* which target residues a design touched. For
inhibitor design the question is sharper: does the binder sit ON the catalytic
epitope? This scores, for each binder-target complex, how well the binder's
interface footprint overlaps the Stage-1 active-site epitope.

Metrics per complex:
  epitope_recall      |interface ∩ epitope| / |epitope|     (coverage of active site)
  interface_precision |interface ∩ epitope| / |interface|   (focus on active site vs elsewhere)
  n_catalytic_contacts number of catalytic / metal-ligand residues contacted
  catalytic_ok         True if it touches the catalytic core (good inhibitor prior)

Numbering is resolved by aligning each complex's target chain to the design domain,
so it works regardless of how the generator numbered residues (BoltzGen re-indexes
from 1; BindCraft/RFdiffusion keep author numbers).

Runs in the base conda env. Usage:
  python active_site_overlap.py --target NDM5 --complexes <dir-or-file> [--cutoff 4.5] [--out out.csv]
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import re
import sys

# BoltzGen prefixes its final-ranked complex files with "rank####_" (e.g.
# rank0001_NDM5_design_1726.cif). The canonical design_id used elsewhere in the
# pipeline (FASTA record id == BoltzGen metrics CSV `id`) does NOT have this
# prefix, so we strip it here to keep design_id consistent across all score CSVs
# and let build_manifest join them on a single key.
_RANK_PREFIX_RE = re.compile(r"^rank\d+_", re.IGNORECASE)


def _design_id_from_path(path: str) -> str:
    stem = os.path.splitext(os.path.basename(path))[0]
    return _RANK_PREFIX_RE.sub("", stem)

import numpy as np
from Bio import SeqIO

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, REPO)
from src.common import structure_io as sio  # noqa: E402


def load_epitope(tkey):
    """Return (set of precursor positions, set of catalytic/metal-ligand precursor positions)."""
    path = os.path.join(REPO, "targets", tkey, "epitope", "hotspots_union.csv")
    epi, catalytic = set(), set()
    for r in csv.DictReader(open(path)):
        p = int(r["precursor_pos"])
        epi.add(p)
        if r["metal_ligand"].strip().lower() == "true":
            catalytic.add(p)
    return epi, catalytic


def load_catalytic_extra(tkey):
    """For KPC (no metals), treat the catalytic Ser/Lys/Glu as the catalytic core."""
    path = os.path.join(REPO, "targets", tkey, "epitope", "hotspots_union.csv")
    core = set()
    rows = list(csv.DictReader(open(path)))
    if any(r["metal_ligand"].strip().lower() == "true" for r in rows):
        return core  # metal enzyme: handled by metal_ligand flag
    # serine enzyme: the closest-contacting Ser plus any Lys/Glu in the epitope
    sers = [r for r in rows if r["resname"] == "SER"]
    if sers:
        nucleophile = min(sers, key=lambda r: float(r["min_dist_A"]))
        core.add(int(nucleophile["precursor_pos"]))
    for r in rows:
        if r["resname"] in ("LYS", "GLU"):
            core.add(int(r["precursor_pos"]))
    return core


def identify_chains(structure, design_seq):
    """Pick target chain = best sequence match to the design domain; binders = the rest."""
    model = next(structure.get_models())
    target_chain, best_id = None, -1.0
    chain_seqs = {}
    for chain in model:
        seq, residues = sio.chain_sequence(structure, chain.id)
        if len(seq) < 5:
            continue
        chain_seqs[chain.id] = (seq, residues)
        # crude identity via alignment score normalised by design length
        from src.common.structure_io import _aligner
        score = _aligner().score(design_seq, seq)
        ident = score / (2.0 * min(len(design_seq), len(seq)))  # match_score=2
        if ident > best_id:
            best_id, target_chain = ident, chain.id
    binder_chains = [c for c in chain_seqs if c != target_chain]
    return target_chain, binder_chains, chain_seqs


def score_complex(path, tkey, precursor, design_seq, epitope, catalytic, cutoff):
    structure = sio.load_structure(path)
    target_chain, binder_chains, chain_seqs = identify_chains(structure, design_seq)
    if not binder_chains:
        return {"design_id": os.path.basename(path), "status": "no binder chain found"}

    # map target chain residues -> precursor positions
    mapping = sio.map_precursor_to_chain(precursor, target_chain, structure)
    r2p = mapping["resnum_to_prec"]

    # binder heavy-atom coords
    model = next(structure.get_models())
    binder_coords = []
    for bc in binder_chains:
        for res in model[bc]:
            for atom in res.get_atoms():
                if atom.element != "H":
                    binder_coords.append(atom.coord)
    binder_coords = np.array(binder_coords, dtype=float)

    # target residues within cutoff of binder = interface footprint (precursor positions)
    contacts = sio.residues_near(structure, target_chain, binder_coords, cutoff)
    interface = {r2p[rn] for rn in contacts if rn in r2p}

    inter_epi = interface & epitope
    cat_hit = interface & catalytic
    binder_len = len(chain_seqs[binder_chains[0]][0]) if binder_chains else 0
    return {
        "design_id": _design_id_from_path(path),
        "binder_len": binder_len,
        "n_interface": len(interface),
        "n_epitope_hit": len(inter_epi),
        "epitope_recall": round(len(inter_epi) / len(epitope), 3) if epitope else 0.0,
        "interface_precision": round(len(inter_epi) / len(interface), 3) if interface else 0.0,
        "n_catalytic_contacts": len(cat_hit),
        "catalytic_ok": len(cat_hit) > 0,
        "contacted_epitope": "|".join(map(str, sorted(inter_epi))),
        "target_chain": target_chain,
        "status": "ok",
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True, choices=["NDM5", "KPC3"])
    ap.add_argument("--complexes", required=True, help="binder-target complex file or directory (.pdb/.cif)")
    ap.add_argument("--cutoff", type=float, default=4.5)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    tkey = args.target
    precursor = str(next(SeqIO.parse(os.path.join(REPO, "targets", tkey, "precursor.fasta"), "fasta")).seq).upper()
    design_seq = str(next(SeqIO.parse(os.path.join(REPO, "targets", tkey, "design_domain.fasta"), "fasta")).seq).upper()
    epitope, catalytic = load_epitope(tkey)
    catalytic |= load_catalytic_extra(tkey)

    if os.path.isdir(args.complexes):
        files = sorted(glob.glob(os.path.join(args.complexes, "*.pdb")) +
                       glob.glob(os.path.join(args.complexes, "*.cif")))
    else:
        files = [args.complexes]
    if not files:
        print("No complex files found.")
        return 1

    rows = []
    for f in files:
        try:
            rows.append(score_complex(f, tkey, precursor, design_seq, epitope, catalytic, args.cutoff))
        except Exception as e:
            rows.append({"design_id": _design_id_from_path(f), "status": f"ERROR {type(e).__name__}: {e}"})

    out = args.out or os.path.join(REPO, "results", "stage3_validation", f"{tkey}_active_site_overlap.csv")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    cols = ["design_id", "binder_len", "n_interface", "n_epitope_hit", "epitope_recall",
            "interface_precision", "n_catalytic_contacts", "catalytic_ok",
            "contacted_epitope", "target_chain", "status"]
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in cols})
    ok = [r for r in rows if r.get("status") == "ok"]
    print(f"[{tkey}] scored {len(ok)}/{len(rows)} complexes (epitope={len(epitope)} residues, "
          f"catalytic core={len(catalytic)}). Cutoff {args.cutoff} A.")
    if ok:
        print(f"  catalytic_ok: {sum(r['catalytic_ok'] for r in ok)}/{len(ok)} ; "
              f"median epitope_recall {np.median([r['epitope_recall'] for r in ok]):.2f}")
    print(f"  -> {os.path.relpath(out, REPO)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
