#!/usr/bin/env python3
"""
Stage 0 - Target preparation.

For each target:
  1. Align the precursor sequence to the design-scaffold crystal structure to find
     the folded (crystallographically-resolved) catalytic domain.
  2. Emit the design-domain FASTA (what sequence-based tools index) and a
     precursor<->scaffold numbering map.
  3. Emit a cleaned design-target PDB: scaffold protein chain + catalytic metals
     (Zn for NDM), with the bound inhibitor and waters/buffer removed. This apo+metal
     structure is the actual target handed to Stage 2 generation.
  4. Annotate the predicted signal peptide (literature/UniProt convention) so the
     biological cleavage is documented separately from the crystallographic start.

Mutations vs template (NDM-5: V88L/M154L; KPC-3: H274Y) are already present in the
target *sequence*; they are conservative and distal to the epitope, so the template
backbone is used as the design scaffold. apply_mutations.py (stub) documents how to
introduce them into the *structure* with PyMOL/PDBFixer when that becomes available.

Runs in the base conda env (Biopython + numpy).
"""
from __future__ import annotations

import json
import os
import sys

from Bio import SeqIO
from Bio.PDB import PDBIO, Select

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, REPO)
from src.common import structure_io as sio  # noqa: E402

CACHE = os.path.join(REPO, "targets", "_pdb_cache")

# Predicted signal-peptide length (UniProt/SignalP convention). Documented, not used
# to define the design domain (the crystallographic domain is used for that).
PREDICTED_SIGNAL_PEPTIDE = {
    "NDM5": {"end": 28, "note": "NDM-1 lipoprotein signal peptide ~1-28 (lipidation near Cys26); "
                                 "soluble crystallized constructs typically start ~residue 36-43."},
    "KPC3": {"end": 24, "note": "KPC-2 Sec signal peptide ~1-24; mature enzyme 25-293."},
}


class _TargetSelect(Select):
    def __init__(self, chain_id, keep_metals):
        self.chain_id = chain_id
        self.keep_metals = {m.upper() for m in keep_metals}

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, res):
        name = res.get_resname().strip().upper()
        if res.id[0] == " " and name in sio.THREE_TO_ONE:
            return True
        return name in self.keep_metals  # metals kept; inhibitor/waters dropped

    def accept_atom(self, atom):
        return (not atom.is_disordered()) or atom.get_altloc() in ("", "A")


def run_target(tkey, tdef):
    print(f"\n=== Stage 0: {tkey} ({tdef['name']}) ===")
    precursor = str(next(SeqIO.parse(os.path.join(REPO, tdef["precursor_fasta"]), "fasta")).seq).upper()
    scaffold_pdb = tdef["design_scaffold_pdb"]
    chain = tdef.get("scaffold_chain", "A")
    keep_metals = tdef.get("keep_metals", [])
    out_struct = os.path.join(REPO, "targets", tkey, "structures")
    os.makedirs(out_struct, exist_ok=True)

    path = sio.download_pdb(scaffold_pdb, CACHE)
    structure = sio.load_structure(path)
    smap = sio.map_precursor_to_chain(precursor, chain, structure)
    mature_start = smap["mature_start_prec"]
    mature_end = max(smap["prec_to_resnum"].keys())
    design_seq = precursor[mature_start - 1: mature_end]

    # --- design-domain FASTA ---
    fasta = os.path.join(REPO, "targets", tkey, "design_domain.fasta")
    with open(fasta, "w") as fh:
        fh.write(f">{tkey}_design_domain precursor {mature_start}-{mature_end} "
                 f"(folded domain resolved in {scaffold_pdb}); mutations {','.join(tdef.get('mutations_vs_template', []))}\n")
        for i in range(0, len(design_seq), 60):
            fh.write(design_seq[i:i + 60] + "\n")

    # --- numbering map CSV ---
    nmap = os.path.join(REPO, "targets", tkey, "numbering_map.csv")
    with open(nmap, "w") as fh:
        fh.write("precursor_pos,residue,scaffold_resnum\n")
        for ppos in sorted(smap["prec_to_resnum"]):
            fh.write(f"{ppos},{precursor[ppos-1]},{smap['prec_to_resnum'][ppos]}\n")

    # --- cleaned design-target structure (protein + metals, no inhibitor/water) ---
    io = PDBIO()
    io.set_structure(structure)
    target_pdb = os.path.join(out_struct, f"{tkey}_target.pdb")
    io.save(target_pdb, _TargetSelect(chain, keep_metals))
    n_metal = sum(1 for m in sio.detect_metals(structure, keep_metals)
                  if m.chain == chain) if keep_metals else 0

    # --- prep report ---
    sp = PREDICTED_SIGNAL_PEPTIDE.get(tkey, {})
    report = os.path.join(REPO, "targets", tkey, "prep_report.md")
    with open(report, "w") as fh:
        fh.write(
            f"# Stage 0 target prep - {tkey} ({tdef['name']})\n\n"
            f"- Class: {tdef['class']}\n"
            f"- Precursor length: {len(precursor)} aa\n"
            f"- Predicted signal peptide (literature): residues 1-{sp.get('end','?')}. {sp.get('note','')}\n"
            f"- Design scaffold: {scaffold_pdb} chain {chain} "
            f"(seq identity to precursor {round(smap['identity'],3)})\n"
            f"- Folded/design domain (crystallographically resolved): precursor "
            f"{mature_start}-{mature_end} ({len(design_seq)} aa)\n"
            f"- Catalytic metals retained in target structure: {keep_metals or 'none'} "
            f"({n_metal} ion(s))\n"
            f"- Mutations vs template (already in target sequence; apply to structure with "
            f"apply_mutations.py if needed): {', '.join(tdef.get('mutations_vs_template', [])) or 'none'}\n\n"
            f"## Artifacts\n"
            f"- `design_domain.fasta` - design-domain sequence\n"
            f"- `numbering_map.csv` - precursor <-> scaffold residue map\n"
            f"- `structures/{tkey}_target.pdb` - cleaned design target (protein + metals, inhibitor/waters removed)\n"
        )
    print(f"  design domain: precursor {mature_start}-{mature_end} ({len(design_seq)} aa); "
          f"metals kept: {keep_metals or 'none'} ({n_metal})")
    print(f"  wrote {os.path.relpath(fasta, REPO)}, {os.path.relpath(target_pdb, REPO)}")


def main():
    with open(os.path.join(REPO, "targets", "targets.json")) as fh:
        targets = {k: v for k, v in json.load(fh).items() if not k.startswith("_")}
    for tkey, tdef in targets.items():
        run_target(tkey, tdef)


if __name__ == "__main__":
    main()
