#!/usr/bin/env python3
"""
Stage 2 (BoltzGen) - turn the Stage-1 epitope into a BoltzGen design spec.

IMPORTANT numbering note
------------------------
BoltzGen's YAML `res_index` is 1-based *along the chain as it appears in the input
file* - NOT the author/PDB residue number. Our hotspots are in PDB/scaffold author
numbering (e.g. His120). This script reads the cleaned target structure, builds the
author->ordinal map, and writes `binding:` using ordinals, with the author numbers and
residue names kept as YAML comments so the spec stays auditable.

Produces, per target, under this directory:
  - <T>_design.yaml     BoltzGen design specification (epitope as binding site)
  - <T>_hotspot_map.csv ordinal <-> author(scaffold) <-> residue name

Generation itself needs a GPU + BoltzGen weights (see README.md in this folder); it is
NOT run here. This step only builds correct, ready-to-run inputs.
"""
from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, REPO)
from src.common import structure_io as sio  # noqa: E402

BINDER_LEN = (12, 45)      # peptide-binder length range (residues)
NUM_DESIGNS = 5000         # intermediate candidates (GPU-heavy; tune to budget)
BUDGET = 200               # final ranked designs to keep


def ranges(nums):
    """Compress a sorted int list into 'a..b,c,d..e' notation for YAML readability."""
    nums = sorted(nums)
    out, i = [], 0
    while i < len(nums):
        j = i
        while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
            j += 1
        out.append(f"{nums[i]}..{nums[j]}" if j > i else f"{nums[i]}")
        i = j + 1
    return ",".join(out)


def build(tkey, tdef):
    epi = json.load(open(os.path.join(REPO, "targets", tkey, "epitope", "boltzgen_hotspots.json")))
    chain = epi["scaffold_chain"]
    target_pdb = os.path.join(REPO, "targets", tkey, "structures", f"{tkey}_target.pdb")
    structure = sio.load_structure(target_pdb)
    _, residues = sio.chain_sequence(structure, chain)

    # author resnum -> 1-based ordinal within the chain (BoltzGen res_index frame)
    author_to_ord = {res.id[1]: i + 1 for i, res in enumerate(residues)}
    ord_to_res = {i + 1: (res.id[1], res.get_resname()) for i, res in enumerate(residues)}

    hot_authors = [int(h[len(chain):]) for h in epi["hotspot_residues_scaffold_numbering"]]
    hot_ord = sorted(author_to_ord[a] for a in hot_authors if a in author_to_ord)
    missing = [a for a in hot_authors if a not in author_to_ord]

    # hotspot map CSV
    map_csv = os.path.join(HERE, f"{tkey}_hotspot_map.csv")
    with open(map_csv, "w") as fh:
        fh.write("boltzgen_ordinal,scaffold_resnum,residue\n")
        for o in hot_ord:
            rn, nm = ord_to_res[o]
            fh.write(f"{o},{rn},{nm}\n")

    # comment lines mapping ordinals back to author numbering
    comment = "  # epitope hotspots (BoltzGen ordinal = author/scaffold resnum : residue):\n"
    comment += "\n".join(f"  #   {o} = {ord_to_res[o][0]} {ord_to_res[o][1]}" for o in hot_ord)

    rel_target = os.path.relpath(target_pdb, HERE)
    n_chain = len(residues)
    metals = " (catalytic Zn retained in the target PDB)" if tdef.get("keep_metals") else ""

    yaml = f"""# BoltzGen design spec - {tkey} ({tdef['name']})
# Auto-generated from targets/{tkey}/epitope/boltzgen_hotspots.json
# Target: cleaned design domain{metals}. Binder = de novo peptide.
# Run:  boltzgen run {tkey}_design.yaml --output out_{tkey} --num_designs {NUM_DESIGNS} --budget {BUDGET}
#
# res_index below is 1-based along chain {chain} of {os.path.basename(target_pdb)} (see {tkey}_hotspot_map.csv).
{comment}
entities:
  - file:
      path: {rel_target}
      include:
        - chain:
            id: {chain}
            res_index: 1..{n_chain}
      binding_types:
        - chain:
            id: {chain}
            binding: {ranges(hot_ord)}
  - protein:
      id: P
      sequence: {BINDER_LEN[0]}..{BINDER_LEN[1]}
"""
    out_yaml = os.path.join(HERE, f"{tkey}_design.yaml")
    with open(out_yaml, "w") as fh:
        fh.write(yaml)
    print(f"[{tkey}] {len(hot_ord)} hotspots -> {os.path.relpath(out_yaml, REPO)}"
          + (f"  (WARNING: {missing} not in target chain)" if missing else ""))


def main():
    targets = {k: v for k, v in json.load(open(os.path.join(REPO, "targets", "targets.json"))).items()
               if not k.startswith("_")}
    for tkey, tdef in targets.items():
        build(tkey, tdef)


if __name__ == "__main__":
    main()
