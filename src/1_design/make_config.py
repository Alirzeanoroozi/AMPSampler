#!/usr/bin/env python3
"""
Stage 2 (BoltzGen) - turn the Stage-1 epitope into a BoltzGen design spec.

Sequence-based (no target PDB): the design-domain FASTA is the fixed target chain,
and epitope hotspots are 1-based positions along that sequence (`mature_pos` from
Stage 1). BoltzGen will fold the target together with the designed binder.

Produces, per target, under boltzgen_inputs/ (repo root):
  - <T>_design.yaml     BoltzGen design specification (epitope as binding site)
  - <T>_hotspot_map.csv sequence_pos <-> scaffold_resnum <-> residue name

Generation itself needs a GPU + BoltzGen weights; it is NOT run here.
"""
from __future__ import annotations

import csv
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
OUT_DIR = os.path.join(REPO, "boltzgen_inputs")

BINDER_LEN = (12, 45)      # peptide-binder length range (residues)
NUM_DESIGNS = 5000         # intermediate candidates (GPU-heavy; tune to budget)
BUDGET = 200               # final ranked designs to keep


def read_seq(path: str) -> str:
    lines = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            lines.append(line.strip())
    return "".join(lines).upper()


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


def load_hotspots(tkey: str):
    """Return list of (seq_pos, scaffold_resnum, resname) from Stage-1 union CSV."""
    union_csv = os.path.join(REPO, "targets", tkey, "epitope", "hotspots_union.csv")
    rows = []
    with open(union_csv) as fh:
        for r in csv.DictReader(fh):
            rows.append((int(r["mature_pos"]), int(r["scaffold_resnum"]), r["resname"]))
    return sorted(rows, key=lambda x: x[0])


def build(tkey, tdef):
    seq_path = os.path.join(REPO, "targets", tkey, "design_domain.fasta")
    seq = read_seq(seq_path)
    hotspots = load_hotspots(tkey)
    hot_ord = [pos for pos, _, _ in hotspots]
    missing = [pos for pos in hot_ord if pos < 1 or pos > len(seq)]
    if missing:
        raise ValueError(f"{tkey}: hotspot positions outside sequence: {missing}")

    os.makedirs(OUT_DIR, exist_ok=True)

    map_csv = os.path.join(OUT_DIR, f"{tkey}_hotspot_map.csv")
    with open(map_csv, "w") as fh:
        fh.write("sequence_pos,scaffold_resnum,residue\n")
        for pos, scaf, name in hotspots:
            fh.write(f"{pos},{scaf},{name}\n")

    comment = "  # epitope hotspots (sequence pos = scaffold_resnum : residue):\n"
    comment += "\n".join(
        f"  #   {pos} = {scaf} {name}" for pos, scaf, name in hotspots
    )

    yaml = f"""# BoltzGen design spec - {tkey} ({tdef['name']})
# Auto-generated from targets/{tkey}/design_domain.fasta + epitope/hotspots_union.csv
# Target = fixed sequence (no structure file). Binder = de novo peptide.
# Run:  boltzgen run {tkey}_design.yaml --output out_{tkey} --num_designs {NUM_DESIGNS} --budget {BUDGET}
#
# binding indices are 1-based along the design-domain sequence (see {tkey}_hotspot_map.csv).
{comment}
entities:
  - protein:
      id: A
      sequence: {seq}
      binding_types:
        binding: {ranges(hot_ord)}
  - protein:
      id: P
      sequence: {BINDER_LEN[0]}..{BINDER_LEN[1]}
"""
    out_yaml = os.path.join(OUT_DIR, f"{tkey}_design.yaml")
    with open(out_yaml, "w") as fh:
        fh.write(yaml)
    print(f"[{tkey}] {len(hot_ord)} hotspots / {len(seq)} aa -> {os.path.relpath(out_yaml, REPO)}")


def main():
    targets = {
        k: v
        for k, v in json.load(open(os.path.join(REPO, "targets", "targets.json"))).items()
        if not k.startswith("_")
    }
    for tkey, tdef in targets.items():
        build(tkey, tdef)


if __name__ == "__main__":
    main()
