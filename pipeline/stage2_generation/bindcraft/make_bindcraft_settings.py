#!/usr/bin/env python3
"""
Stage 2 (BindCraft) - emit a BindCraft target-settings JSON from the Stage-1 epitope.

BindCraft uses author/PDB residue numbering for `target_hotspot_residues` (unlike
BoltzGen), so we pass the scaffold author numbers directly. Generation needs a GPU +
the BindCraft install (see ../README.md); not run here.
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))

BINDER_LEN = [12, 45]
N_FINAL = 100


def build(tkey, tdef):
    epi = json.load(open(os.path.join(REPO, "targets", tkey, "epitope", "boltzgen_hotspots.json")))
    chain = epi["scaffold_chain"]
    authors = sorted(int(h[len(chain):]) for h in epi["hotspot_residues_scaffold_numbering"])
    target_pdb = os.path.relpath(os.path.join(REPO, "targets", tkey, "structures", f"{tkey}_target.pdb"), HERE)
    settings = {
        "design_path": f"./out_{tkey}/",
        "binder_name": f"{tkey}_binder",
        "starting_pdb": target_pdb,
        "chains": chain,
        "target_hotspot_residues": ",".join(str(a) for a in authors),
        "lengths": BINDER_LEN,
        "number_of_final_designs": N_FINAL,
        "_note": f"hotspots are the {tkey} active-site epitope (PDB author numbering); "
                 f"catalytic metals {'retained' if tdef.get('keep_metals') else 'n/a'}.",
    }
    out = os.path.join(HERE, f"{tkey}_bindcraft.json")
    json.dump(settings, open(out, "w"), indent=2)
    print(f"[{tkey}] {len(authors)} hotspots -> {os.path.relpath(out, REPO)}")


def main():
    targets = {k: v for k, v in json.load(open(os.path.join(REPO, "targets", "targets.json"))).items()
               if not k.startswith("_")}
    for tkey, tdef in targets.items():
        build(tkey, tdef)


if __name__ == "__main__":
    main()
