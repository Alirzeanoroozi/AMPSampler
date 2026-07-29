#!/usr/bin/env python3
"""
Apply the NDM-5 / KPC-3 point mutations to the design-target *structure*.

The target SEQUENCES already carry these substitutions (the README sequences are
NDM-5 / KPC-3). The design *scaffold structures* are the parent enzymes
(NDM-1 = 5YPM, KPC-2 = 3DW0), so to make the structures exactly NDM-5 / KPC-3 you
apply, in scaffold numbering:

  NDM5:  V88L, M154L
  KPC3:  H274Y

All three are conservative and >15 A from the epitope/active site, so for binder
design against the active-site epitope the parent backbone is an adequate scaffold
and this step is optional. Proper side-chain rebuild needs a rotamer tool that is
NOT installed in this environment (no PyMOL / PDBFixer / Rosetta). This script
performs the mutation if PyMOL is importable, otherwise prints the exact commands.
"""
import os
import sys

MUTATIONS = {
    "NDM5": [("A", 88, "LEU"), ("A", 154, "LEU")],
    "KPC3": [("A", 274, "TYR")],
}
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def with_pymol(target):
    import pymol
    from pymol import cmd
    pymol.finish_launching(["pymol", "-qc"])
    src = os.path.join(REPO, "targets", target, "structures", f"{target}_target.pdb")
    cmd.load(src, "t")
    cmd.wizard("mutagenesis")
    for chain, resi, new in MUTATIONS[target]:
        cmd.get_wizard().set_mode(new)
        cmd.get_wizard().do_select(f"t and chain {chain} and resi {resi}")
        cmd.get_wizard().apply()
    cmd.set_wizard()
    out = os.path.join(REPO, "targets", target, "structures", f"{target}_target_mutated.pdb")
    cmd.save(out, "t")
    print(f"[{target}] wrote {out}")


def main():
    for target, muts in MUTATIONS.items():
        try:
            with_pymol(target)
        except Exception as e:
            print(f"[{target}] PyMOL unavailable ({type(e).__name__}). To apply manually:")
            print(f"  pymol -qc targets/{target}/structures/{target}_target.pdb")
            for chain, resi, new in muts:
                print(f"    # mutate chain {chain} resi {resi} -> {new}")
            print("  (or use PDBFixer / Rosetta fixbb with a resfile)\n")


if __name__ == "__main__":
    sys.exit(main())
