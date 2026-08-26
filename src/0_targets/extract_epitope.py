#!/usr/bin/env python3
"""
Stage 1 - Define the design epitope from known inhibitor complexes.

For each target (NDM5, KPC3):
  1. Download the inhibitor-bound template structures from RCSB.
  2. Auto-detect the bound inhibitor (largest non-solvent HETATM group) and, for
     metallo-enzymes, the catalytic metal ions.
  3. Collect every protein residue with a heavy atom within `epitope_cutoff_A` of
     the inhibitor (union with metal-coordinating residues for NDM).
  4. Map all residues into one frame (the precursor sequence) by alignment, then
     also onto the design-scaffold PDB numbering used for generation in Stage 2.
  5. Validate that the extracted epitope contains the known catalytic machinery.
  6. Auto-curate candidate inhibitor PDBs (right enzyme + has a ligand) and fold
     them into the epitope union.

Outputs (per target, under targets/<T>/epitope/):
  - epitope_<PDB>.csv         per-complex contacting residues
  - hotspots_union.csv        merged epitope (precursor + scaffold numbering)
  - boltzgen_hotspots.json    hotspot residues in scaffold numbering for Stage 2
  - epitope_report.md         human-readable summary + validation

Runs in the base conda env (Biopython + numpy; no scipy needed).
"""
from __future__ import annotations

import json
import os
import sys
from collections import defaultdict

import numpy as np
from Bio import SeqIO

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, REPO)

import structure_io as sio  # noqa: E402

CACHE = os.path.join(REPO, "targets", "_pdb_cache")


def read_seq(path: str) -> str:
    rec = next(SeqIO.parse(os.path.join(REPO, path), "fasta"))
    return str(rec.seq).replace("\n", "").strip().upper()


def pick_binding_chain(structure, ligand: "sio.LigandHit", cutoff: float) -> str:
    """Protein chain with the most contacts to this ligand."""
    model = next(structure.get_models())
    best, best_n = None, -1
    for chain in model:
        prot = [r for r in chain if r.id[0] == " " and r.get_resname() in sio.THREE_TO_ONE]
        if not prot:
            continue
        n = len(sio.residues_near(structure, chain.id, ligand.coords(), cutoff))
        if n > best_n:
            best, best_n = chain.id, n
    return best


def epitope_from_complex(pdb_id, precursor, cutoff, keep_metals, only_resname=None):
    """Return (epitope dict keyed by precursor pos, info dict) for one PDB complex."""
    path = sio.download_pdb(pdb_id, CACHE)
    structure = sio.load_structure(path)
    title = sio.get_title(path)

    ligands = sio.detect_ligands(structure, only_resname=only_resname)
    metals = sio.detect_metals(structure, keep_metals) if keep_metals else []

    info = {"pdb": pdb_id, "title": title, "path": path,
            "ligands": [l.label for l in ligands[:5]],
            "n_metals": len(metals)}
    if not ligands and not metals:
        info["status"] = "no inhibitor ligand / metal found"
        return {}, info

    # choose the protein chain that the top ligand (or first metal) binds
    probe_ligand = ligands[0] if ligands else metals[0]
    chain = pick_binding_chain(structure, probe_ligand, cutoff)
    info["binding_chain"] = chain

    mapping = sio.map_precursor_to_chain(precursor, chain, structure)
    info["identity_to_precursor"] = round(mapping["identity"], 3)
    info["mature_start_prec"] = mapping["mature_start_prec"]

    # contacts to organic ligand(s) within cutoff
    contacts = {}
    if ligands:
        lig_coords = np.vstack([ligands[0].coords()])
        contacts.update(sio.residues_near(structure, chain, lig_coords, cutoff))
        info["inhibitor"] = ligands[0].label
    # metal-coordinating residues (tight cutoff) for di-zinc enzymes
    metal_ligands = {}
    if metals:
        mcoords = np.vstack([m.coords() for m in metals])
        metal_ligands = sio.residues_near(structure, chain, mcoords, 2.8)
        # also widen epitope to the metal pocket
        contacts.update(sio.residues_near(structure, chain, mcoords, cutoff))

    # key by precursor position
    epitope = {}
    r2p = mapping["resnum_to_prec"]
    for resnum, d in contacts.items():
        ppos = r2p.get(resnum)
        if ppos is None:
            continue
        epitope[ppos] = {"resname": d["resname"], "resnum": resnum,
                         "min_dist": d["min_dist"],
                         "metal_ligand": resnum in metal_ligands}
    info["n_epitope"] = len(epitope)
    info["metal_ligand_residues"] = sorted(
        {r2p.get(rn) for rn in metal_ligands} - {None}
    )
    info["metal_ligand_types"] = [metal_ligands[rn]["resname"] for rn in metal_ligands]
    return epitope, info


def validate_ndm(union, perfile_info):
    """Check the di-zinc machinery: >=3 distinct His + Asp + Cys among metal ligands."""
    names = [v["resname"] for v in union.values() if v.get("metal_ligand")]
    his, asp, cys = names.count("HIS"), names.count("ASP"), names.count("CYS")
    ok = his >= 3 and asp >= 1 and cys >= 1
    return ok, (f"distinct metal-coordinating residues in epitope: {his}xHIS, {asp}xASP, "
                f"{cys}xCYS (expect >=3 HIS, >=1 ASP, >=1 CYS for the B1 di-Zn site)")


def validate_kpc(union, perfile_info):
    """Check a Ser nucleophile + supporting Lys/Glu/Ser are in the epitope."""
    names = [v["resname"] for v in union.values()]
    has_ser = names.count("SER") >= 1
    has_lys = "LYS" in names
    has_glu = "GLU" in names
    ok = has_ser and (has_lys or has_glu)
    return ok, f"active-site residues in epitope: SER={names.count('SER')}, LYS={'LYS' in names}, GLU={'GLU' in names} (expect catalytic Ser + Lys/Glu)"


def run_target(tkey, tdef):
    print(f"\n{'='*70}\n{tkey}: {tdef['name']}\n{'='*70}")
    precursor = read_seq(tdef["precursor_fasta"])
    cutoff = tdef.get("epitope_cutoff_A", 5.0)
    keep_metals = tdef.get("keep_metals", [])
    out_dir = os.path.join(REPO, "targets", tkey, "epitope")
    os.makedirs(out_dir, exist_ok=True)

    # --- declared inhibitor complexes ---
    union = defaultdict(lambda: {"resname": None, "sources": set(), "min_dist": 99.9,
                                 "metal_ligand": False})
    perfile_info = []
    confirmed = []
    for spec in tdef["inhibitor_complexes"]:
        epi, info = epitope_from_complex(spec["pdb"], precursor, cutoff, keep_metals)
        perfile_info.append(info)
        print(f"  [{spec['pdb']}] {info.get('inhibitor','-')} | chain {info.get('binding_chain','?')}"
              f" | id={info.get('identity_to_precursor','?')} | epitope={info.get('n_epitope',0)}"
              f" | title='{info['title'][:50]}'")
        if not epi:
            print(f"      !! {info.get('status','no epitope')}")
            continue
        confirmed.append(spec["pdb"])
        _write_perfile(out_dir, spec["pdb"], epi)
        _merge(union, epi, spec["pdb"])

    # --- auto-curate candidate inhibitor PDBs (right enzyme + has a ligand) ---
    for cand in tdef.get("inhibitor_candidates_to_verify", []):
        try:
            epi, info = epitope_from_complex(cand, precursor, cutoff, keep_metals)
        except Exception as e:
            print(f"  [{cand}] candidate: FETCH/PARSE FAILED ({type(e).__name__})")
            continue
        ident = info.get("identity_to_precursor", 0) or 0
        if ident < 0.6:
            print(f"  [{cand}] candidate REJECTED: identity {ident} < 0.6 (different protein) "
                  f"| title='{info['title'][:50]}'")
            continue
        if not epi:
            print(f"  [{cand}] candidate REJECTED: no inhibitor ligand | title='{info['title'][:50]}'")
            continue
        print(f"  [{cand}] candidate ACCEPTED: {info.get('inhibitor','-')} id={ident} "
              f"epitope={info['n_epitope']} | title='{info['title'][:50]}'")
        perfile_info.append(info)
        confirmed.append(cand)
        _write_perfile(out_dir, cand, epi)
        _merge(union, epi, cand)

    # --- map union onto the design-scaffold numbering (for Stage 2 / BoltzGen) ---
    scaffold_pdb = tdef["design_scaffold_pdb"]
    scaffold_path = sio.download_pdb(scaffold_pdb, CACHE)
    scaffold = sio.load_structure(scaffold_path)
    scaffold_chain = tdef.get("scaffold_chain", "A")
    smap = sio.map_precursor_to_chain(precursor, scaffold_chain, scaffold)
    p2s = smap["prec_to_resnum"]

    # --- validation ---
    if keep_metals:
        ok, msg = validate_ndm(union, perfile_info)
    else:
        ok, msg = validate_kpc(union, perfile_info)

    # --- write union + boltzgen hotspots ---
    rows = []
    boltz_hotspots = []
    for ppos in sorted(union):
        v = union[ppos]
        s_resnum = p2s.get(ppos)
        rows.append({
            "precursor_pos": ppos,
            "mature_pos": (ppos - smap["mature_start_prec"] + 1) if smap["mature_start_prec"] else "",
            "scaffold_resnum": s_resnum if s_resnum is not None else "",
            "resname": v["resname"],
            "min_dist_A": round(v["min_dist"], 2),
            "metal_ligand": v["metal_ligand"],
            "sources": "|".join(sorted(v["sources"])),
        })
        if s_resnum is not None:
            boltz_hotspots.append(f"{scaffold_chain}{s_resnum}")

    union_csv = os.path.join(out_dir, "hotspots_union.csv")
    with open(union_csv, "w") as fh:
        cols = ["precursor_pos", "mature_pos", "scaffold_resnum", "resname",
                "min_dist_A", "metal_ligand", "sources"]
        fh.write(",".join(cols) + "\n")
        for r in rows:
            fh.write(",".join(str(r[c]) for c in cols) + "\n")

    boltz_json = os.path.join(out_dir, "boltzgen_hotspots.json")
    with open(boltz_json, "w") as fh:
        json.dump({
            "target": tkey,
            "design_scaffold_pdb": scaffold_pdb,
            "scaffold_chain": scaffold_chain,
            "scaffold_identity_to_precursor": round(smap["identity"], 3),
            "first_ordered_precursor_pos": smap["mature_start_prec"],
            "_note": "first_ordered_precursor_pos = first residue resolved in the crystal; "
                     "biological signal-peptide cleavage is annotated separately in Stage 0.",
            "epitope_cutoff_A": cutoff,
            "n_hotspots": len(boltz_hotspots),
            "hotspot_residues_scaffold_numbering": boltz_hotspots,
            "source_pdbs": confirmed,
        }, fh, indent=2,
            default=lambda o: int(o) if isinstance(o, np.integer)
            else float(o) if isinstance(o, np.floating) else str(o))

    _write_report(out_dir, tkey, tdef, precursor, smap, union, rows, perfile_info, confirmed, ok, msg, cutoff)

    print(f"  -> {len(union)} epitope residues; first ordered residue at precursor pos "
          f"{smap['mature_start_prec']}")
    print(f"  -> VALIDATION {'PASS' if ok else 'FAIL'}: {msg}")
    print(f"  -> wrote {union_csv}, {boltz_json}")
    return ok


def _merge(union, epi, pdb_id):
    for ppos, d in epi.items():
        u = union[ppos]
        u["resname"] = d["resname"]
        u["sources"].add(pdb_id)
        u["min_dist"] = min(u["min_dist"], d["min_dist"])
        u["metal_ligand"] = u["metal_ligand"] or d.get("metal_ligand", False)


def _write_perfile(out_dir, pdb_id, epi):
    p = os.path.join(out_dir, f"epitope_{pdb_id}.csv")
    with open(p, "w") as fh:
        fh.write("precursor_pos,resnum,resname,min_dist_A,metal_ligand\n")
        for ppos in sorted(epi):
            d = epi[ppos]
            fh.write(f"{ppos},{d['resnum']},{d['resname']},{d['min_dist']},{d['metal_ligand']}\n")


def _write_report(out_dir, tkey, tdef, precursor, smap, union, rows, perfile_info,
                  confirmed, ok, msg, cutoff):
    sp = smap["mature_start_prec"]
    lines = [
        f"# Stage 1 epitope report - {tkey} ({tdef['name']})", "",
        f"- Class: {tdef['class']}",
        f"- Precursor length: {len(precursor)} aa",
        f"- First crystallographically-ordered residue: precursor pos {sp} "
        f"(biological signal-peptide cleavage annotated in Stage 0)",
        f"- Design scaffold: {tdef['design_scaffold_pdb']} chain {tdef.get('scaffold_chain','A')} "
        f"(identity to precursor {round(smap['identity'],3)})",
        f"- Mutations vs template (already in target sequence): {', '.join(tdef.get('mutations_vs_template', []))}",
        f"- Epitope cutoff: {cutoff} A around bound inhibitor (+2.8 A around catalytic metal)",
        f"- Inhibitor complexes used: {', '.join(confirmed)}",
        "",
        f"## Validation: {'PASS' if ok else 'FAIL'}",
        f"{msg}", "",
        "## Source complexes", "",
        "| PDB | inhibitor | chain | identity | #epitope | title |",
        "|-----|-----------|-------|----------|----------|-------|",
    ]
    for inf in perfile_info:
        lines.append(f"| {inf['pdb']} | {inf.get('inhibitor','-')} | {inf.get('binding_chain','?')} "
                     f"| {inf.get('identity_to_precursor','?')} | {inf.get('n_epitope',0)} "
                     f"| {inf['title'][:60]} |")
    lines += ["", f"## Epitope / hotspot residues ({len(union)} total)", "",
              "Scaffold numbering is what Stage 2 (BoltzGen) consumes.", "",
              "| precursor | mature | scaffold | residue | min_dist (A) | metal-ligand | sources |",
              "|-----------|--------|----------|---------|--------------|--------------|---------|"]
    for r in rows:
        lines.append(f"| {r['precursor_pos']} | {r['mature_pos']} | {r['scaffold_resnum']} | "
                     f"{r['resname']} | {r['min_dist_A']} | {r['metal_ligand']} | {r['sources']} |")
    with open(os.path.join(out_dir, "epitope_report.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


def main():
    with open(os.path.join(REPO, "targets", "targets.json")) as fh:
        targets = {k: v for k, v in json.load(fh).items() if not k.startswith("_")}
    results = {}
    for tkey, tdef in targets.items():
        try:
            results[tkey] = run_target(tkey, tdef)
        except Exception as e:
            import traceback
            traceback.print_exc()
            results[tkey] = False
    print(f"\n{'='*70}\nSUMMARY: " + ", ".join(f"{k}={'PASS' if v else 'FAIL'}" for k, v in results.items()))
    return 0 if all(results.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
