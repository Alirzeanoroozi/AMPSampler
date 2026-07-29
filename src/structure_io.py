"""
Shared structure I/O for the inhibitor-binder pipeline.

Dependency-light on purpose: Biopython + numpy only (no scipy / no network libs
beyond urllib), so it runs in the base conda env present on this machine.

Provides:
  - download_pdb(pdb_id)                : fetch from RCSB (cached on disk)
  - load_structure(path)                : Biopython Structure
  - chain_sequence(structure, chain)    : (seq_str, [residue, ...]) for amino acids
  - map_precursor_to_chain(...)         : align precursor seq <-> template chain,
                                          return resnum<->precursor-position maps
  - detect_ligands(structure)           : organic ligands (non-solvent HETATM groups)
  - detect_metals(structure, elements)  : metal ions (e.g. ZN)
  - residues_near(structure, chain, probe_atoms, cutoff): contacting residues
"""
from __future__ import annotations

import os
import urllib.request
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser, Selection
from Bio.Align import PairwiseAligner

THREE_TO_ONE = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
    # common modified residues seen near active sites
    "MSE": "M", "SEP": "S", "TPO": "T", "PTR": "Y", "CSO": "C", "KCX": "K",
}

# HETATM residue names that are solvent / buffer / crystallization additives / ions.
# Anything in here is NOT treated as an inhibitor ligand.
SOLVENT_AND_IONS = {
    "HOH", "DOD", "WAT",
    "NA", "K", "CL", "CA", "MG", "MN", "FE", "FE2", "CD", "CU", "NI", "CO",
    "ZN",  # metals handled separately via detect_metals
    "SO4", "PO4", "PI", "PO3", "VO4", "WO4", "MOO",
    "GOL", "EDO", "PEG", "PG4", "PGE", "1PE", "2PE", "P6G", "MPD", "BU3",
    "ACT", "FMT", "ACY", "EPE", "MES", "TRS", "BTB", "BIS", "TAR", "CIT", "FLC",
    "DMS", "BME", "MRD", "IMD", "IPA", "DTT", "TLA", "MLA", "OXL",
    "IOD", "BR", "F", "NO3", "NH4", "CO3", "CAC", "AZI", "SCN", "PER",
    "GLC", "NAG", "MAN", "BMA", "FUC",  # generic sugars (rarely the inhibitor here)
}

ELEMENT_FROM_RESNAME = {"ZN": "ZN", "MG": "MG", "MN": "MN", "FE": "FE", "CA": "CA"}

CACHE_DIR_DEFAULT = "targets/_pdb_cache"


@dataclass
class LigandHit:
    resname: str
    chain: str
    resseq: int
    icode: str
    n_atoms: int
    atoms: list = field(default_factory=list)  # list of Bio.PDB Atom

    @property
    def label(self) -> str:
        ic = self.icode.strip()
        return f"{self.resname}/{self.chain}{self.resseq}{ic}"

    def coords(self) -> np.ndarray:
        return np.array([a.coord for a in self.atoms], dtype=float)


def download_pdb(pdb_id: str, cache_dir: str = CACHE_DIR_DEFAULT) -> str:
    """Download <pdb_id>.pdb from RCSB into cache_dir (skip if present). Falls back to .cif."""
    pdb_id = pdb_id.strip().upper()
    os.makedirs(cache_dir, exist_ok=True)
    pdb_path = os.path.join(cache_dir, f"{pdb_id}.pdb")
    if os.path.exists(pdb_path) and os.path.getsize(pdb_path) > 0:
        return pdb_path
    for ext, url in (
        ("pdb", f"https://files.rcsb.org/download/{pdb_id}.pdb"),
        ("cif", f"https://files.rcsb.org/download/{pdb_id}.cif"),
    ):
        out = os.path.join(cache_dir, f"{pdb_id}.{ext}")
        try:
            urllib.request.urlretrieve(url, out)
            if os.path.getsize(out) > 0:
                return out
        except Exception:
            continue
    raise RuntimeError(f"Could not download {pdb_id} from RCSB.")


def load_structure(path: str):
    parser = MMCIFParser(QUIET=True) if path.lower().endswith((".cif", ".mmcif")) else PDBParser(QUIET=True)
    return parser.get_structure(os.path.basename(path), path)


def get_title(path: str) -> str:
    """Best-effort structure title (for verifying which inhibitor a PDB actually contains)."""
    title = []
    try:
        with open(path) as fh:
            for line in fh:
                if line.startswith("TITLE"):
                    title.append(line[10:].rstrip())
                elif line.startswith(("ATOM", "HETATM", "_")):
                    break
                elif line.startswith("_struct.title"):
                    title.append(line.split(None, 1)[-1].strip().strip("'\""))
    except Exception:
        pass
    return " ".join(t.strip() for t in title).strip()


def _first_model(structure):
    return next(structure.get_models())


def chain_sequence(structure, chain_id: str) -> Tuple[str, List]:
    """Return (one-letter sequence, [residue,...]) of standard amino acids in chain order."""
    model = _first_model(structure)
    if chain_id not in model:
        raise KeyError(f"Chain {chain_id} not in structure (have {[c.id for c in model]}).")
    seq, residues = [], []
    for res in model[chain_id]:
        if res.id[0] != " ":  # skip HETATM/water
            continue
        one = THREE_TO_ONE.get(res.get_resname())
        if one is None:
            continue
        seq.append(one)
        residues.append(res)
    return "".join(seq), residues


def _aligner() -> PairwiseAligner:
    a = PairwiseAligner()
    a.mode = "global"
    a.match_score = 2.0
    a.mismatch_score = -1.0
    a.open_gap_score = -5.0
    a.extend_gap_score = -0.5
    # End-gap-free (semiglobal). Set the four side sub-scores directly at both
    # open and extend. Don't probe aggregate names with hasattr(): once open/extend
    # differ, the getter for properties like end_insertion_score raises ValueError
    # ("gap scores are different"), which hasattr() does not suppress.
    for side in ("target_left", "target_right", "query_left", "query_right"):
        for kind in ("open_gap_score", "extend_gap_score"):
            setattr(a, f"{side}_{kind}", 0.0)
    return a


def map_precursor_to_chain(precursor_seq: str, chain_id: str, structure) -> Dict:
    """
    Align the full precursor sequence to a template chain's modeled sequence.

    Returns dict with:
      resnum_to_prec : {pdb_resnum (int) -> precursor position (1-based)}
      prec_to_resnum : reverse
      mature_start_prec : 1-based precursor position of the first modeled residue
                          (== inferred signal-peptide length + 1)
      identity : fraction identical over aligned columns
    """
    chain_seq, residues = chain_sequence(structure, chain_id)
    aln = _aligner().align(precursor_seq, chain_seq)[0]
    # aln.aligned -> blocks of (start,end) index pairs into (precursor, chain)
    resnum_to_prec, prec_to_resnum = {}, {}
    n_match = n_cols = 0
    for (p0, p1), (c0, c1) in zip(aln.aligned[0], aln.aligned[1]):
        for k in range(int(p1) - int(p0)):
            p_idx = int(p0) + k     # 0-based into precursor
            c_idx = int(c0) + k     # 0-based into chain residues list
            res = residues[c_idx]
            resnum = int(res.id[1])
            prec_pos = p_idx + 1    # 1-based
            resnum_to_prec[resnum] = prec_pos
            prec_to_resnum[prec_pos] = resnum
            n_cols += 1
            if precursor_seq[p_idx] == chain_seq[c_idx]:
                n_match += 1
    mature_start_prec = min(resnum_to_prec.values()) if resnum_to_prec else None
    # mature_start_prec above is wrong direction; recompute as smallest precursor pos mapped
    mature_start_prec = min(prec_to_resnum.keys()) if prec_to_resnum else None
    return {
        "resnum_to_prec": resnum_to_prec,
        "prec_to_resnum": prec_to_resnum,
        "mature_start_prec": mature_start_prec,
        "identity": (n_match / n_cols) if n_cols else 0.0,
        "n_aligned": n_cols,
        "chain_seq": chain_seq,
    }


def detect_ligands(structure, min_atoms: int = 8, only_resname: Optional[str] = None) -> List[LigandHit]:
    """Organic (non-solvent, non-metal, non-amino-acid) HETATM groups, largest first."""
    model = _first_model(structure)
    hits: List[LigandHit] = []
    for chain in model:
        for res in chain:
            hetflag, resseq, icode = res.id
            if hetflag == " ":
                continue
            name = res.get_resname().strip()
            if name in SOLVENT_AND_IONS or name in THREE_TO_ONE:
                continue
            if only_resname and name != only_resname:
                continue
            atoms = [a for a in res.get_atoms() if a.element != "H"]
            if len(atoms) < min_atoms:
                continue
            hits.append(LigandHit(name, chain.id, resseq, icode, len(atoms), atoms))
    hits.sort(key=lambda h: h.n_atoms, reverse=True)
    return hits


def detect_metals(structure, elements: Sequence[str] = ("ZN",)) -> List[LigandHit]:
    elements = {e.upper() for e in elements}
    model = _first_model(structure)
    hits: List[LigandHit] = []
    for chain in model:
        for res in chain:
            if res.id[0] == " ":
                continue
            name = res.get_resname().strip().upper()
            if name in elements or ELEMENT_FROM_RESNAME.get(name) in elements:
                atoms = list(res.get_atoms())
                hits.append(LigandHit(name, chain.id, res.id[1], res.id[2], len(atoms), atoms))
    return hits


def residues_near(structure, chain_id: str, probe_coords: np.ndarray, cutoff: float) -> Dict[int, dict]:
    """
    Protein residues in chain_id with any heavy atom within `cutoff` of any probe coord.
    Returns {pdb_resnum: {resname, min_dist, n_contacts}}.
    """
    _, residues = chain_sequence(structure, chain_id)
    probe = np.asarray(probe_coords, dtype=float)
    if probe.ndim == 1:
        probe = probe[None, :]
    out: Dict[int, dict] = {}
    for res in residues:
        atoms = np.array([a.coord for a in res.get_atoms() if a.element != "H"], dtype=float)
        if atoms.size == 0:
            continue
        # pairwise distances (small structures -> brute force is fine)
        d = np.sqrt(((atoms[:, None, :] - probe[None, :, :]) ** 2).sum(-1))
        dmin = float(d.min())
        if dmin <= cutoff:
            out[res.id[1]] = {
                "resname": res.get_resname(),
                "min_dist": round(dmin, 2),
                "n_contacts": int((d <= cutoff).sum()),
            }
    return out
