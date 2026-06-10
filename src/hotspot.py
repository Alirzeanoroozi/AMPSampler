import numpy as np
from Bio.PDB import PDBParser, MMCIFParser, Selection
from scipy.spatial import cKDTree
import pandas as pd
import os

three_to_one_map = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

def hotspot_residues(trajectory_file, binder_chain="B", atom_distance_cutoff=4.0):
    # Detect file type for parsing
    if trajectory_file.lower().endswith('.cif') or trajectory_file.lower().endswith('.mmcif'):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", trajectory_file)

    # Get the specified binder chain
    binder_atoms = Selection.unfold_entities(structure[0][binder_chain], 'A')
    binder_coords = np.array([atom.coord for atom in binder_atoms])

    # Get atoms and coords for the target chain
    target_atoms = Selection.unfold_entities(structure[0]['A'], 'A')
    target_coords = np.array([atom.coord for atom in target_atoms])

    # Build KD trees for both chains
    binder_tree = cKDTree(binder_coords)
    target_tree = cKDTree(target_coords)

    # Prepare to collect interacting residues
    interacting_residues = {}

    # Query the tree for pairs of atoms within the distance cutoff
    pairs = binder_tree.query_ball_tree(target_tree, atom_distance_cutoff)

    # Process each binder atom's interactions
    for binder_idx, close_indices in enumerate(pairs):
        binder_residue = binder_atoms[binder_idx].get_parent()
        binder_resname = binder_residue.get_resname()

        # Convert three-letter code to single-letter code using the manual dictionary
        if binder_resname in three_to_one_map:
            aa_single_letter = three_to_one_map[binder_resname]
            for close_idx in close_indices:
                target_residue = target_atoms[close_idx].get_parent()
                interacting_residues[binder_residue.id[1]] = aa_single_letter

    return interacting_residues

def get_file_path(x):
    for file in config:
        if x in file:
            print(file)
            return f"../workbench/NDM5_run/final_ranked_designs/final_10000_designs/{file}"
    return None

config = os.listdir("../workbench/NDM5_run/final_ranked_designs/final_10000_designs")
df = pd.read_csv("../workbench/NDM5_run/final_ranked_designs/final_designs_metrics_10000.csv")
new_df = df[['id', 'final_rank', 'designed_sequence', 'quality_score']]

new_df['hotspot_residues_A'] = new_df['id'].apply(lambda x: list(hotspot_residues(get_file_path(x), "A").keys()))
new_df['hotspot_residues_B'] = new_df['id'].apply(lambda x: list(hotspot_residues(get_file_path(x), "B").keys()))
new_df.to_csv("AMPS.csv", index=False)
