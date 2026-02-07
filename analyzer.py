from Bio.PDB import PDBParser, PDBIO, Select
import logging
import os
import numpy as np
from scipy.spatial import KDTree

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_pdbqt_as_pdb(pdbqt_file, output_pdb):
    """
    Converts a PDBQT file to a minimal PDB format for BioPython parsing.
    Basically just keeps ATOM/HETATM lines and strips the charges/types.
    """
    try:
        with open(pdbqt_file, 'r') as f:
            lines = f.readlines()
        
        with open(output_pdb, 'w') as f:
            for line in lines:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    # PDBQT is compatible with PDB parser if we ignore the extra columns
                    # But cleaning it up is safer.
                    # We just copy the first 66 chars (up to bfactor)
                    if len(line) > 66:
                        f.write(line[:66] + "\n")
                    else:
                        f.write(line)
        return True
    except Exception as e:
        logging.error(f"Error converting PDBQT to PDB: {e}")
        return False

def analyze_docking(receptor_pdb, docked_pdbqt, residue_offset=0):
    """
    Analyzes docking results to find Hydrogen Bonds and Metal Coordination.
    
    Args:
        receptor_pdb (str): Path to receptor PDB.
        docked_pdbqt (str): Path to docked ligand PDBQT.
        residue_offset (int): Offset to apply to receptor residue numbers (default: 0).
        
    Returns:
        dict: Analysis results including bond count and details.
    """
    try:
        # Load receptor
        parser = PDBParser(QUIET=True)
        receptor = parser.get_structure("receptor", receptor_pdb)
        
        # Load ligand
        # PDBQT is close enough to PDB format for BioPython to parse standard ATOM/HETATM lines
        # We need to convert PDBQT to a temp PDB for BioPython
        temp_ligand_pdb = docked_pdbqt + ".temp.pdb"
        if not parse_pdbqt_as_pdb(docked_pdbqt, temp_ligand_pdb):
            return {"h_bond_count": 0, "metal_bond_count": 0, "details": []}
            
        ligand = parser.get_structure("ligand", temp_ligand_pdb)
        
        # 1. Build KDTree for Receptor Atoms
        receptor_atoms_coords = []
        receptor_atom_data = [] # Stores (residue object, atom object, is_metal_residue)
        
        for model in receptor:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        # We are interested in polar atoms for H-bonds and metals for coordination
                        element = atom.element.upper() if atom.element else atom.name[0].upper()
                        
                        # H-bond acceptors/donors (approximate): N, O, F, S
                        # Metals for coordination: ZN, MG, MN, FE, CA, CO
                        is_polar = element in ['N', 'O', 'F', 'S']
                        is_metal_residue = residue.get_resname().strip().upper() in ["ZN", "MG", "MN", "FE", "CA", "CO"]
                        
                        if is_polar or is_metal_residue:
                            receptor_atoms_coords.append(atom.get_coord())
                            receptor_atom_data.append((residue, atom, is_metal_residue))
                            
        if not receptor_atoms_coords:
            if os.path.exists(temp_ligand_pdb):
                os.remove(temp_ligand_pdb)
            return {"h_bond_count": 0, "metal_bond_count": 0, "details": []}
            
        receptor_tree = KDTree(receptor_atoms_coords)
        
        # 2. Check Ligand Atoms
        interactions = []
        h_bond_count = 0
        metal_bond_count = 0
        
        for model in ligand:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        element = atom.element.upper() if atom.element else atom.name[0].upper()
                        if element in ['N', 'O', 'F', 'S']: # Ligand atoms that can participate in H-bonds or metal coordination
                            # Query close receptor atoms
                            # Search slightly wider than 3.5 to catch potential interactions
                            indices = receptor_tree.query_ball_point(atom.get_coord(), 3.5)
                            
                            for idx in indices:
                                rec_res, rec_atom, rec_is_metal_residue = receptor_atom_data[idx]
                                distance = float(np.linalg.norm(atom.get_coord() - rec_atom.get_coord()))
                                
                                # Classification Logic
                                interaction_type = "Hydrogen Bond"
                                is_valid = False
                                
                                if rec_is_metal_residue:
                                    # Metal Coordination: typically shorter, strong interaction (< 2.5 A)
                                    if distance < 2.5: # Strict cutoff per user request
                                        interaction_type = "Metal Coordination"
                                        is_valid = True
                                else:
                                    # Standard H-Bond
                                    if distance < 3.5:
                                        is_valid = True
                                
                                if is_valid:
                                    # Apply Residue Offset
                                    res_num = rec_res.get_id()[1] + residue_offset
                                    
                                    # Format residue name: e.g. THR199 (O)
                                    res_label = f"{rec_res.get_resname()}{res_num}"
                                    
                                    # Check for duplicate interactions (same ligand atom -> same receptor atom)
                                    # We keep them, frontend can filter if needed.
                                    
                                    interaction = {
                                        "ligand_atom": f"{atom.get_name()} ({element})",
                                        "receptor_atom": f"{rec_atom.get_name()} ({rec_atom.element})",
                                        "residue": res_label,
                                        "distance": distance,
                                        "type": interaction_type
                                    }
                                    
                                    interactions.append(interaction)
                                    
                                    if interaction_type == "Metal Coordination":
                                        metal_bond_count += 1
                                    else:
                                        h_bond_count += 1
        
        # Cleanup
        if os.path.exists(temp_ligand_pdb):
            os.remove(temp_ligand_pdb)

        # Sort by distance
        interactions.sort(key=lambda x: x['distance'])
        
        return {
            "h_bond_count": h_bond_count,
            "metal_bond_count": metal_bond_count,
            "details": interactions
        }

    except Exception as e:
        logging.error(f"Error analyzing interactions: {e}")
        if os.path.exists(temp_ligand_pdb):
            os.remove(temp_ligand_pdb)
        return {"h_bond_count": 0, "metal_bond_count": 0, "details": []}

if __name__ == "__main__":
    # Test
    rec = "./pdbs/P00533.pdb"
    dock = "./results/test_aspirin_out.pdbqt"
    if os.path.exists(rec) and os.path.exists(dock):
        print("Running analysis...")
        res = analyze_docking(rec, dock)
        print("Analysis Result:", res)
    else:
        print("Test files not found.")
