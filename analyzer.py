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
        receptor_atom_data = [] 
        
        for model in receptor:
            for chain in model:
                for residue in chain:
                    resname = residue.get_resname().strip().upper()
                    for atom in residue:
                        atom_name = atom.get_name().strip().upper()
                        # Extract element carefully
                        element = atom.element.upper() if atom.element else "".join(filter(str.isalpha, atom_name))
                        
                        is_polar = element in ['N', 'O', 'F', 'S']
                        is_metal = resname in ["ZN", "MG", "MN", "FE", "CA", "CO"]
                        is_acidic = resname in ["ASP", "GLU"] and atom_name in ["OD1", "OD2", "OE1", "OE2"]
                        is_hydrophobic = resname in ["VAL", "LEU", "ILE", "MET", "PHE", "TRP", "ALA"] and element in ["C", "A"]
                        
                        if is_polar or is_metal or is_acidic or is_hydrophobic:
                            receptor_atoms_coords.append(atom.get_coord())
                            receptor_atom_data.append({
                                'residue': residue,
                                'atom': atom,
                                'is_metal': is_metal,
                                'is_polar': is_polar,
                                'is_acidic': is_acidic,
                                'is_hydrophobic': is_hydrophobic,
                                'element': element,
                                'resname': resname
                            })
                            
        if not receptor_atoms_coords:
            if os.path.exists(temp_ligand_pdb):
                os.remove(temp_ligand_pdb)
            return {"h_bond_count": 0, "metal_bond_count": 0, "salt_bridge_count": 0, "halogen_bond_count": 0, "hydrophobic_saturation": 0, "details": []}
            
        receptor_tree = KDTree(receptor_atoms_coords)
        
        # 2. Check Ligand Atoms
        interactions = []
        h_bond_count = 0
        metal_bond_count = 0
        salt_bridge_count = 0
        halogen_bond_count = 0
        
        unique_ligand_carbons = set()
        total_ligand_carbons = 0

        for model in ligand:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        atom_name = atom.get_name().strip().upper()
                        element = atom.element.upper() if atom.element else "".join(filter(str.isalpha, atom_name))
                        
                        if element in ['C', 'A']:
                            total_ligand_carbons += 1
                            
                        # Search receptor atoms up to 4.0A
                        indices = receptor_tree.query_ball_point(atom.get_coord(), 4.0)
                        
                        for idx in indices:
                            data = receptor_atom_data[idx]
                            rec_res = data['residue']
                            rec_atom = data['atom']
                            distance = float(np.linalg.norm(atom.get_coord() - rec_atom.get_coord()))
                            
                            # HYDROPHOBIC SATURATION
                            if element in ['C', 'A'] and data['is_hydrophobic']:
                                if distance < 4.0:
                                    unique_ligand_carbons.add(id(atom))

                            # Determine interaction type
                            interaction_type = None
                            
                            if element == 'N' and data['is_acidic']:
                                if distance < 4.0:
                                    interaction_type = "Salt Bridge"
                            elif data['is_metal']:
                                if distance <= 3.2:
                                    interaction_type = "Metal Coordination"
                                elif distance <= 3.5:
                                    interaction_type = "Weak Coordination"
                            elif element in ['F', 'CL', 'BR', 'I'] and data['element'] in ['N', 'O']:
                                if distance < 3.5:
                                    interaction_type = "Halogen Bond"
                            elif element in ['N', 'O', 'F', 'S'] and data['is_polar']:
                                if distance < 3.5:
                                    interaction_type = "Hydrogen Bond"

                            if interaction_type:
                                res_num = rec_res.get_id()[1] + residue_offset
                                res_label = f"{data['resname']}{res_num}:{rec_atom.get_name().strip()}"
                                
                                interaction = {
                                    "ligand_atom": atom.get_name().strip(),
                                    "receptor_atom": rec_atom.get_name().strip(),
                                    "residue": res_label,
                                    "distance": distance,
                                    "type": interaction_type
                                }
                                interactions.append(interaction)
                                
                                if interaction_type == "Metal Coordination":
                                    metal_bond_count += 1
                                elif interaction_type == "Salt Bridge":
                                    salt_bridge_count += 1
                                elif interaction_type == "Halogen Bond":
                                    halogen_bond_count += 1
                                elif interaction_type == "Hydrogen Bond":
                                    h_bond_count += 1
        
        hydrophobic_saturation = 0
        if total_ligand_carbons > 0:
            hydrophobic_saturation = (len(unique_ligand_carbons) / total_ligand_carbons) * 100

        # Cleanup
        if os.path.exists(temp_ligand_pdb):
            os.remove(temp_ligand_pdb)

        # Sort by distance
        interactions.sort(key=lambda x: x['distance'])
        
        return {
            "h_bond_count": h_bond_count,
            "metal_bond_count": metal_bond_count,
            "salt_bridge_count": salt_bridge_count,
            "halogen_bond_count": halogen_bond_count,
            "hydrophobic_saturation": round(hydrophobic_saturation, 1),
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
