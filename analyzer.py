from Bio.PDB import PDBParser, NeighborSearch, PDBIO, Select
import logging
import os

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

def analyze_docking(receptor_pdb, docked_pdbqt, distance_cutoff=3.5):
    """
    Analyzes the docking result for Hydrogen Bonds.
    
    Args:
        receptor_pdb: Path to the original receptor PDB.
        docked_pdbqt: Path to the docked ligand PDBQT.
        distance_cutoff: Max distance for H-bond (default 3.5A).
        
    Returns:
        dict: Analysis results containing bond count and details.
    """
    try:
        # 1. Parse Receptor
        parser = PDBParser(QUIET=True)
        receptor_structure = parser.get_structure("receptor", receptor_pdb)
        receptor_atoms = list(receptor_structure.get_atoms())
        
        # 2. Parse Ligand
        # We need to convert PDBQT to a temp PDB for BioPython
        temp_ligand_pdb = docked_pdbqt + ".temp.pdb"
        if not parse_pdbqt_as_pdb(docked_pdbqt, temp_ligand_pdb):
            return None
            
        ligand_structure = parser.get_structure("ligand", temp_ligand_pdb)
        ligand_atoms = list(ligand_structure.get_atoms())
        
        # 3. Neighbor Search
        # Put receptor atoms in KDTree
        ns = NeighborSearch(receptor_atoms)
        
        h_bonds = []
        
        for latom in ligand_atoms:
            # We are looking for interactions between Heavy Atoms (N, O, F, S)
            # Simple geometric criteria: distance < 3.5A
            # We don't check angles here (simplified).
            
            l_elem = latom.element.upper()
            if l_elem not in ['N', 'O', 'F', 'S']:
                continue
                
            # Find neighbors
            neighbors = ns.search(latom.get_coord(), distance_cutoff)
            
            for ratom in neighbors:
                r_elem = ratom.element.upper()
                if r_elem not in ['N', 'O', 'F', 'S']:
                    continue
                    
                # Store interaction
                res = ratom.get_parent()
                res_id = f"{res.get_resname()}{res.get_id()[1]}"
                
                bond_info = {
                    "ligand_atom": f"{latom.get_name().strip()} ({l_elem})",
                    "receptor_atom": f"{ratom.get_name().strip()} ({r_elem})",
                    "residue": res_id,
                    "distance": float(latom - ratom)
                }
                
                # Filter duplicates (e.g. if multiple models or whatever)
                h_bonds.append(bond_info)
        
        # Cleanup
        if os.path.exists(temp_ligand_pdb):
            os.remove(temp_ligand_pdb)
            
        return {
            "h_bond_count": len(h_bonds),
            "details": h_bonds
        }
        
    except Exception as e:
        logging.error(f"Error during analysis: {e}")
        return None

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
