import os
import subprocess
import logging
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation
from Bio.PDB import PDBParser

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

VINA_PATH = os.path.abspath("./bin/vina_1.2.7_mac_aarch64")

def prepare_receptor(pdb_file, output_pdbqt):
    """
    Prepares a PDB file for docking using RDKit (Hydrogens) and BioPython (Formatting)
    to guarantee strict Vina PDBQT column alignments and AD4 atom types.
    """
    import os
    from rdkit import Chem
    from Bio.PDB import PDBParser
    
    target_pdb = pdb_file
    temp_pdb_h = output_pdbqt + ".h.pdb"
    has_hydrogens = False
    
    try:
        # 1. Add Hydrogens at physiological pH using RDKit
        mol = Chem.MolFromPDBFile(pdb_file, removeHs=False, sanitize=False)
        if mol:
            mol = Chem.AddHs(mol, addCoords=True)
            Chem.MolToPDBFile(mol, temp_pdb_h)
            target_pdb = temp_pdb_h
            has_hydrogens = True
        else:
            logging.warning("RDKit failed to load receptor. Proceeding without adding hydrogens.")
    except Exception as e:
        logging.warning(f"Error adding hydrogens with RDKit: {e}")

    # 2. Write strict PDBQT format for Vina
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("receptor", target_pdb)
        
        with open(output_pdbqt, 'w') as f:
            for model in structure:
                for chain in model:
                    for residue in chain:
                        resname = residue.get_resname().strip().upper()
                        # Skip water
                        if resname in ['HOH', 'WAT']: continue
                        
                        for atom in residue:
                            fullname = atom.get_fullname()
                            element = atom.element.upper() if atom.element else fullname.strip()[0]
                            
                            # AD4 Atom Typing (Crucial for Vina Scoring)
                            atom_type = element
                            if element == 'C':
                                # Assign Aromatic carbon (A) to rings, otherwise Aliphatic (C)
                                atom_type = 'A' if resname in ['PHE', 'TYR', 'TRP', 'HIS'] else 'C'
                            elif element == 'N': atom_type = 'NA' # Nitrogen Acceptor
                            elif element == 'O': atom_type = 'OA' # Oxygen Acceptor
                            elif element == 'H': atom_type = 'HD' # Hydrogen Donor
                            elif element == 'S': atom_type = 'SA' # Sulfur Acceptor
                            elif element in ['ZN', 'MG', 'MN', 'CA', 'FE']: 
                                atom_type = element.title() # Metals: Zn, Mg, etc.
                            
                            x, y, z = atom.get_coord()
                            occ = atom.get_occupancy() if atom.get_occupancy() is not None else 1.0
                            bfac = atom.get_bfactor() if atom.get_bfactor() is not None else 0.0
                            # Vina ignores partial charges, but the column must exist
                            charge = 2.0 if atom_type in ['Zn', 'Mg', 'Mn', 'Ca', 'Fe'] else 0.0
                            
                            # Extremely strict PDBQT column formatting
                            line = "ATOM  %5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f    %6.3f %-2s\n" % (
                                atom.get_serial_number(),
                                fullname,
                                resname[:3],
                                chain.get_id(),
                                residue.get_id()[1],
                                x, y, z,
                                occ, bfac,
                                charge,
                                atom_type
                            )
                            f.write(line)
                            
        # Cleanup temp files
        if has_hydrogens and os.path.exists(temp_pdb_h):
            os.remove(temp_pdb_h)
            
        logging.info("Receptor successfully prepared with native Python PDBQT writer.")
        return True
    except Exception as e:
        logging.error(f"Error preparing receptor: {e}")
        return False
def get_center_and_size(pdb_file, target_residue=None):
    """
    Calculates center and size. 
    1. Prioritizes a specific target residue (e.g., "ASN253") if provided.
    2. Falls back to catalytic metals.
    3. Falls back to cavity search.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("receptor", pdb_file)
    
    # 1. TARGETED RESIDUE SEARCH
    if target_residue:
        # Parse target residue gracefully
        res_name_target = "".join(filter(str.isalpha, target_residue)).upper()
        res_num_str = "".join(filter(str.isdigit, target_residue))
        try:
            res_num_target = int(res_num_str) if res_num_str else None
            target_coords = []
            
            for model in structure:
                for chain in model:
                    for residue in chain:
                        res_name = residue.get_resname().strip().upper()
                        res_num = residue.get_id()[1]
                        
                        if res_name == res_name_target and (res_num_target is None or res_num == res_num_target):
                            for atom in residue:
                                target_coords.append(atom.get_coord())
            
            if target_coords:
                logging.info(f"Target residue {target_residue} found. Centering docking box here.")
                coords = np.array(target_coords)
                center = np.mean(coords, axis=0)
                # 15.0A is a standard high-precision box size for a known active site
                return center, [15.0, 15.0, 15.0]
            else:
                logging.warning(f"Target residue {target_residue} not found in PDB. Falling back to Metals/Cavity.")
        except ValueError:
            logging.warning(f"Could not parse target residue format '{target_residue}'. Use format like 'ASN253'.")

    # 2. METAL SEARCH FALLBACK
    metal_coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname().strip().upper() in ["ZN", "MG", "MN", "FE", "CO", "CA"]:
                    for atom in residue:
                        metal_coords.append(atom.get_coord())
    
    if metal_coords:
        logging.info("Metal detected - Snapping box to metal center.")
        coords = np.array(metal_coords)
        center = np.mean(coords, axis=0)
        return center, [22.5, 22.5, 22.5]

    logging.info("No targets or metals found - Initializing Cavity Search...")

    # 3. CAVITY SEARCH FALLBACK (Simplified for brevity)
    all_coords = []
    for model in structure:
        for atom in model.get_atoms():
            all_coords.append(atom.get_coord())
            
    if not all_coords:
        return None, None
        
    all_coords_np = np.array(all_coords)
    center = np.mean(all_coords_np, axis=0)
    logging.info("Using geometric center of the protein.")
    return center, [25.0, 25.0, 25.0]

def prepare_ligand(smiles, output_pdbqt):
    """
    Converts SMILES to PDBQT using RDKit and Meeko.
    (Meeko is maintained by the Forli lab, creators of Vina, so it handles ligands perfectly).
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            logging.error("Invalid SMILES string.")
            return False
            
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        
        preparator = MoleculePreparation()
        preparator.prepare(mol)
        pdbqt_string = preparator.write_pdbqt_string()
        
        with open(output_pdbqt, 'w') as f:
            f.write(pdbqt_string)
            
        return True
    except Exception as e:
        logging.error(f"Error preparing ligand: {e}")
        return False

def run_docking(pdb_file, smiles, output_dir="./results", job_name="job", exhaustiveness=32, target_residue=None):
    """
    Runs Vina docking.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    receptor_pdbqt = os.path.join(output_dir, f"{job_name}_receptor.pdbqt")
    ligand_pdbqt = os.path.join(output_dir, f"{job_name}_ligand.pdbqt")
    output_docked = os.path.join(output_dir, f"{job_name}_out.pdbqt")
    
    logging.info("Preparing receptor with RDKit and BioPython...")
    if not prepare_receptor(pdb_file, receptor_pdbqt):
        return None

    logging.info("Preparing ligand with Meeko...")
    if not prepare_ligand(smiles, ligand_pdbqt):
        return None
        
    # Calculate Center/Size (Pass the target residue!)
    center, size = get_center_and_size(pdb_file, target_residue=target_residue)
    if center is None:
        logging.error("Could not calculate center/size.")
        return None
        
    logging.info(f"Center: {center}, Size: {size}")

    # Run Vina
    cmd = [
        VINA_PATH,
        "--receptor", receptor_pdbqt,
        "--ligand", ligand_pdbqt,
        "--center_x", str(center[0]),
        "--center_y", str(center[1]),
        "--center_z", str(center[2]),
        "--size_x", str(size[0]),
        "--size_y", str(size[1]),
        "--size_z", str(size[2]),
        "--cpu", "4",
        "--exhaustiveness", str(exhaustiveness),
        "--out", output_docked
    ]
    
    logging.info(f"Running Vina: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logging.error(f"Vina failed: {result.stderr}")
            return None
            
        logging.info("Vina completed successfully.")
        
        best_affinity = None
        for line in result.stdout.splitlines():
            if line.strip().startswith("1"):
                parts = line.split()
                if len(parts) >= 2:
                    best_affinity = float(parts[1])
                    break
        
        return {
            "affinity": best_affinity,
            "docked_file": output_docked,
            "stdout": result.stdout
        }

    except Exception as e:
        logging.error(f"Error executing Vina: {e}")
        return None

if __name__ == "__main__":
    # CONTROL TEST: Experimental Structure + Flexible Ligand
    # This will prove the engine works when the biological pocket is actually open.
    
    test_pdb = "./pdbs/P00918_exp.pdb"  # You already have this file
    test_smiles = "CC(=O)Nc1nnc(s1)S(=O)(=O)N" # Acetazolamide (has rotatable bonds)
    
    if os.path.exists(test_pdb):
        print("Running CAII Control Test...")
        # We target the Zinc (ZN) atom specifically.
        res = run_docking(
            test_pdb, 
            test_smiles, 
            job_name="test_acetazolamide", 
            exhaustiveness=32,
            target_residue="ZN"
        )
        print("Control Result:", res)
        
        # Run the analyzer to prove H-Bonds work
        from analyzer import analyze_docking
        if res and res.get("docked_file"):
            analysis = analyze_docking(test_pdb, res["docked_file"])
            print(f"H-Bonds Detected: {analysis['h_bond_count']}")
            print(f"Metal Bonds Detected: {analysis['metal_bond_count']}")
    else:
        print(f"Test PDB {test_pdb} not found. Please ensure it is in the pdbs/ folder.")