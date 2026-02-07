import os
import subprocess
import logging
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation
from Bio.PDB import PDBParser, PDBIO, Select

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

VINA_PATH = os.path.abspath("./bin/vina_1.2.7_mac_aarch64")

class NotHet(Select):
    def accept_residue(self, residue):
        # Exclude water (HOH) and other hetatoms (usually start with H_)
        # Keep standard residues.
        # This is a basic filter.
        return residue.id[0] == " "

def prepare_receptor(pdb_file, output_pdbqt):
    """
    Prepares a PDB file for docking by cleaning and formatting as PDBQT.
    """
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("receptor", pdb_file)
        
        # Filter non-standard residues/water
        # (For AlphaFold models, they are usually clean standard residues)
        
        with open(output_pdbqt, 'w') as f:
            for model in structure:
                for chain in model:
                    for residue in chain:
                        if residue.id[0] != " ":
                            continue # Skip hetatoms/water
                            
                        for atom in residue:
                            # Format line as PDBQT
                            # Columns:
                            # 1-6   Record name (ATOM  )
                            # 7-11  Serial number
                            # 13-16 Atom name
                            # 17    AltLoc
                            # 18-20 Residue name
                            # 22    Chain ID
                            # 23-26 ResSeq
                            # 27    iCode
                            # 31-38 X
                            # 39-46 Y
                            # 47-54 Z
                            # 55-60 Occupancy
                            # 61-66 TempFactor
                            # 71-76 Partial Charge
                            # 78-79 Atom Type
                            
                            # Get element for atom type
                            element = atom.element.upper()
                            name = atom.get_name()
                            # Basic mapping for Vina types (AD4 types)
                            # Simplified: uncharged carbon is C, aromatics A, etc.
                            # For this simulator, using Element is often a safe fallback for Vina 1.2
                            # if we don't have full typing. 
                            # However, Vina expects specific types.
                            # Let's use the element as the type.
                            
                            # Charge: 0.00 (dummy)
                            # Atom Type: Element
                            
                            # Note: Vina is strict about columns. 
                            # Charge is 71-76 (6 chars). 
                            # Atom Type is 78-79 (2 chars).
                            # There should be a space at 77? Or just contiguous? 
                            # Meeko adds a space.
                            
                            line = (
                                f"ATOM  {atom.get_serial_number():5d}  {name:<4s}{residue.get_resname():3s} {chain.get_id():1s}{residue.get_id()[1]:4d}    "
                                f"{atom.get_coord()[0]:8.3f}{atom.get_coord()[1]:8.3f}{atom.get_coord()[2]:8.3f}"
                                f"{atom.get_occupancy():6.2f}{atom.get_bfactor():6.2f}    {0.0:6.3f} {element:<2s}\n"
                            )
                            f.write(line)
        return True
    except Exception as e:
        logging.error(f"Error preparing receptor: {e}")
        return False

def get_center_and_size(pdb_file):
    """
    Calculates the geometric center and a box size covering the receptor.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("receptor", pdb_file)
    coords = []
    for model in structure:
        for atom in model.get_atoms():
             coords.append(atom.get_coord())
    
    if not coords:
        return None, None

    import numpy as np
    coords = np.array(coords)
    center = np.mean(coords, axis=0) # [x, y, z]
    
    # Calculate box size (min/max range + padding)
    min_coord = np.min(coords, axis=0)
    max_coord = np.max(coords, axis=0)
    size = (max_coord - min_coord) + 10.0 # 10A padding
    
    return center, size

def prepare_ligand(smiles, output_pdbqt):
    """
    Converts SMILES to PDBQT using RDKit and Meeko.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            logging.error("Invalid SMILES string.")
            return False
            
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        
        # Meeko preparation
        preparator = MoleculePreparation()
        preparator.prepare(mol)
        pdbqt_string = preparator.write_pdbqt_string()
        
        with open(output_pdbqt, 'w') as f:
            f.write(pdbqt_string)
            
        return True
    except Exception as e:
        logging.error(f"Error preparing ligand: {e}")
        return False

def run_docking(pdb_file, smiles, output_dir="./results", job_name="job", exhaustiveness=8):
    """
    Runs Vina docking.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    receptor_pdbqt = os.path.join(output_dir, f"{job_name}_receptor.pdbqt")
    ligand_pdbqt = os.path.join(output_dir, f"{job_name}_ligand.pdbqt")
    output_docked = os.path.join(output_dir, f"{job_name}_out.pdbqt")
    
    # 1. Prepare Receptor
    logging.info("Preparing receptor...")
    if not prepare_receptor(pdb_file, receptor_pdbqt):
        return None

    # 2. Prepare Ligand
    logging.info("Preparing ligand...")
    if not prepare_ligand(smiles, ligand_pdbqt):
        return None
        
    # 3. Calculate Center/Size
    center, size = get_center_and_size(pdb_file)
    if center is None:
        logging.error("Could not calculate center/size.")
        return None
        
    logging.info(f"Center: {center}, Size: {size}")
    
    # 4. Run Vina
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
        
        # Extract best affinity from stdout
        # Output format usually contains:
        # mode |   affinity | dist from best mode
        #      | (kcal/mol) | rmsd l.b.| rmsd u.b.
        # -----+------------+----------+----------
        #    1 |      -6.4  |      0.000 |      0.000
        
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
    # Test
    # Using the pre-downloaded P00533.pdb if available
    test_pdb = "./pdbs/P00533.pdb"
    test_smiles = "CC(=O)Oc1ccccc1C(=O)O" # Aspirin
    
    if os.path.exists(test_pdb):
        print("Running test docking...")
        # Use low exhaustiveness for quick verification
        res = run_docking(test_pdb, test_smiles, job_name="test_aspirin", exhaustiveness=1)
        print("Result:", res)
    else:
        print(f"Test PDB {test_pdb} not found. Run fetcher.py first.")
