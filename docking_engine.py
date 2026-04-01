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
    Prepares a PDB file for docking by adding hydrogens (via RDKit) and formatting as PDBQT.
    """
    import os
    from rdkit import Chem
    
    target_pdb = pdb_file
    temp_pdb_h = output_pdbqt + ".h.pdb"
    has_hydrogens = False
    
    try:
        # Try adding hydrogens with RDKit
        mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)
        if mol is None:
             # Try unsanitized if first attempt fails
             mol = Chem.MolFromPDBFile(pdb_file, removeHs=False, sanitize=False)
             
        if mol:
            mol = Chem.AddHs(mol, addCoords=True)
            Chem.MolToPDBFile(mol, temp_pdb_h)
            target_pdb = temp_pdb_h
            has_hydrogens = True
        else:
            logging.warning("RDKit failed to load receptor. Proceeding without adding hydrogens.")
            
    except Exception as e:
        logging.warning(f"Error adding hydrogens with RDKit: {e}. Proceeding without them.")

    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("receptor", target_pdb)
        
        with open(output_pdbqt, 'w') as f:
            for model in structure:
                for chain in model:
                    for residue in chain:
                        # Keep standard residues and specific metals (Zn, Mg, etc)
                        # Carbonic Anhydrase needs ZN.
                        resname = residue.get_resname().strip().upper()
                        is_metal = resname in ["ZN", "MG", "MN", "FE", "CA"]
                        
                        # Keep standard residues and specific metals (Zn, Mg, etc)
                        # Carbonic Anhydrase needs ZN.
                        resname = residue.get_resname().strip().upper()
                        is_metal = resname in ["ZN", "MG", "MN", "FE", "CA"]
                        
                        # Debug logging
                        # logging.info(f"Checking residue {resname} {residue.id}")

                        if residue.id[0] != " " and not is_metal:
                            # logging.info(f"Skipping {resname}")
                            continue 
                            
                        # logging.info(f"Keeping {resname}")
                        
                        for atom in residue:
                            element = atom.element.upper()
                            name = atom.get_name()
                            
                            # Atom Type should be title case for two-letter elements (Zn, Mg, Ca)
                            # and upper case for single letter (C, N, O).
                            # Actually AD4 types are: C, A, N, O, S, H, Zn, Mg, Mn, Fe, Ca...
                            atom_type = element.title() if len(element) > 1 else element
                            
                            # Assign +2.0 charge to metals for better electrostatics (heuristic)
                            charge = 2.0 if is_metal else 0.0
                            
                            line = (
                                f"ATOM  {atom.get_serial_number():5d}  {name:<4s}{residue.get_resname():3s} {chain.get_id():1s}{residue.get_id()[1]:4d}    "
                                f"{atom.get_coord()[0]:8.3f}{atom.get_coord()[1]:8.3f}{atom.get_coord()[2]:8.3f}"
                                f"{atom.get_occupancy():6.2f}{atom.get_bfactor():6.2f}    {charge:6.3f} {atom_type:<2s}\n"
                            )
                            f.write(line)
        
        # Cleanup
        if has_hydrogens and os.path.exists(temp_pdb_h):
            os.remove(temp_pdb_h)
            
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

def get_center_and_size(pdb_file):
    """
    Calculates center and size. Prioritizes metal ions (ZN, etc) for targeted docking.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("receptor", pdb_file)
    
    # 1. Look for Catalytic Metals (ZN, MG, MN, FE, CO, CA) to center the grid
    metal_coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname().strip().upper() in ["ZN", "MG", "MN", "FE", "CO", "CA"]:
                    for atom in residue:
                        metal_coords.append(atom.get_coord())
    
    import numpy as np
    
    if metal_coords:
        # TARGETED DOCKING: Center on the Metal (Active Site)
        logging.info("Metal (ZN) detected - Snapping box to metal center.")
        coords = np.array(metal_coords)
        center = np.mean(coords, axis=0)
        
        # For targeted docking, we don't need a huge box.
        # 22.5A is a standard high-precision box size for Vina.
        size = [22.5, 22.5, 22.5]
        return center, size

    logging.info("No metals found - Initializing Cavity Search...")

    # 2. Fallback: Cavity Search
    all_coords = []
    is_polar = []
    
    for model in structure:
        for atom in model.get_atoms():
            all_coords.append(atom.get_coord())
            element = atom.element.upper() if atom.element else atom.name[0].upper()
            name = atom.get_name().strip()
            # 1. Backbone Filter: Exclude N, CA, C, O to prevent distraction by structural core
            is_backbone = name in ['N', 'CA', 'C', 'O']
            # 2. Side-Chain Bias: Give 2x weight only to polar side-chain atoms 
            is_polar.append(1 if (element in ['N', 'O'] and not is_backbone) else 0)
            
    if not all_coords:
        return None, None
        
    all_coords_np = np.array(all_coords)
    is_polar_np = np.array(is_polar)
    
    # Grid Buffer of 5.0 A
    min_coord = np.min(all_coords_np, axis=0) - 5.0
    max_coord = np.max(all_coords_np, axis=0) + 5.0
    
    from scipy.spatial import cKDTree
    tree = cKDTree(all_coords_np)
    
    def calculate_scores(valid_pts, radius=10.0):
        neighbors = tree.query_ball_point(valid_pts, r=radius)
        scores = []
        polar_5a_counts = []
        
        # Another query for 5.0 A to check polar proximity
        neighbors_5a = tree.query_ball_point(valid_pts, r=5.0)
        
        for i, n in enumerate(neighbors):
            # 2x weight for polar side-chain, 1x weight for non-polar/backbone
            score = len(n) + np.sum(is_polar_np[n]) 
            scores.append(score)
            
            p5 = np.sum(is_polar_np[neighbors_5a[i]])
            polar_5a_counts.append(p5)
            
        final_scores = np.array(scores, dtype=float)
        # 3. Z-Axis Search: Bias ties toward extracellular top (highest Z)
        final_scores += valid_pts[:, 2] * 0.5
        
        return final_scores, np.array(polar_5a_counts)
    
    # Phase 1: Coarse 3.0A search
    x = np.arange(min_coord[0], max_coord[0], 3.0)
    y = np.arange(min_coord[1], max_coord[1], 3.0)
    z = np.arange(min_coord[2], max_coord[2], 3.0)
    xv, yv, zv = np.meshgrid(x, y, z, indexing='ij')
    grid_points = np.vstack([xv.ravel(), yv.ravel(), zv.ravel()]).T
    
    # Exclude clashing points (within 2.5A of protein)
    dist, _ = tree.query(grid_points, k=1)
    valid_points = grid_points[dist > 2.5]
    
    if len(valid_points) == 0:
        # Fallback to geometric center if something weird happens
        logging.warning("Cavity search failed to find valid points. Using geometric center.")
        center = np.mean(all_coords_np, axis=0)
        return center, [25.0, 25.0, 25.0]
        
    # Rank by polar-weighted score
    scores, polar_5a = calculate_scores(valid_points)
    sorted_indices = np.argsort(scores)[::-1]
    
    coarse_best = None
    for idx in sorted_indices:
        if polar_5a[idx] > 0:
            coarse_best = valid_points[idx]
            break
            
    if coarse_best is None:
        coarse_best = valid_points[sorted_indices[0]]
    
    # Phase 2: Refinement 1.0A search around coarse_best (e.g. within a 3.0A radius)
    rx = np.arange(coarse_best[0] - 3.0, coarse_best[0] + 3.001, 1.0)
    ry = np.arange(coarse_best[1] - 3.0, coarse_best[1] + 3.001, 1.0)
    rz = np.arange(coarse_best[2] - 3.0, coarse_best[2] + 3.001, 1.0)
    rxv, ryv, rzv = np.meshgrid(rx, ry, rz, indexing='ij')
    refine_points = np.vstack([rxv.ravel(), ryv.ravel(), rzv.ravel()]).T
    
    # Exclude clashing points
    rdist, _ = tree.query(refine_points, k=1)
    rvalid_points = refine_points[rdist > 2.5]
    
    if len(rvalid_points) > 0:
        rscores, rpolar_5a = calculate_scores(rvalid_points)
        rsorted_indices = np.argsort(rscores)[::-1]
        
        center = None
        for idx in rsorted_indices:
            if rpolar_5a[idx] > 0:
                center = rvalid_points[idx]
                break
                
        if center is None:
            center = rvalid_points[rsorted_indices[0]]
    else:
        center = coarse_best
        
    logging.info(f"Cavity search complete. Found deep pocket at [{center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f}]")
    return center, [25.0, 25.0, 25.0]

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

def run_docking(pdb_file, smiles, output_dir="./results", job_name="job", exhaustiveness=4):
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
    
    # 4. Check for Zinc and Optimize
    # specific_metals = ["ZN", "MG", "MN", "FE", "CA"]
    has_zinc = False
    with open(receptor_pdbqt, 'r') as f:
        content = f.read()
        if "ZN" in content:
            has_zinc = True
            
    if has_zinc:
        logging.info("Zinc detected in receptor. Boosting exhaustiveness to 32 for deep pocket search.")
        exhaustiveness = max(exhaustiveness, 32)

    # 5. Run Vina
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
