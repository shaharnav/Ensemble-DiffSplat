import logging
from fetcher import fetch_alphafold_structure
from docking_engine import run_docking
from analyzer import analyze_docking

logging.basicConfig(level=logging.INFO)

tests = [
    {
        "name": "Streptavidin-Biotin",
        "uniprot": "P22629",
        "smiles": "O=C(O)CCCC1SC2C(NC(=O)N2)C1", 
        "target_residue": "TRP120"
    },
    {
        "name": "Trypsin-Benzamidine",
        "uniprot": "P00760",
        "smiles": "NC(=N)c1ccccc1", 
        "target_residue": "SER195" 
    },
    {
        "name": "EGFR-Gefitinib",
        "uniprot": "P00533",
        "smiles": "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN4CCOCC4", 
        "target_residue": "MET793"
    }
]

for t in tests:
    print(f"\n--- TESTING {t['name']} ---")
    pdb_path = fetch_alphafold_structure(t["uniprot"])
    if not pdb_path:
        print(f"Failed to fetch {t['uniprot']}")
        continue
        
    job_name = f"test_{t['uniprot']}"
    res = run_docking(pdb_path, t["smiles"], output_dir="./results", job_name=job_name, exhaustiveness=8, target_residue=t["target_residue"])
    
    if res:
        analysis = analyze_docking(pdb_path, res["docked_file"])
        print(f"[{t['name']}] Affinity: {res['affinity']}")
        print(f"  Hydrophobic Saturation: {analysis.get('hydrophobic_saturation', 0)}%")
        print(f"  Salt Bridges: {analysis.get('salt_bridge_count', 0)}")
        print(f"  Halogen Bonds: {analysis.get('halogen_bond_count', 0)}")
        print(f"  H-Bonds: {analysis.get('h_bond_count', 0)}")
        
        # Display details summary
        for bond in analysis.get('details', []):
            print(f"    - {bond['type']}: {bond['ligand_atom']} <-> {bond['residue']} at {bond['distance']:.2f} Å")
