import logging
import sys
import os

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from fetcher import fetch_alphafold_structure
from docking_engine import run_docking

logging.basicConfig(level=logging.INFO)

# Test Metal Priority with P00918 (Carbonic Anhydrase)
print("--- TESTING P00918 ---")
pdb_path1 = fetch_alphafold_structure("P00918")
res1 = run_docking(pdb_path1, "CC", job_name="test_p00918", exhaustiveness=1)
if res1:
    print("Success P00918")

# Test Cavity Finding with P29274 (Caffeine / A2A)
print("\n--- TESTING P29274 ---")
pdb_path2 = fetch_alphafold_structure("P29274")
res2 = run_docking(pdb_path2, "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", job_name="test_p29274", exhaustiveness=1)
if res2:
    print("Success P29274")
