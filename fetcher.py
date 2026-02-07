import os
import requests
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

from Bio.PDB import PDBParser, PDBIO, Select

class NotWater(Select):
    """
    Select class to filter out water residues (hetero flag 'W').
    """
    def accept_residue(self, residue):
        # residue.id is a tuple (hetero_flag, sequence_identifier, insertion_code)
        # Water residues typically have 'W' as the hetero flag.
        return residue.id[0] != 'W'

def clean_pdb(pdb_path):
    """
    Removes water molecules from a PDB file in-place.
    """
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("structure", pdb_path)
        
        io = PDBIO()
        io.set_structure(structure)
        io.save(pdb_path, select=NotWater())
        logging.info(f"Cleaned {pdb_path}: Removed water molecules.")
        return True
    except Exception as e:
        logging.error(f"Error cleaning PDB {pdb_path}: {e}")
        return False

def fetch_alphafold_structure(uniprot_id, output_dir="./pdbs"):
    """
    Fetches the predicted structure for a given UniProt ID from AlphaFold DB.
    
    Args:
        uniprot_id (str): The UniProt ID of the protein (e.g., 'P00533').
        output_dir (str): Directory where the PDB file should be saved.
        
    Returns:
        str: Path to the downloaded PDB file, or None if the download failed.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    # AlphaFold DB API endpoint
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    
    try:
        logging.info(f"Querying AlphaFold DB API for {uniprot_id}...")
        api_response = requests.get(api_url)
        
        if api_response.status_code != 200:
            logging.error(f"API query failed. Status code: {api_response.status_code}")
            return None
            
        data = api_response.json()
        if not data or not isinstance(data, list) or len(data) == 0:
            logging.error(f"No structure found for {uniprot_id} in AlphaFold DB.")
            return None
            
        # Get the PDB URL from the first entry (usually the most relevant one)
        pdb_url = data[0].get("pdbUrl")
        if not pdb_url:
            logging.error("PDB URL not found in API response.")
            return None
            
        output_path = os.path.join(output_dir, f"{uniprot_id}.pdb")
        
        # Check if file already exists
        if os.path.exists(output_path):
            logging.info(f"File already exists at {output_path}. Skipping download.")
            return output_path

        logging.info(f"Downloading structure from {pdb_url}...")
        
        pdb_response = requests.get(pdb_url)
        if pdb_response.status_code == 200:
            with open(output_path, "wb") as f:
                f.write(pdb_response.content)
            logging.info(f"Successfully downloaded to {output_path}")
            
            # Clean the PDB (Remove Waters)
            clean_pdb(output_path)
            
            return output_path
        else:
            logging.error(f"Failed to download PDB file. Status code: {pdb_response.status_code}")
            return None
            
    except Exception as e:
        logging.error(f"An error occurred: {e}")
        return None

if __name__ == "__main__":
    # verification
    test_id = "P00533" # EGFR
    path = fetch_alphafold_structure(test_id)
    if path and os.path.exists(path):
        print(f"VERIFICATION SUCCESS: File found at {path}")
    else:
        print("VERIFICATION FAILED")
