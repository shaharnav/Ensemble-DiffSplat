from flask import Flask, request, jsonify, send_from_directory
from flask_cors import CORS
import os
import logging
from fetcher import fetch_alphafold_structure
from docking_engine import run_docking
from analyzer import analyze_docking

app = Flask(__name__, static_folder="frontend/dist", static_url_path="")
CORS(app) # Allow cross-origin requests for dev

# Configure logging
logging.basicConfig(level=logging.INFO)

RESULT_DIR = "./results"
if not os.path.exists(RESULT_DIR):
    os.makedirs(RESULT_DIR)

@app.route("/api/dock", methods=["POST"])
def dock():
    data = request.json
    uniprot_id = data.get("uniprot_id")
    smiles = data.get("smiles")
    
    if not uniprot_id or not smiles:
        return jsonify({"error": "Missing UniProt ID or SMILES"}), 400

    job_id = f"{uniprot_id}_{abs(hash(smiles))}"
    
    # 1. Fetch Structure
    pdb_path = fetch_alphafold_structure(uniprot_id)
    if not pdb_path:
        return jsonify({"error": "Failed to fetch structure from AlphaFold DB"}), 404
        
    # 2. Run Docking
    # Use standard exhaustiveness for real runs, maybe 8? 
    # For speed in this demo, let's keep it moderate (4).
    docking_result = run_docking(pdb_path, smiles, output_dir=RESULT_DIR, job_name=job_id, exhaustiveness=4)
    
    if not docking_result:
        return jsonify({"error": "Docking failed"}), 500
        
    # 3. Analyze
    analysis = analyze_docking(pdb_path, docking_result["docked_file"])
    
    response = {
        "affinity": docking_result["affinity"],
        "h_bonds": analysis["h_bond_count"] if analysis else 0,
        "h_bond_details": analysis["details"] if analysis else [],
        "docked_file": docking_result["docked_file"],
        "stdout": docking_result["stdout"]
    }
    
    return jsonify(response)

@app.route("/", defaults={'path': ''})
@app.route("/<path:path>")
def serve(path):
    if path != "" and os.path.exists(app.static_folder + '/' + path):
        return send_from_directory(app.static_folder, path)
    else:
        return send_from_directory(app.static_folder, 'index.html')

if __name__ == "__main__":
    app.run(port=5001, debug=True)
