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
    
    # 1. Grab the target residue from the frontend, or default to ASN253 for our benchmark
    target_residue = data.get("target_residue", "ASN253") 
    
    if not uniprot_id or not smiles:
        return jsonify({"error": "Missing UniProt ID or SMILES"}), 400

    job_id = f"{uniprot_id}_{abs(hash(smiles))}"
    
    # Fetch Structure
    pdb_path = fetch_alphafold_structure(uniprot_id)
    if not pdb_path:
        return jsonify({"error": "Failed to fetch structure from AlphaFold DB"}), 404
        
    # 2. Pass target_residue into run_docking!
    docking_result = run_docking(
        pdb_path, 
        smiles, 
        output_dir=RESULT_DIR, 
        job_name=job_id, 
        exhaustiveness=32,
        target_residue=target_residue  # <--- THIS WAS MISSING
    )
    
    if not docking_result:
        return jsonify({"error": "Docking failed"}), 500
        
    # Analyze
    residue_offset = int(data.get("residue_offset", 0))
    analysis = analyze_docking(pdb_path, docking_result["docked_file"], residue_offset=residue_offset)
    
    docked_file = docking_result["docked_file"]
    if docked_file.startswith("./"):
        docked_file = docked_file[2:]

    response = {
        "affinity": docking_result["affinity"],
        "docked_file": docked_file,
        "interaction_summary": {
            "hydrophobic_saturation": analysis.get("hydrophobic_saturation", 0),
            "salt_bridge_count": analysis.get("salt_bridge_count", 0),
            "halogen_bond_count": analysis.get("halogen_bond_count", 0)
        },
        "h_bond_count": analysis.get("h_bond_count", 0),
        "metal_bond_count": analysis.get("metal_bond_count", 0),
        "interaction_details": analysis.get("details", [])
    }
    
    return jsonify(response)
@app.route("/results/<path:filename>")
def serve_results(filename):
    return send_from_directory(RESULT_DIR, filename)

@app.route("/pdbs/<path:filename>")
def serve_pdbs(filename):
    return send_from_directory("./pdbs", filename)

@app.route("/", defaults={'path': ''})
@app.route("/<path:path>")
def serve(path):
    if path != "" and os.path.exists(app.static_folder + '/' + path):
        return send_from_directory(app.static_folder, path)
    else:
        return send_from_directory(app.static_folder, 'index.html')

if __name__ == "__main__":
    app.run(port=5001, debug=True)
