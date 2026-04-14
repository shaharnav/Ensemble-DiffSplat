import json
import glob
import os
from collections import Counter

def summarize_all_results():
    # Gather all results JSONs (handles generic results.json or specific ones like results/payload_results.json)
    json_files = glob.glob("**/*results.json", recursive=True)
    if os.path.exists("results.json") and "results.json" not in json_files:
        json_files.append("results.json")
        
    if not json_files:
        print("No results.json files found.")
        return

    os.makedirs("result_summaries", exist_ok=True)

    for file_path in json_files:
        print(f"Processing {file_path}...")
        
        # Extract a protein name (e.g., "3PTB_results.json" -> "3PTB")
        basename = os.path.basename(file_path)
        protein_name = basename.replace("_results.json", "").replace("results.json", "baseline").replace(".json", "")
        if not protein_name:
            protein_name = "protein"

        # 1. Load the JSON data
        try:
            with open(file_path, "r") as f:
                data = json.load(f)
        except Exception as e:
            print(f"Failed to process {file_path}: {e}")
            continue

        if not data or not isinstance(data, list):
            print(f"No valid list data in {file_path}. Skipping.")
            continue

        # Prepare output file
        out_path = os.path.join("result_summaries", f"{protein_name}_ligand_rank.txt")
        
        with open(out_path, "w") as out:
            # 2. The ConforMix MVP Analysis
            out.write("--- Winning Conformations Breakdown ---\n")
            conformations = [row.get('winning_conformation') for row in data if row.get('winning_conformation')]
            if conformations:
                counts = Counter(conformations)
                for conf, count in counts.most_common():
                    out.write(f"{conf:<30} {count}\n")
            else:
                out.write("No 'winning_conformation' data available.\n")
            out.write("\n\n")

            # 3. Define the "Goldilocks" Thresholds
            # Strong affinity, very easy to synthesize, and highly drug-like
            elite_candidates = []
            lowest_affinity = None
            
            for row in data:
                aff = row.get('true_affinity')
                sa = row.get('sa_score')
                qed = row.get('qed')

                # Track lowest affinity for debugging context if no candidates are found
                if aff is not None:
                    if lowest_affinity is None or aff < lowest_affinity:
                        lowest_affinity = aff

                # Safely check conditions
                passes_aff = (aff is not None and aff <= -7.0)
                passes_sa = (sa is not None and sa <= 4.0)
                passes_qed = (qed is not None and qed >= 0.5)

                # Check if elements are missing entirely vs failing threshold
                # If SA or QED are missing, we still filter on true_affinity as minimal fallback
                has_filter_data = ('sa_score' in row and 'qed' in row)
                
                if has_filter_data:
                    is_elite = passes_aff and passes_sa and passes_qed
                else:
                    is_elite = passes_aff

                if is_elite:
                    elite_candidates.append(row)

            # Output elite candidates
            out.write(f"--- Top {len(elite_candidates)} Elite Candidates ---\n")
            
            if elite_candidates:
                # Sort the elite candidates by affinity (lowest first)
                elite_candidates.sort(key=lambda r: (r.get('true_affinity') is None, r.get('true_affinity')))

                # Define the columns to print
                header = f"{'Idx':<4} {'SMILES':<60} {'Affinity':<10} {'QED':<6} {'SA Score':<8}"
                out.write(header + "\n")
                out.write("-" * len(header) + "\n")
                
                for row in elite_candidates:
                    idx = row.get('original_index', 'N/A')
                    smi = row.get('smiles', '')
                    if len(smi) > 58:
                        smi = smi[:55] + "..."
                    
                    aff = f"{row.get('true_affinity'):.2f}" if row.get('true_affinity') is not None else "N/A"
                    qed = f"{row.get('qed'):.3f}" if row.get('qed') is not None else "N/A"
                    sa = f"{row.get('sa_score'):.3f}" if row.get('sa_score') is not None else "N/A"
                    
                    out.write(f"{idx:<4} {smi:<60} {aff:<10} {qed:<6} {sa:<8}\n")
            else:
                out.write("No elite candidates found matching the affinity/QED/SA_Score filters.\n")
                if lowest_affinity is not None:
                    out.write(f"\n[Info] Lowest observed affinity in dataset was {lowest_affinity:.2f} kcal/mol\n")
                    
            out.write("\n")
            print(f"Saved ranking to {out_path}")

if __name__ == "__main__":
    summarize_all_results()