# Ensemble-DiffSplat

Ensemble-DiffSplat is a high-throughput computational biology pipeline integrating Generative Structure-Based Drug Design (DiffSBDD) with dynamic protein conformational breathing (ConforMix & Boltz). It systematically generates *de novo* ligands and evaluates them against an induced-fit structural ensemble via AutoDock Vina—saving robust cross-docking matrices for dynamic 4D rendering via Gaussian Splatting interfaces.

## Usage: $M \times N$ Ensemble Pipeline

### 1. Cloud Generation (Colab)
To generate your candidate libraries and protein structural states, open the `diffsbdd_generation.ipynb` notebook in a Google Colab GPU-runtime instance.

1. **Configure Parameters (Cell 1):** Set your target `PDB_ID`, select your critical `subset_residues`, and define `N_SAMPLES` to control how many *de novo* molecular geometries DiffSBDD synthesizes.
2. **Tune ConforMix Flexibility:** Scroll down to the `run_conformixrmsd_boltz` cell. 
   - `--num-twist-targets`: Controls the quantity of unique structural variations (e.g., set to 5 for a robust ensemble).
   - `--twist-target-stop`: Controls the maximal root-mean-square deviation (RMSD) opening limit (highly recommended to stay at `2.0` Å to capture natural "loop breathing" without breaking the pocket).
3. **Execute:** Run all cells top-to-bottom. The pipeline will pull crystal structures, generate the matrices, structurally align all variants to the native active-site origin, and package the payload.
4. **Download Payload:** Once execution completes, download `ensemble_payload.zip` locally from the Colab filesystem directory.

### 2. Local Auditor Docking
Your local machine manages the thermodynamic evaluation natively using AutoDock Vina, fully backing up results to disk.

1. **Stage Payload:** Place `ensemble_payload.zip` directly into the `ensemble_input/` directory of this local repository.
2. **Execute Auditor:** Run the following command exactly:
   ```bash
   ./venv/bin/python ensemble_auditor.py --payload ensemble_input/ensemble_payload.zip
   ```
   *The auditor securely checks Vina cache points, auto-constructs ID trackers for SMILES linkage, and parses all states.*

### 3. Interpreting Results
When computation halts, the original payload safely relocates into the timestamp-sealed `archive/` folder. Analyzing the outputs is trivial:

- **`results.csv`**: Contains your full combinatorial matrix. Drag this into any DataFrame to track binding variance or calculate aggregate stability across the protein's conformational sweep.
- **The CLI Rankings:**
  ```text
  --- Top candidates by True Affinity ---
  Rank   ID              Affinity   Baseline   H-Bonds    Winning Conformation     
  ---------------------------------------------------------------------------------
  1      Cmpd-0083       -4.94      -2.85      4          conformix_var_0.pdb     
  ```
  - **Affinity**: The "induced-fit" energy peak achieved against the flexible states.
  - **Baseline**: What that identical drug scored against the rigid immobile crystal default. Use the differential here to prove binding plasticity!
  - **ID**: Look up `Cmpd-0083` strictly within your `results/payload_unpacked/valid_trajectories/` directory. Fetching `mol_0083.xyz` lets you pipe the perfect generated coordinates into your downstream visual render engines linearly with no ambiguity.
