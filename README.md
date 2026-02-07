# Enzyme-Substrate Simulator

A local enzyme-substrate docking simulator with a Swiss-style React frontend and Vina/AlphaFold backend.

## Features
- **AlphaFold Integration**: Fetches protein structures automatically.
- **Vina Docking**: Runs local molecular docking simulations using AutoDock Vina.
- **Metal Awareness**: Detecting and centering on catalytic metals (Zn, Mg, etc.) with improved exhaustiveness.
- **Analysis**: Calculates binding affinity and detects hydrogen bonds/metal coordination.

## Usage
1. Install Python dependencies: `pip install -r requirements.txt` (Create requirements.txt if needed)
2. Run: `python start.py`
3. Open `http://localhost:5001`
