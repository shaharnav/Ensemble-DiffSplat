# Enzyme-Substrate Simulator

A local enzyme-substrate docking simulator with a Swiss-style React frontend and Vina/AlphaFold backend.

## Features
- **AlphaFold Integration**: Fetches protein structures automatically.
- **Vina Docking**: Runs local molecular docking simulations using AutoDock Vina.
- **Metal Awareness**: Detecting and centering on catalytic metals (Zn, Mg, etc.) with improved exhaustiveness.
- **Analysis**: Calculates binding affinity and detects hydrogen bonds/metal coordination.

## Usage
1. Create and activate a virtual environment:
   ```bash
   python3 -m venv venv
   source venv/bin/activate
   ```
2. Install Python dependencies:
   ```bash
   pip install -r requirements.txt
   ```
3. Run the application:
   ```bash
   python start.py
   ```
4. Open `http://localhost:5001`
