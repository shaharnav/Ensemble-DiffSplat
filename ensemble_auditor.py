"""
ensemble_auditor.py — De Novo Drug Pipeline Validation Engine
=============================================================
Implements an N×M ensemble cross-docking matrix:
  - N  SMILES drug candidates (--smiles)
  - M  protein conformational variants from ConforMix (--receptors)

For each SMILES, the ligand is cross-docked against every .pdb receptor.
The lowest (most negative) affinity across all M conformations is treated
as the "True Affinity" (induced-fit scoring proxy).

Outputs a ranked results.json with per-candidate fields:
  smiles, true_affinity, winning_conformation, h_bond_count

Active-site targeting is delegated entirely to docking_engine.run_docking,
which applies the catalytic-triad / metal / cavity fallback cascade defined
in get_center_and_size — guaranteeing Vina never defaults to blind docking.
"""

import argparse
import glob
import json
import logging
import os
import sys
import tempfile

from docking_engine import run_docking
from analyzer import analyze_docking

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="ensemble_auditor",
        description=(
            "N×M ensemble cross-docking auditor for de novo drug candidates. "
            "Ranks every SMILES by its best (most negative) affinity across "
            "all M protein conformational variants."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--smiles",
        required=True,
        metavar="SMILES_FILE",
        help="Path to a text file containing N SMILES strings, one per line.",
    )
    parser.add_argument(
        "--receptors",
        required=True,
        metavar="RECEPTORS_DIR",
        help="Directory containing M ConforMix .pdb receptor conformations.",
    )
    parser.add_argument(
        "--output",
        default="results.json",
        metavar="OUTPUT_JSON",
        help="Path for the ranked JSON output file.",
    )
    parser.add_argument(
        "--exhaustiveness",
        type=int,
        default=16,
        metavar="INT",
        help="AutoDock Vina exhaustiveness (lower = faster, higher = more accurate).",
    )
    parser.add_argument(
        "--target-residue",
        default=None,
        metavar="RESIDUE",
        help=(
            "Optional active-site residue hint forwarded to get_center_and_size "
            "(e.g. 'ASN253'). Falls back through metal → cavity search if omitted."
        ),
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def load_smiles(path: str) -> list[str]:
    """Return non-empty, stripped SMILES lines from *path*."""
    if not os.path.isfile(path):
        logger.error(f"SMILES file not found: {path}")
        sys.exit(1)

    smiles_list: list[str] = []
    with open(path, "r") as fh:
        for lineno, raw in enumerate(fh, start=1):
            smi = raw.strip()
            if smi and not smi.startswith("#"):
                smiles_list.append(smi)
            else:
                logger.debug(f"Skipping blank/comment line {lineno}")

    if not smiles_list:
        logger.error("No valid SMILES found in the input file.")
        sys.exit(1)

    logger.info(f"Loaded {len(smiles_list)} SMILES candidate(s) from '{path}'.")
    return smiles_list


def load_receptors(directory: str) -> list[str]:
    """Return sorted list of .pdb file paths inside *directory*."""
    if not os.path.isdir(directory):
        logger.error(f"Receptors directory not found: {directory}")
        sys.exit(1)

    pdbs = sorted(glob.glob(os.path.join(directory, "*.pdb")))
    if not pdbs:
        logger.error(f"No .pdb files found in '{directory}'.")
        sys.exit(1)

    logger.info(f"Found {len(pdbs)} receptor conformation(s) in '{directory}'.")
    return pdbs


def make_job_name(smiles: str, pdb_path: str, smiles_idx: int) -> str:
    """Build a filesystem-safe job identifier."""
    pdb_stem = os.path.splitext(os.path.basename(pdb_path))[0]
    return f"ens_s{smiles_idx:04d}_{pdb_stem}"


# ---------------------------------------------------------------------------
# Core: single SMILES × single receptor
# ---------------------------------------------------------------------------
def dock_one(
    smiles: str,
    pdb_path: str,
    job_name: str,
    work_dir: str,
    exhaustiveness: int,
    target_residue: str | None,
) -> tuple[float | None, int, str]:
    """
    Dock *smiles* against *pdb_path*.

    Returns
    -------
    (affinity, h_bond_count, docked_pdbqt_path)
        affinity is None on failure.
    """
    result = run_docking(
        pdb_file=pdb_path,
        smiles=smiles,
        output_dir=work_dir,
        job_name=job_name,
        exhaustiveness=exhaustiveness,
        target_residue=target_residue,
    )

    if result is None or result.get("affinity") is None:
        logger.warning(f"  ✗ Docking failed for job '{job_name}'.")
        return None, 0, ""

    affinity: float = result["affinity"]
    docked_pdbqt: str = result.get("docked_file", "")

    # Interaction analysis to retrieve h_bond_count
    h_bond_count = 0
    if docked_pdbqt and os.path.isfile(docked_pdbqt):
        try:
            analysis = analyze_docking(pdb_path, docked_pdbqt)
            h_bond_count = analysis.get("h_bond_count", 0)
        except Exception as exc:
            logger.warning(f"  ⚠ Interaction analysis failed for '{job_name}': {exc}")

    logger.info(
        f"  ✓ Job '{job_name}' → affinity={affinity:.2f} kcal/mol, "
        f"h_bonds={h_bond_count}"
    )
    return affinity, h_bond_count, docked_pdbqt


# ---------------------------------------------------------------------------
# Core: N×M ensemble matrix
# ---------------------------------------------------------------------------
def run_ensemble(
    smiles_list: list[str],
    receptor_paths: list[str],
    work_dir: str,
    exhaustiveness: int,
    target_residue: str | None,
) -> list[dict]:
    """
    Execute the full N×M cross-docking matrix.

    Returns a list of result dicts (one per SMILES), sorted by true_affinity
    ascending (most negative first).
    """
    n_smiles = len(smiles_list)
    m_receptors = len(receptor_paths)
    total_jobs = n_smiles * m_receptors
    logger.info(
        f"Starting {n_smiles}×{m_receptors} ensemble matrix "
        f"({total_jobs} total docking jobs)."
    )

    candidates: list[dict] = []

    for si, smiles in enumerate(smiles_list):
        logger.info(
            f"\n[{si + 1}/{n_smiles}] SMILES: {smiles[:60]}"
            f"{'…' if len(smiles) > 60 else ''}"
        )

        best_affinity: float | None = None
        best_h_bonds: int = 0
        best_conformation: str = ""

        for pdb_path in receptor_paths:
            pdb_basename = os.path.basename(pdb_path)
            job_name = make_job_name(smiles, pdb_path, si)

            logger.info(f"  Docking vs. {pdb_basename} …")
            affinity, h_bonds, docked_file = dock_one(
                smiles=smiles,
                pdb_path=pdb_path,
                job_name=job_name,
                work_dir=work_dir,
                exhaustiveness=exhaustiveness,
                target_residue=target_residue,
            )

            # Keep the most negative (lowest) affinity score — induced-fit winner
            if affinity is not None:
                if best_affinity is None or affinity < best_affinity:
                    best_affinity = affinity
                    best_h_bonds = h_bonds
                    best_conformation = pdb_basename

        if best_affinity is None:
            logger.warning(
                f"  All docking attempts failed for SMILES index {si}. "
                "Candidate will still appear in output with null affinity."
            )

        candidates.append(
            {
                "smiles": smiles,
                "true_affinity": best_affinity,
                "winning_conformation": best_conformation,
                "h_bond_count": best_h_bonds,
            }
        )

    # Global ranking: most negative affinity first; failures sink to the bottom
    candidates.sort(
        key=lambda c: (
            c["true_affinity"] is None,          # None → sort last
            c["true_affinity"] if c["true_affinity"] is not None else 0.0,
        )
    )

    return candidates


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def main() -> None:
    args = parse_args()

    smiles_list = load_smiles(args.smiles)
    receptor_paths = load_receptors(args.receptors)

    # Intermediate docking artefacts go into a dedicated subdirectory so they
    # don't pollute the working directory and can be inspected afterwards.
    work_dir = os.path.join("results", "ensemble_audit")
    os.makedirs(work_dir, exist_ok=True)
    logger.info(f"Intermediate docking files will be written to '{work_dir}'.")

    ranked = run_ensemble(
        smiles_list=smiles_list,
        receptor_paths=receptor_paths,
        work_dir=work_dir,
        exhaustiveness=args.exhaustiveness,
        target_residue=args.target_residue,
    )

    # Serialize results
    output_path = args.output
    with open(output_path, "w") as fh:
        json.dump(ranked, fh, indent=2)

    logger.info(f"\n{'='*60}")
    logger.info(f"Ensemble audit complete. {len(ranked)} candidate(s) ranked.")
    logger.info(f"Results saved → {os.path.abspath(output_path)}")
    logger.info(f"{'='*60}\n")

    # Pretty-print top-5 to stdout
    top = ranked[:5]
    print("\n--- Top candidates by True Affinity ---")
    print(f"{'Rank':<6} {'Affinity (kcal/mol)':<22} {'H-Bonds':<10} {'Winning Conformation'}")
    print("-" * 80)
    for rank, entry in enumerate(top, start=1):
        aff = f"{entry['true_affinity']:.2f}" if entry["true_affinity"] is not None else "N/A"
        print(
            f"{rank:<6} {aff:<22} {entry['h_bond_count']:<10} "
            f"{entry['winning_conformation']}"
        )
    print()


if __name__ == "__main__":
    main()
