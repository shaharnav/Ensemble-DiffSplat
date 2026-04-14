"""
ensemble_auditor.py — De Novo Drug Pipeline Validation Engine
=============================================================
Implements an N×M ensemble cross-docking matrix:
  - N  SMILES drug candidates (--smiles)
  - M  protein conformational variants from ConforMix (--receptors)

For each SMILES, the ligand is cross-docked against every .pdb receptor.
The lowest (most negative) affinity across all M conformations is treated
as the "True Affinity" (induced-fit scoring proxy).

Supports two input modes:
  1. Classic:  --smiles FILE --receptors DIR
  2. Payload:  --payload ZIP  (auto-extracts SMILES from SDF, pocket from PDB,
               and docking box from metadata.json)

Outputs a ranked results.json with per-candidate fields:
  smiles, true_affinity, winning_conformation, h_bond_count

Active-site targeting is delegated entirely to docking_engine.run_docking,
which applies the catalytic-triad / metal / cavity fallback cascade defined
in get_center_and_size — guaranteeing Vina never defaults to blind docking.
When running in payload mode, the explicit pocket center from metadata.json
is forwarded to the docking engine, bypassing auto-detection.
"""

import argparse
import glob
import json
import logging
import os
import shutil
import sys
import tempfile
import zipfile

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

    # ── Payload mode (new) ────────────────────────────────────────────────
    parser.add_argument(
        "--payload",
        default=None,
        metavar="PAYLOAD_ZIP",
        help=(
            "Path to an ensemble_payload.zip produced by the DiffSBDD Colab "
            "notebook.  Contains valid_candidates.sdf, pocket.pdb, "
            "valid_trajectories/, and metadata.json.  When provided, "
            "--smiles and --receptors are not required."
        ),
    )

    # ── Classic mode ──────────────────────────────────────────────────────
    parser.add_argument(
        "--smiles",
        default=None,
        metavar="SMILES_FILE",
        help="Path to a text file containing N SMILES strings, one per line.",
    )
    parser.add_argument(
        "--receptors",
        default=None,
        metavar="RECEPTORS_DIR",
        help="Directory containing M ConforMix .pdb receptor conformations.",
    )

    # ── Shared options ────────────────────────────────────────────────────
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
    args = parser.parse_args()

    # Validate: either --payload OR (--smiles AND --receptors) must be given
    if args.payload is None and (args.smiles is None or args.receptors is None):
        parser.error(
            "Either --payload OR both --smiles and --receptors must be provided."
        )

    return args


# ---------------------------------------------------------------------------
# Payload unpacking
# ---------------------------------------------------------------------------
def unpack_payload(zip_path: str, dest_dir: str) -> dict:
    """
    Unzip *zip_path* into *dest_dir* and return a dict with resolved paths
    and metadata.

    Expected zip layout:
        metadata.json
        pocket.pdb
        valid_candidates.sdf
        valid_trajectories/
    """
    if not os.path.isfile(zip_path):
        logger.error(f"Payload zip not found: {zip_path}")
        sys.exit(1)

    logger.info(f"Unpacking payload: {zip_path} → {dest_dir}")
    with zipfile.ZipFile(zip_path, "r") as zf:
        zf.extractall(dest_dir)

    # ── Read metadata ─────────────────────────────────────────────────────
    meta_path = os.path.join(dest_dir, "metadata.json")
    if not os.path.isfile(meta_path):
        logger.error("metadata.json not found in payload.")
        sys.exit(1)

    with open(meta_path) as fh:
        metadata = json.load(fh)

    logger.info(
        f"  Payload metadata: pdb_id={metadata.get('pdb_id')}, "
        f"center={metadata.get('pocket_center')}, "
        f"radius={metadata.get('pocket_radius')}, "
        f"candidates={metadata.get('n_valid_candidates')}"
    )

    # ── Locate expected files ─────────────────────────────────────────────
    pocket_pdb = os.path.join(dest_dir, "pocket.pdb")
    sdf_path = os.path.join(dest_dir, "valid_candidates.sdf")
    traj_dir = os.path.join(dest_dir, "valid_trajectories")

    for label, path in [("pocket.pdb", pocket_pdb),
                        ("valid_candidates.sdf", sdf_path)]:
        if not os.path.isfile(path):
            logger.error(f"{label} not found in unpacked payload at {path}")
            sys.exit(1)

    return {
        "metadata": metadata,
        "pocket_pdb": pocket_pdb,
        "sdf_path": sdf_path,
        "traj_dir": traj_dir if os.path.isdir(traj_dir) else None,
    }


def extract_smiles_from_sdf(sdf_path: str) -> list[dict]:
    """
    Parse every molecule from an SDF file and return a list of dicts, each
    containing 'smiles' and any SD properties (QED, SA_Score, OriginalIndex).

    Uses RDKit; molecules that fail to parse are skipped with a warning.
    """
    from rdkit import Chem

    suppl = Chem.SDMolSupplier(sdf_path)
    candidates = []

    for idx, mol in enumerate(suppl):
        if mol is None:
            logger.warning(f"  Skipping molecule {idx} — RDKit parse failure.")
            continue

        smiles = Chem.MolToSmiles(mol)
        if not smiles:
            logger.warning(f"  Skipping molecule {idx} — empty SMILES.")
            continue

        entry = {"smiles": smiles}

        # Carry forward SDF properties
        for prop in mol.GetPropsAsDict():
            entry[prop] = mol.GetPropsAsDict()[prop]

        candidates.append(entry)

    logger.info(f"  Extracted {len(candidates)} SMILES from {sdf_path}")
    return candidates


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
    center_coords=None,
    box_size=None,
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
        center_coords=center_coords,
        box_size=box_size,
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
    center_coords=None,
    box_size=None,
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
                center_coords=center_coords,
                box_size=box_size,
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
# Payload mode orchestrator
# ---------------------------------------------------------------------------
def run_payload_mode(
    payload_zip: str,
    output_path: str,
    exhaustiveness: int,
    target_residue: str | None,
) -> None:
    """
    End-to-end pipeline: unzip → extract SMILES → dock → rank → save.
    """
    # ── 1. Unpack ─────────────────────────────────────────────────────────
    unpack_dir = os.path.join("results", "payload_unpacked")
    os.makedirs(unpack_dir, exist_ok=True)
    payload = unpack_payload(payload_zip, unpack_dir)

    metadata = payload["metadata"]
    pocket_pdb = payload["pocket_pdb"]
    sdf_path = payload["sdf_path"]

    # ── 2. Extract SMILES from SDF ────────────────────────────────────────
    sdf_entries = extract_smiles_from_sdf(sdf_path)
    if not sdf_entries:
        logger.error("No valid molecules found in the payload SDF.")
        sys.exit(1)

    smiles_list = [e["smiles"] for e in sdf_entries]

    # ── 3. Resolve docking box from metadata ──────────────────────────────
    pocket_center = metadata.get("pocket_center")
    pocket_radius = metadata.get("pocket_radius", 10.0)

    if pocket_center is not None:
        center_coords = tuple(pocket_center)
        # Box size: use 2× radius as the search box dimension (capped at 25Å)
        box_dim = min(pocket_radius * 2.0, 25.0)
        box_size = [box_dim, box_dim, box_dim]
        logger.info(
            f"Using metadata pocket center {center_coords}, "
            f"box {box_size} (radius {pocket_radius} Å)"
        )
    else:
        center_coords = None
        box_size = None
        logger.warning(
            "No pocket_center in metadata — falling back to auto-detection."
        )

    # ── 4. Receptors: single pocket.pdb for now ───────────────────────────
    receptor_paths = [pocket_pdb]
    logger.info(f"Receptor: {pocket_pdb}")

    # ── 5. Run the ensemble docking matrix ────────────────────────────────
    work_dir = os.path.join("results", "ensemble_audit")
    os.makedirs(work_dir, exist_ok=True)
    logger.info(f"Intermediate docking files → '{work_dir}'")

    ranked = run_ensemble(
        smiles_list=smiles_list,
        receptor_paths=receptor_paths,
        work_dir=work_dir,
        exhaustiveness=exhaustiveness,
        target_residue=target_residue,
        center_coords=center_coords,
        box_size=box_size,
    )

    # ── 6. Enrich results with SDF properties (QED, SA) ───────────────────
    # Build a SMILES→SDF-entry lookup for merging QED/SA scores
    sdf_lookup = {e["smiles"]: e for e in sdf_entries}
    for entry in ranked:
        sdf_data = sdf_lookup.get(entry["smiles"], {})
        entry["qed"] = sdf_data.get("QED")
        entry["sa_score"] = sdf_data.get("SA_Score")
        entry["original_index"] = sdf_data.get("OriginalIndex")

    # ── 7. Save results ───────────────────────────────────────────────────
    with open(output_path, "w") as fh:
        json.dump(ranked, fh, indent=2)

    logger.info(f"\n{'='*60}")
    logger.info(
        f"Payload audit complete. {len(ranked)} candidate(s) ranked."
    )
    logger.info(f"Target: {metadata.get('pdb_id', 'unknown')}")
    logger.info(f"Results saved → {os.path.abspath(output_path)}")
    logger.info(f"{'='*60}\n")

    # Pretty-print top-5
    _print_rankings(ranked)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def main() -> None:
    args = parse_args()

    if args.payload is not None:
        # ── Payload mode ──────────────────────────────────────────────────
        run_payload_mode(
            payload_zip=args.payload,
            output_path=args.output,
            exhaustiveness=args.exhaustiveness,
            target_residue=args.target_residue,
        )
    else:
        # ── Classic mode ──────────────────────────────────────────────────
        smiles_list = load_smiles(args.smiles)
        receptor_paths = load_receptors(args.receptors)

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
        with open(args.output, "w") as fh:
            json.dump(ranked, fh, indent=2)

        logger.info(f"\n{'='*60}")
        logger.info(f"Ensemble audit complete. {len(ranked)} candidate(s) ranked.")
        logger.info(f"Results saved → {os.path.abspath(args.output)}")
        logger.info(f"{'='*60}\n")

        _print_rankings(ranked)


def _print_rankings(ranked: list[dict]) -> None:
    """Pretty-print the top-5 ranked candidates to stdout."""
    top = ranked[:5]
    print("\n--- Top candidates by True Affinity ---")

    # Determine which columns to show based on available data
    has_qed = any(e.get("qed") is not None for e in top)

    header = f"{'Rank':<6} {'Affinity (kcal/mol)':<22} {'H-Bonds':<10} {'Winning Conformation':<25}"
    if has_qed:
        header += f" {'QED':<8} {'SA':<8}"
    print(header)
    print("-" * len(header))

    for rank, entry in enumerate(top, start=1):
        aff = f"{entry['true_affinity']:.2f}" if entry["true_affinity"] is not None else "N/A"
        line = (
            f"{rank:<6} {aff:<22} {entry['h_bond_count']:<10} "
            f"{entry['winning_conformation']:<25}"
        )
        if has_qed:
            qed_str = f"{entry['qed']}" if entry.get("qed") is not None else "N/A"
            sa_str = f"{entry['sa_score']}" if entry.get("sa_score") is not None else "N/A"
            line += f" {qed_str:<8} {sa_str:<8}"
        print(line)
    print()


if __name__ == "__main__":
    main()
