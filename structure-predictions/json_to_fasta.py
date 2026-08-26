#!/usr/bin/env python3
"""
json_to_fasta.py
----------------
Converts a smuprot_N.json file into a Chai-1-compatible FASTA file.

The SMILES for each ligand is written EXACTLY as it appears in the JSON —
including the Mn metal centre.  No atom stripping or charge modification is
performed.

Why the Mn is kept intact
--------------------------
Chai-1 uses RDKit internally to build an initial 3-D conformer of the ligand
before the neural network runs.  RDKit's UFF force-field typer does not have
parameters for transition-metal atom types such as Mn5+2, so it prints a
warning:

    UFFTYPER: Unrecognised atom type: Mn5+2 (12)

This is a WARNING, not an error.  RDKit first builds the conformer with its
distance-geometry method (ETKDG), which does NOT require force-field atom
types and succeeds for metal-containing molecules.  The UFF step is an
optional post-optimisation that fails silently and is skipped.  The Mn atom
and all its bonds are preserved in the conformer that Chai-1 actually uses.

The warning is suppressed here with RDKit's own logger so it does not alarm
the user during batch runs, but it can be re-enabled with --verbose.

Usage
-----
    # Convert one JSON, writing into its default output folder
    python json_to_fasta.py --json smuprot_0.json

    # Convert all three at once
    python json_to_fasta.py --json smuprot_0.json smuprot_1.json smuprot_2.json

    # Specify a custom output folder name
    python json_to_fasta.py --json smuprot_1.json --outdir my_folder

    # Show the UFF warning (useful for debugging)
    python json_to_fasta.py --json smuprot_0.json --verbose
"""

import argparse
import json
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def suppress_rdkit_warnings():
    """Silence RDKit's UFF atom-type warning for metal centres."""
    try:
        from rdkit import RDLogger
        RDLogger.DisableLog('rdApp.warning')
        RDLogger.DisableLog('rdApp.error')
    except ImportError:
        pass   # RDKit not installed; warnings will just appear normally


def load_entry(json_path: Path) -> dict:
    """Return the first (and typically only) entry from a smuprot JSON."""
    with open(json_path) as f:
        data = json.load(f)
    return data[0] if isinstance(data, list) else data


def entry_to_fasta_lines(entry: dict) -> tuple[str, str, str]:
    """
    Extract protein sequence and ligand SMILES from a JSON entry.

    Returns
    -------
    name : str
        The entry's 'name' field.
    protein_seq : str
        Single-letter amino-acid sequence.
    ligand_smiles : str
        SMILES string exactly as stored in the JSON (Mn and all charges intact).
    """
    name = entry['name']
    protein_seq   = None
    ligand_smiles = None

    for seq in entry['sequences']:
        if 'proteinChain' in seq:
            protein_seq = seq['proteinChain']['sequence']
        elif 'ligand' in seq:
            ligand_smiles = seq['ligand']['smiles']

    if protein_seq is None:
        raise ValueError(f"No 'proteinChain' key found in entry '{name}'.")
    if ligand_smiles is None:
        raise ValueError(f"No 'ligand' key found in entry '{name}'.")

    return name, protein_seq, ligand_smiles


def validate_smiles(smiles: str, verbose: bool = False) -> bool:
    """
    Use RDKit to check that the SMILES is parseable.

    If RDKit is not installed the check is skipped and True is returned.
    The UFF warning is suppressed unless verbose=True.
    """
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            print(f"  WARNING: RDKit could not parse SMILES: {smiles[:80]}...")
            return False
        if verbose:
            metal_syms = [a.GetSymbol() for a in mol.GetAtoms()
                          if a.GetAtomicNum() not in (1, 2, 5, 6, 7, 8, 9, 15, 16, 17, 35, 53)]
            if metal_syms:
                print(f"  INFO: Metal atom(s) detected in SMILES: {metal_syms}")
                print("        The UFF typer will warn about these — this is expected and "
                      "non-fatal.\n        Chai-1 uses ETKDG (not UFF) for the initial conformer.")
        return True
    except ImportError:
        return True   # RDKit absent — skip validation


def write_fasta(name: str, protein_seq: str, ligand_smiles: str,
                json_stem: str, out_dir: Path) -> Path:
    """
    Write a Chai-1-format FASTA file.

    File layout
    -----------
    >protein|name=<name>_protein
    <amino-acid sequence>
    >ligand|name=<name>_ligand
    <SMILES>

    The file is named <json_stem>.fasta and placed in out_dir.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    fasta_path = out_dir / f"{json_stem}.fasta"

    lines = [
        f">protein|name={name}_protein",
        protein_seq,
        f">ligand|name={name}_ligand",
        ligand_smiles,
    ]
    fasta_path.write_text('\n'.join(lines) + '\n')
    return fasta_path


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def convert(json_path: Path, out_dir: Path | None, verbose: bool) -> Path:
    """Convert one JSON file to a FASTA. Returns the path to the FASTA."""
    if not verbose:
        suppress_rdkit_warnings()

    entry = load_entry(json_path)
    name, protein_seq, ligand_smiles = entry_to_fasta_lines(entry)

    # Default output folder: same stem as JSON file (e.g. smuprot_0/ for smuprot_0.json)
    if out_dir is None:
        out_dir = json_path.parent / json_path.stem

    print(f"  JSON   : {json_path.name}")
    print(f"  Name   : {name}")
    print(f"  Protein: {len(protein_seq)} residues")
    print(f"  Ligand : {ligand_smiles[:90]}{'...' if len(ligand_smiles) > 90 else ''}")

    validate_smiles(ligand_smiles, verbose=verbose)

    fasta_path = write_fasta(name, protein_seq, ligand_smiles,
                             json_path.stem, out_dir)
    print(f"  FASTA  : {fasta_path}")
    return fasta_path


def main():
    parser = argparse.ArgumentParser(
        description="Convert smuprot JSON file(s) to Chai-1 FASTA format.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        '--json', nargs='+', required=True,
        help='One or more smuprot_N.json files to convert.',
    )
    parser.add_argument(
        '--outdir', default=None,
        help='Output directory. Only valid when converting a single JSON. '
             'Defaults to a folder named after the JSON stem.',
    )
    parser.add_argument(
        '--verbose', action='store_true',
        help='Show RDKit metal-atom info and UFF warnings.',
    )
    args = parser.parse_args()

    if args.outdir and len(args.json) > 1:
        print("Error: --outdir can only be used with a single --json file.")
        sys.exit(1)

    out_dir = Path(args.outdir) if args.outdir else None

    for json_str in args.json:
        json_path = Path(json_str)
        if not json_path.exists():
            print(f"Error: {json_path} not found — skipping.")
            continue
        print(f"\nConverting {json_path.name} ...")
        fasta = convert(json_path, out_dir, args.verbose)
        print(f"  Done → {fasta}\n")


if __name__ == '__main__':
    main()
