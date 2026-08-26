#!/usr/bin/env python3
"""
Process Chai-1 prediction results for the protein_mn_ligand system.

Mirrors R/process_best.py.  The only substantive difference is the ligand
relabelling: the original relabelled chain B → PMP and chain C → TLN because
those proteins always had PMP as a co-factor.  Here the single ligand is the
Mn-porphyrin (MNP), so chain B is relabelled MNP instead.
"""
from imports import (
    os,
    re,
    sys,
    glob,
    argparse,
    json,
    shutil,
    np,
    Path,
)


def find_latest_output_dir(pattern="*"):
    dirs = sorted(glob.glob(pattern))
    if not dirs:
        print(f"Error: No directories matching '{pattern}' found!")
        return None
    latest = dirs[-1]
    print(f"Using latest output directory: {latest}")
    return latest


def load_score_from_npz(npz_file):
    try:
        data = np.load(npz_file)
        for key in ['confidence', 'score', 'plddt', 'ptm', 'overall_confidence']:
            if key in data:
                score = data[key]
                if isinstance(score, np.ndarray):
                    return float(score.item()) if score.size == 1 else float(score.mean())
                return float(score)
        for key in data.files:
            value = data[key]
            if isinstance(value, (int, float, np.number)):
                return float(value)
            elif isinstance(value, np.ndarray) and value.size == 1:
                return float(value.item())
    except Exception as e:
        print(f"Error loading {npz_file}: {e}")
    return None


def find_best_structure(output_dir):
    output_path = Path(output_dir)
    if not output_path.exists():
        print(f"Error: '{output_dir}' not found!")
        return None

    cif_files = list(output_path.glob("pred.model_idx_*.cif"))
    if not cif_files:
        print("Error: No CIF files found!")
        return None

    print(f"Found {len(cif_files)} structure files")
    best_score = -float('inf')
    best_file = None
    all_scores = []

    for cif_file in cif_files:
        match = re.search(r'model_idx_(\d+)', cif_file.name)
        if match:
            idx = match.group(1)
            score_file = cif_file.parent / f"scores.model_idx_{idx}.npz"
            if score_file.exists():
                score = load_score_from_npz(score_file)
                if score is not None:
                    all_scores.append((cif_file, score))
                    print(f"  {cif_file.name}: score = {score:.4f}")
                    if score > best_score:
                        best_score = score
                        best_file = cif_file

    if best_file:
        print(f"\nBest: {best_file.name}  (score={best_score:.4f})")
        if len(all_scores) > 1:
            scores_only = [s[1] for s in all_scores]
            print(f"  Range {min(scores_only):.4f} – {max(scores_only):.4f}  "
                  f"mean {np.mean(scores_only):.4f}  std {np.std(scores_only):.4f}")
        return best_file, best_score, all_scores
    else:
        print("Error: Could not determine best structure!")
        return None


def convert_cif_to_pdb(cif_file, output_name=None):
    """Convert CIF → PDB via biotite (preferred) or gemmi."""
    out_name = output_name or (cif_file.stem + ".pdb")
    out_path = cif_file.parent / out_name

    # biotite
    try:
        import biotite.structure.io.pdbx as pdbx
        import biotite.structure.io.pdb as pdb
        cif_obj = pdbx.CIFFile.read(str(cif_file))
        atoms = pdbx.get_structure(cif_obj, model=1)
        pdb_file = pdb.PDBFile()
        pdb.set_structure(pdb_file, atoms)
        pdb_file.write(str(out_path))
        print(f"CIF → PDB (biotite): {out_path}")
        return out_path
    except Exception as e:
        print(f"biotite failed: {e}")

    # gemmi fallback
    try:
        import gemmi
        doc = gemmi.cif.read(str(cif_file))
        st = gemmi.make_structure_from_block(doc.sole_block())
        st.write_pdb(str(out_path))
        print(f"CIF → PDB (gemmi): {out_path}")
        return out_path
    except Exception as e:
        print(f"gemmi failed: {e}")

    return None


def relabel_ligands_in_pdb(pdb_path, chain_to_resname=None):
    """Relabel HETATM residues named 'LIG' by chain ID.

    Default for protein_mn_ligand:
      chain B → MNP  (the Mn-porphyrin ligand)
    Adjust chain_to_resname if Chai-1 places the ligand on a different chain.
    """
    if chain_to_resname is None:
        chain_to_resname = {"B": "MNP"}

    try:
        with open(pdb_path, 'r') as f:
            lines = f.readlines()

        changed = 0
        new_lines = []
        for line in lines:
            if (line.startswith('HETATM') or line.startswith('ATOM  ')) and len(line) >= 22:
                chain_id = line[21]
                if chain_id in chain_to_resname and line[17:20].strip() == 'LIG':
                    new_res = chain_to_resname[chain_id][:3].rjust(3)
                    line = line[:17] + new_res + line[20:]
                    changed += 1
            new_lines.append(line)

        with open(pdb_path, 'w') as f:
            f.writelines(new_lines)

        print(f"Relabelled {changed} atoms in {Path(pdb_path).name} "
              f"(chains: {list(chain_to_resname.keys())})")
    except Exception as e:
        print(f"Warning: relabelling failed: {e}")


def derive_complex_name(output_dir: str, base_name: str | None = None) -> str:
    if base_name:
        return base_name
    out = Path(output_dir)
    meta_path = out / "metadata.json"
    if meta_path.exists():
        try:
            with open(meta_path) as f:
                meta = json.load(f)
            if isinstance(meta, dict):
                return str(meta.get("directory_name") or meta.get("complex_name") or "")
        except Exception:
            pass
    try:
        cif_files = sorted(out.glob("pred.model_idx_*.cif"))
        if cif_files:
            return cif_files[0].stem
    except Exception:
        pass
    return Path(output_dir).name


def process_output_dir(output_dir: str, base_name: str | None = None, final_dir: str | None = None):
    """Pick best CIF, convert to PDB, write summary. Returns (cif_path, pdb_path, score)."""
    result = find_best_structure(output_dir)
    if result is None:
        raise RuntimeError("Failed to find best structure")
    best_file, best_score, all_scores = result

    complex_name = derive_complex_name(output_dir, base_name)
    out_dir  = Path(output_dir)
    dest_dir = Path(final_dir) if final_dir else out_dir
    os.makedirs(dest_dir, exist_ok=True)

    new_cif_path = dest_dir / f"{complex_name}.cif"
    shutil.copy2(best_file, new_cif_path)
    print(f"\nSaved CIF: {new_cif_path}")

    pdb_path = convert_cif_to_pdb(new_cif_path, f"{complex_name}.pdb")
    if pdb_path:
        # Relabel the Mn-porphyrin ligand (chain B → MNP by default)
        relabel_ligands_in_pdb(pdb_path, {"B": "MNP"})

    # Summary
    timestamp = Path(output_dir).name.split('_')[-1]
    summary_file = Path(output_dir) / f"summary_{timestamp}.txt"
    with open(summary_file, 'w') as f:
        f.write(f"Chai-1 Prediction Summary — protein_mn_ligand\n")
        f.write("=" * 50 + "\n")
        f.write(f"Complex      : {complex_name}\n")
        f.write(f"Engine config: 500 steps, 5 samples\n")
        f.write(f"Best CIF     : {Path(best_file).name}\n")
        f.write(f"Best score   : {best_score:.4f}\n")
        f.write(f"Output CIF   : {new_cif_path.name}\n")
        f.write(f"Output PDB   : {pdb_path.name if pdb_path else 'conversion failed'}\n")
        if len(all_scores) > 1:
            scores_only = [s[1] for s in all_scores]
            f.write(f"\nScore stats ({len(scores_only)} samples):\n")
            f.write(f"  mean {np.mean(scores_only):.4f}  std {np.std(scores_only):.4f}  "
                    f"min {np.min(scores_only):.4f}  max {np.max(scores_only):.4f}\n")

    return new_cif_path, (pdb_path if pdb_path else None), best_score


def cli_main():
    parser = argparse.ArgumentParser(description="Process Chai-1 output for protein_mn_ligand")
    parser.add_argument("--output-dir", required=False)
    parser.add_argument("--base-name",  required=False)
    parser.add_argument("--final-dir",  required=False)
    args = parser.parse_args()

    output_dir = args.output_dir or find_latest_output_dir()
    if not output_dir:
        print("Error: No output directory found.")
        sys.exit(1)

    final_dir = args.final_dir or str(Path(output_dir).parent)

    try:
        cif_path, pdb_path, best_score = process_output_dir(output_dir, args.base_name, final_dir)
    except Exception as e:
        import traceback
        print(f"Processing failed: {e}")
        traceback.print_exc()
        sys.exit(1)

    print(f"\nFINAL_CIF={cif_path}")
    print(f"FINAL_PDB={pdb_path or ''}")
    print(f"BEST_SCORE={best_score:.4f}")


if __name__ == "__main__":
    cli_main()
