#!/usr/bin/env python3
"""
Batch runner for the protein_mn_ligand Chai-1 pipeline.

Mirrors R/run_batch.py exactly; iterates all *.fasta files inside a protein
folder, runs ULTRA-HIGH inference on each, picks the best structure, and saves
the final CIF/PDB back into that folder.

Usage:
  python run_batch.py --protein protein_mn_ligand
  python run_batch.py --protein all
"""
from imports import (
    os,
    sys,
    argparse,
    time,
    tqdm,
    Path,
    parse_fasta_stem_to_dirname,
)
from run_prediction import run_ultrahigh
from process_best import process_output_dir


def run_prediction(fasta_path: Path, output_root: Path, engine: str) -> Path:
    return run_ultrahigh(str(fasta_path), str(output_root), engine)


def process_output(output_dir: Path, base_name: str, final_dir: Path):
    return process_output_dir(str(output_dir), base_name, str(final_dir))


def main():
    parser = argparse.ArgumentParser(description="Batch Chai-1 predictions for protein_mn_ligand")
    parser.add_argument("--protein",   required=True,
                        help="Protein folder name (e.g. protein_mn_ligand) or 'all'")
    parser.add_argument("--workspace", required=False, default=str(Path.cwd()),
                        help="Workspace root containing the protein folders")
    parser.add_argument("--engine",    required=False, choices=["chai1", "chai2"],
                        default="chai1", help="Engine to use (default: chai1)")
    args = parser.parse_args()

    workspace = Path(args.workspace)

    def process_protein_folder(protein_dir: Path):
        fasta_files = sorted(protein_dir.glob("*.fasta"))
        if not fasta_files:
            print(f"No FASTA files found in {protein_dir}")
            return

        print("=" * 70)
        print(f"Running predictions for: {protein_dir.name}  [{args.engine}]")
        print(f"Found {len(fasta_files)} FASTA file(s)")
        print("=" * 70)

        outputs_root = protein_dir / args.engine
        os.makedirs(outputs_root, exist_ok=True)

        for i, fasta in enumerate(tqdm(fasta_files, desc=protein_dir.name, unit="fasta"), 1):
            base_name = parse_fasta_stem_to_dirname(fasta.stem)
            tqdm.write(f"[{i}/{len(fasta_files)}] Predicting: {base_name}")

            t0 = time.time()
            out_dir = run_prediction(fasta, outputs_root, args.engine)
            tqdm.write(f"  ✓ Prediction done in {time.time()-t0:.1f}s  →  {out_dir}")

            tqdm.write(f"[{i}/{len(fasta_files)}] Processing: {base_name}")
            t0 = time.time()
            process_output(out_dir, base_name, protein_dir)
            tqdm.write(f"  ✓ Processing done in {time.time()-t0:.1f}s")

        tqdm.write(f"✓ Completed {protein_dir.name}")

    if args.protein.lower() == "all":
        candidates = [c for c in sorted(workspace.iterdir())
                      if c.is_dir() and any(c.glob("*.fasta"))]
        if not candidates:
            print(f"No protein folders with FASTA files found under {workspace}")
            sys.exit(1)

        print("Folders to process:")
        for c in candidates:
            print(f"  - {c.name}")

        for i, folder in enumerate(tqdm(candidates, desc="All folders", unit="folder"), 1):
            tqdm.write(f"[{i}/{len(candidates)}] {folder.name}")
            process_protein_folder(folder)
            engine_folder = folder / args.engine
            if engine_folder.exists():
                import shutil
                shutil.rmtree(engine_folder)
                tqdm.write(f"  ✓ Cleaned up {engine_folder}")

        tqdm.write(f"\n✓ All {len(candidates)} folder(s) processed.")
    else:
        protein_dir = workspace / args.protein
        if not protein_dir.exists():
            print(f"Error: Folder not found: {protein_dir}")
            sys.exit(1)
        process_protein_folder(protein_dir)
        engine_folder = protein_dir / args.engine
        if engine_folder.exists():
            import shutil
            shutil.rmtree(engine_folder)
            print(f"Cleaned up {engine_folder}")
        print("\nDone.")


if __name__ == "__main__":
    main()
