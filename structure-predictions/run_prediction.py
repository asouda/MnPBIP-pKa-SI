#!/usr/bin/env python3
"""
ULTRA-HIGH quality Chai-1 prediction script — protein_mn_ligand variant.
500 diffusion steps, 5 samples, maximum quality configuration.

Mirrors R/run_prediction.py with no logic changes; imports from local imports.py.
"""
from imports import (
    os,
    json,
    sys,
    time,
    argparse,
    inspect,
    datetime,
    Path,
    get_chai_run_inference,
    parse_fasta_stem_to_dirname,
)

os.environ['DISABLE_PANDERA_IMPORT_WARNING'] = 'True'


def get_timestamp():
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def detect_fasta_and_complex(explicit_fasta: str | None = None):
    """Detect the FASTA file and complex name.

    If explicit_fasta is provided, use it. Otherwise, pick the first *.fasta in CWD.
    """
    if explicit_fasta is not None:
        fasta_path = Path(explicit_fasta)
        if not fasta_path.exists():
            raise FileNotFoundError(f"FASTA file not found: {explicit_fasta}")
    else:
        fasta_files = list(Path.cwd().glob("*.fasta"))
        if not fasta_files:
            raise FileNotFoundError("No fasta files found in current directory")
        fasta_path = fasta_files[0]

    complex_name = fasta_path.stem
    return complex_name, str(fasta_path)


def run_ultrahigh(fasta: str, output_root: str | None = None, engine: str = "chai1") -> Path:
    """Run ULTRA-HIGH quality prediction on a FASTA. Returns the output directory."""
    run_inference = get_chai_run_inference(engine)

    complex_name, fasta_file = detect_fasta_and_complex(fasta)
    if not os.path.exists(fasta_file):
        raise FileNotFoundError(f"FASTA file '{fasta_file}' not found!")

    fasta_input = Path(fasta_file)
    dir_name = parse_fasta_stem_to_dirname(fasta_input.stem)

    base_dir = Path(output_root) if output_root else Path.cwd()
    output_dir = base_dir / dir_name
    os.makedirs(output_dir, exist_ok=True)

    print("=" * 70)
    print(f"ULTRA-HIGH Quality Chai-1 Prediction for {complex_name}")
    print("=" * 70)
    print(f"  Complex  : {complex_name}")
    print(f"  FASTA    : {fasta_file}")
    print(f"  Engine   : {engine}")
    print(f"  Steps    : 500  |  Samples : 5")
    print(f"  Output   : {output_dir}")
    print("-" * 70)

    all_kwargs = {
        "output_dir":          output_dir,
        "num_trunk_recycles":  3,
        "num_diffn_timesteps": 500,
        "num_diffn_samples":   5,
        "use_esm_embeddings":  True,
    }
    try:
        sig = inspect.signature(run_inference)
        valid_kwargs = {k: v for k, v in all_kwargs.items() if k in sig.parameters}
    except Exception:
        valid_kwargs = all_kwargs

    start_time = time.time()
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Starting prediction ...")
    run_inference(fasta_input, **valid_kwargs)

    meta = {
        "complex_name":        complex_name,
        "fasta_file":          str(fasta_input.resolve()),
        "output_dir":          str(output_dir.resolve()),
        "engine":              engine,
        "num_trunk_recycles":  3,
        "num_diffn_timesteps": 500,
        "num_diffn_samples":   5,
        "directory_name":      dir_name,
    }
    try:
        with open(output_dir / "metadata.json", "w") as f:
            json.dump(meta, f, indent=2)
    except Exception:
        pass

    elapsed = time.time() - start_time
    print(f"\nPrediction completed in {int(elapsed//60)}:{int(elapsed%60):02d}")
    print(f"Generated files:")
    for cif in sorted(output_dir.glob("*.cif")):
        print(f"  - {cif.name}")
    print(f"Results saved to: {output_dir}")

    return output_dir


def main():
    parser = argparse.ArgumentParser(description="Run ULTRA-HIGH Chai-1 prediction")
    parser.add_argument("--fasta",  required=False, help="Path to FASTA (defaults to first *.fasta in CWD)")
    parser.add_argument("--outdir", required=False, help="Directory for output (defaults to CWD)")
    parser.add_argument("--engine", required=False, choices=["chai1", "chai2"], default="chai1")
    args = parser.parse_args()

    try:
        output_dir = run_ultrahigh(args.fasta, args.outdir, args.engine)
    except ImportError as e:
        print(f"Error: Chai-1 not properly installed\n  {e}")
        sys.exit(1)
    except Exception as e:
        import traceback
        print(f"Prediction failed: {e}")
        traceback.print_exc()
        sys.exit(1)

    print(f"OUTPUT_DIR={output_dir}")


if __name__ == "__main__":
    main()
