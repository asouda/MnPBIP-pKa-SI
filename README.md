# MnP-BIP pKa Supplementary Information

Supplementary computational materials for a study of the redox-linked pKa shift of a Mn-porphyrin/bis-imidazolylporphyrin (MnP-BIP) complex bound to a 115-residue protein.

The repository covers the two workflows of the calculation pipeline:

1. **`structure-predictions/`** - AI structure prediction (Chai-1) of the protein bound to each of several Mn-porphyrin ligand variants, producing candidate 3-D complex structures.

2. **`pb-electrostatics/`** - continuum electrostatics (APBS) on one such complex, at four combinations of {Mn(II)/Mn(III)} porphyrin with {protonated/deprotonated} benzimidazol phenol (BIP} group, combined in a thermodynamic double-difference cycle to estimate the pKa shift upon the oxidation of Mn(II) to Mn(III).

Each folder has its own detailed `README.md` (script-by-script reference, example results, caveats). This root README gives the overview and ties the two together.

## Repository layout

```
.
├── pyproject.toml                 uv-managed Python environment (pb-electrostatics stack)
├── .python-version                pins Python 3.11
│
├── structure-predictions/         Chai-1 structure prediction of protein + Mn-porphyrin
│   ├── README.md                  -> full script reference, ligand notes, how to run
│   ├── smuprot_0.json … smuprot_4_protonated.json    6 ligand-variant input specs
│   ├── json_to_fasta.py, run_batch.py, run_prediction.py, process_best.py, imports.py
│   ├── run.sh                     HPC/command-line entry point
│   ├── run_protein_mn_ligand.ipynb   Colab notebook version
│   ├── requirements.txt, environment.yml    separate env (chai_lab, GPU) - see below
│   └── smuprot_0/ … smuprot_4_protonated/   predicted output (.fasta/.cif/.pdb) per variant
│
└── pb-electrostatics/              APBS continuum-electrostatics pKa-shift workflow
    ├── README.md                  -> full workflow reference, CLI options, example results
    ├── protein-MnP-BIP-complex.pdb   all-atom complex template (115-residue protein + MNP)
    ├── gaussian/                  4 Gaussian16 MK-charge log files (DFT, one per state)
    ├── apbs/                      4 PQR structures + matching APBS input files
    ├── pka_workflow.py            command-line pipeline: APBS runs -> pKa shift -> PyMOL viz
    ├── pka_workflow.ipynb         same pipeline as an executed, narrated notebook
    └── requirements.txt, environment.yml   alternative (non-uv) environment managers
```

`structure-predictions/` predicts the 3-D structure of the protein bound to the Mn-porphyrin ligand (6 SMILES variants are screened; `smuprot_4_protonated.json` in particular already encodes a protonated peripheral imidazole N as `[NH+]` in its SMILES - the same N-H group whose pKa shift is the subject of `pb-electrostatics/`). The 115-protein-residue + single-MNP-ligand composition of `pb-electrostatics/protein-MnP-BIP-complex.pdb` matches that target exactly, so this file is the geometric template underlying the electrostatics calculations.

**Note:** The intermediate work - the generation of the PQR files with the DFT-derived ligand charges for each of the four states (from the `gaussian/*.log` Merz-Kollman partial charges) - was done outside the workflows described here. 

## Workflow summary

| | Stage | Folder | Tool | Output |
|---|---|---|---|---|
|1| Predict complex structure | `structure-predictions/` | Chai-1 (deep learning) | `smuprot_N/*.pdb` candidate structures |
|2| Build charge states *(see note above)* | - | AMBER ff14SB (protein) + Gaussian16/MK (ligand, `gaussian/*.log`) | `pb-electrostatics/apbs/*.pqr` (4 states) |
|3| Solve electrostatics | `pb-electrostatics/` | APBS (`pka_workflow.py`/`.ipynb`) | `.dx` potential maps, total energies |
|4| Compute pKa shift | `pb-electrostatics/` | isodesmic double-difference cycle | predicted ΔpKa |
|5| Visualize | `pb-electrostatics/` | PyMOL (using the generated script or headless via `pymol2`) | difference-map render (`.png`) |

See `pb-electrostatics/README.md` for the full thermodynamic-cycle derivation and `structure-predictions/README.md` for the ligand-variant details and Chai-1 pipeline reference.

## Environment setup

This repository has **two separate Python environments** - the root `pyproject.toml` covers only the `pb-electrostatics` stack (APBS, PyMOL). `structure-predictions` depends on `chai_lab`, which needs a GPU and is not part of the root environment; set it up separately as described below.

### 1. Install `uv`

[`uv`](https://docs.astral.sh/uv/) is a fast Python package/project manager used here to install and pin the `pb-electrostatics` dependencies.

**macOS / Linux:**
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

**Windows (PowerShell):**
```powershell
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

**Alternative (any platform, if you already have Python/pipx):**
```bash
pip install uv
# or
pipx install uv
```

Verify the install:
```bash
uv --version
```

### 2. Set up a common environment for both workflows

From the repository root (where `pyproject.toml` lives):

```bash
uv sync --all-groups
```

This reads `pyproject.toml` (Python >=3.11, pinned to 3.11 via `.python-version`), creates a `.venv/` in the repo root (if it has not already been created), resolves and installs the packages from all dependency groups defined in `pyproject.toml`, and writes a `uv.lock` file for reproducibility. Re-run `uv sync` any time `pyproject.toml` changes.


### 3. Set up the `pb-electrostatics` environment only

From the repository root (where `pyproject.toml` lives):

```bash
uv sync --group apbs --group dev
```

This reads `pyproject.toml` (Python >=3.11, pinned to 3.11 via `.python-version`), resolves and installs `apbs-binary`  and `pymol-open-source` packages from the dependency group `apbs` (plus default packages and `ipykernel` from the `dev` dependency group, installed by default), and writes a `uv.lock` file for reproducibility. Re-run `uv sync` any time `pyproject.toml` changes.

### 4. Run the electrostatics workflow

Either prefix commands with `uv run` (no manual activation needed):

```bash
cd pb-electrostatics
uv run python pka_workflow.py \
    --mn2-prot apbs/complex_MNP2.pqr \
    --mn2-deprot apbs/complex_MNP2_deprot.pqr \
    --mn3-prot apbs/complex_MNP3.pqr \
    --mn3-deprot apbs/complex_MNP3_deprot.pqr \
    --pdie 2.0 --outdir results

uv run jupyter notebook pka_workflow.ipynb
```

or activate the venv directly and run normally:

```bash
source .venv/bin/activate      # Windows: .venv\Scripts\activate
cd pb-electrostatics
python pka_workflow.py --help
```

See `pb-electrostatics/README.md` for the full CLI reference, the input PQR-state requirements, and example results.

### 5. Set up the `structure-predictions` environment (separate)

From the repository root (where `pyproject.toml` lives):

```bash
uv sync --group chai1
```

This reads `pyproject.toml` (Python >=3.11, pinned to 3.11 via `.python-version`), resolves and installs `chai_lab`, `numpy`, `biopython`, `gemmi`, and `biotite` packages from the dependency group `chai1` (plus default packages and `ipykernel` from the `dev` dependency group, installed by default), and writes a `uv.lock` file for reproducibility. Re-run `uv sync` any time `pyproject.toml` changes.

See `structure-predictions/README.md` for the full script reference,
per-ligand notes, and HPC/Colab usage.

## References

- Bashford, D. & Karplus, M. *pKa's of ionizable groups in proteins: atomic detail from a continuum electrostatic model.* Biochemistry **29**, 10219-10225 (1990).
- Baker, N. A., Sept, D., Joseph, S., Holst, M. J. & McCammon, J. A. *Electrostatics of nanosystems: application to microtubules and the ribosome.* Proc. Natl. Acad. Sci. USA **98**, 10037-10041 (2001).
