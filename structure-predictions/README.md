# R_smuprot — Chai-1 Structure Prediction for Mn-Porphyrin Complexes

This directory predicts the 3-D structure of a 115-residue target protein bound to
**six** different Mn-porphyrin ligand variants (`smuprot_0` … `smuprot_4`, plus a
protonated variant of `smuprot_4`) using **Chai-1**, a deep-learning structure
predictor from Chai Discovery.

---

## Directory layout

```
R_smuprot/
├── smuprot_0.json                 Input: protein + Mn(III)-tetraphenylporphyrin (vinyl)
├── smuprot_1.json                 Input: protein + fluorinated Mn-porphyrin (Mn+2)
├── smuprot_2.json                 Input: protein + fluorinated Mn-porphyrin (Mn-2, same molecule as smuprot_1)
├── smuprot_3.json                 Input: protein + bis-CF3 fluorinated Mn-porphyrin ([Mn] neutral notation)
├── smuprot_4.json                 Input: protein + bis-CF3 fluorinated Mn-porphyrin (Mn-2, same molecule as smuprot_3)
├── smuprot_4_protonated.json      Input: same ligand as smuprot_4, with the peripheral imidazole N protonated ([NH+])
│
├── json_to_fasta.py               JSON → Chai-1 FASTA converter  ← start here
├── run_batch.py                   Iterates FASTA files, calls run_prediction + process_best
├── run_prediction.py              Calls Chai-1 run_inference for one FASTA
├── process_best.py                Picks highest-confidence CIF, converts to PDB, writes summary
├── imports.py                     Shared imports and helpers used by all scripts above
│
├── run.sh                         Shell entry point: runs the full pipeline
├── run_protein_mn_ligand.ipynb    Colab notebook with per-JSON cells and explanatory text
├── requirements.txt               pip dependencies
│
├── smuprot_0/                     Created at runtime by json_to_fasta.py + run.sh
│   ├── smuprot_0.fasta            Generated FASTA (protein + ligand SMILES)
│   ├── smuprot_0.cif              Best predicted structure (mmCIF)
│   └── smuprot_0.pdb              Best predicted structure (PDB)
├── smuprot_1/ … smuprot_4_protonated/   (same structure as smuprot_0/, one folder per JSON above)
│
└── protein_mn_ligand_0/           Leftover manual test run (ligand-0 FASTA only, pipeline
                                    not re-run here — kept for reference, not part of the
                                    smuprot_N pipeline)
```

> **Note:** `process_best.py` writes a `summary_*.txt` file with score statistics
> into the *intermediate* `smuprot_N/chai1/` folder, but `run_batch.py` deletes
> that folder as a cleanup step once the final CIF/PDB have been copied out — so
> the summary file does **not** persist in `smuprot_N/`. Only the `.fasta`,
> `.cif`, and `.pdb` files remain after a full `run.sh` pass.

---

## How to run

### On NERSC / HPC (command line)

```bash
cd R_smuprot/

# Run all smuprot_*.json files
./run.sh

# Run one JSON only
./run.sh smuprot_0

# Use Chai-2 engine instead of Chai-1
./run.sh smuprot_1 chai2
```

### On Google Colab

Open `run_protein_mn_ligand.ipynb`.  The notebook has a cell for each JSON and
explanatory text throughout.  Run all cells top to bottom.

> The notebook's "Results Summary" cell currently only reports `smuprot_0`,
> `smuprot_1`, and `smuprot_2` — it has not been extended to also collect
> results for `smuprot_3`, `smuprot_4`, and `smuprot_4_protonated`, even though
> the notebook does run predictions for all six JSONs.

---

## Script reference

### `json_to_fasta.py`

**Purpose:** Reads one or more `smuprot_N.json` files and writes a Chai-1-compatible
FASTA file into a matching subfolder (`smuprot_N/smuprot_N.fasta`).

The ligand SMILES is copied **verbatim** from the JSON — the Mn metal centre is
**not modified or stripped**.

**The UFF warning explained:** Chai-1 uses RDKit internally to build an initial
3-D conformer of the ligand.  RDKit's UFF force-field typer does not have
parameters for transition metals (Mn, Fe, Cu, …), so it prints:

```
UFFTYPER: Unrecognised atom type: Mn5+2 (12)
```

This is a **non-fatal warning**.  RDKit builds the conformer via its
distance-geometry algorithm (ETKDG) first; this method does not require
force-field atom types and succeeds for metal-containing molecules.  The UFF
optimisation step is skipped.  The Mn and all its bonds are fully present in
the conformer that Chai-1 uses.  The warning is suppressed by default
(`RDLogger.DisableLog`) to avoid alarming users during batch runs.

```bash
# Convert all six JSONs at once
python json_to_fasta.py --json smuprot_0.json smuprot_1.json smuprot_2.json \
                                smuprot_3.json smuprot_4.json smuprot_4_protonated.json

# Show metal info and UFF note (useful for debugging)
python json_to_fasta.py --json smuprot_0.json --verbose

# Write to a custom folder
python json_to_fasta.py --json smuprot_0.json --outdir my_output/
```

---

### `run_batch.py`

**Purpose:** Iterates all `*.fasta` files in a given protein folder, calls
`run_prediction.run_ultrahigh()` on each one, then calls
`process_best.process_output_dir()` to pick the best structure and write the
final CIF/PDB.  Cleans up the intermediate per-engine folder afterwards.

```bash
python run_batch.py --protein smuprot_0               # single folder
python run_batch.py --protein all                     # all folders with *.fasta
python run_batch.py --protein smuprot_1 --engine chai2
```

| Argument | Default | Description |
|---|---|---|
| `--protein` | *(required)* | Protein folder name or `all` |
| `--workspace` | current directory | Root containing protein folders |
| `--engine` | `chai1` | `chai1` or `chai2` |

---

### `run_prediction.py`

**Purpose:** Calls `chai_lab.chai1.run_inference` (or `chai2`) for a single
FASTA file with the ULTRA-HIGH quality configuration:

| Parameter | Value | Meaning |
|---|---|---|
| `num_trunk_recycles` | 3 | Evoformer recycle passes |
| `num_diffn_timesteps` | 500 | Diffusion steps (more = better geometry) |
| `num_diffn_samples` | 5 | Independent structure samples per run |
| `use_esm_embeddings` | True | Use ESM-2 sequence embeddings |

Outputs are written to `<output_root>/<complex_name>/pred.model_idx_N.cif` and
matching `scores.model_idx_N.npz` files.  A `metadata.json` is also written for
downstream processing.

```bash
python run_prediction.py --fasta smuprot_0/smuprot_0.fasta
python run_prediction.py --fasta smuprot_1/smuprot_1.fasta --engine chai2
```

---

### `process_best.py`

**Purpose:** Reads the `scores.model_idx_N.npz` files from a Chai-1 output
directory, selects the CIF with the highest confidence score, converts it to
PDB format, relabels any residues named `LIG` to `MNP` (the Mn-porphyrin), and
writes a plain-text summary.

**CIF → PDB conversion** is attempted first with **biotite**; if that fails,
**gemmi** is used as a fallback.

**Ligand relabelling:** Chai-1 names all small-molecule ligands `LIG` in the
output PDB.  `relabel_ligands_in_pdb()` changes chain B → `MNP` by default, so
the Mn-porphyrin has an informative residue name for downstream tools (PyMOL,
OpenMM, etc.).  Edit the `chain_to_resname` dictionary if Chai-1 places the
ligand on a different chain.

```bash
python process_best.py --output-dir smuprot_0/chai1/protein_mn_ligand \
                       --final-dir smuprot_0/
```

---

### `imports.py`

**Purpose:** Centralises all standard-library and third-party imports plus two
shared helper functions used by `run_batch.py`, `run_prediction.py`, and
`process_best.py`.

| Function | Description |
|---|---|
| `get_chai_run_inference(engine)` | Returns the `run_inference` callable from `chai_lab.chai1` or `chai_lab.chai2` |
| `parse_fasta_stem_to_dirname(stem)` | Converts a FASTA filename stem to an output directory name.  Handles both the legacy `Protein_with_X_and_PMP` naming (from the original R/ pipeline) and simple stems like `smuprot_0` (used here). |

---

## Dependencies

Install with:
```bash
pip install chai_lab gemmi biopython biotite numpy tqdm
```

See `requirements.txt` for the full list.

Chai-1 also requires a GPU (CUDA) for practical runtimes.  On CPU a single
prediction with 500 steps and 5 samples takes several hours.

---

## Notes on the six ligands

| File | Ligand description | Key features |
|---|---|---|
| `smuprot_0.json` | Mn(III)-tetraphenylporphyrin with vinyl groups | Four *meso* phenyl groups; β-pyrrole vinyls; `[n+]` charge notation |
| `smuprot_1.json` | Mn-porphyrin with pentafluorophenyl, –OH, –OMe | Electron-withdrawing C₆F₅; H-bond donors/acceptors; `[Mn+2]` notation |
| `smuprot_2.json` | Same molecule as smuprot_1, different charge convention | Identical compound; `[Mn-2]` / `[N+]` notation; useful consistency check |
| `smuprot_3.json` | Same scaffold as smuprot_1/2, plus two extra –CF₃ groups | Bis-trifluoromethyl substitution; neutral `[Mn]` with two `[N+]` pyrroles |
| `smuprot_4.json` | Same molecule as smuprot_3, different charge convention | Identical bis-CF₃ compound; `[Mn-2]` / three `[N+]` notation; consistency check vs. smuprot_3 |
| `smuprot_4_protonated.json` | Same molecule as smuprot_4, peripheral N protonated | The imidazole-type N on the *meso*-aryl substituent is written `[NH+]` instead of neutral `N`, modelling a protonated state of that group |

`smuprot_1`/`smuprot_2` and `smuprot_3`/`smuprot_4` are each intended as internal
consistency checks (same molecule, different SMILES charge bookkeeping) — see the
"Consistency check" discussion in `run_protein_mn_ligand.ipynb`.
