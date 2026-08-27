# Redox-Linked pKa Shift Workflow (APBS + PyMOL)

Estimates how the pKa of a titratable group shifts when a nearby metal
center changes oxidation state, using continuum electrostatics (APBS) and
an isodesmic-style double-difference thermodynamic cycle. Also generates
the electrostatic potential difference map between the two protonated
states and a PyMOL visualization of it.

Built for a Mn(II)/Mn(III) porphyrin complex and its axial imidazolium
N-H, but the script takes any four PQR files as input and is not specific
to this system.

## Files

| File | Description |
|---|---|
| `pka_workflow.py` | Command-line script. Runs the full pipeline end to end. |
| `pka_workflow.ipynb` | Tutorial notebook version, same underlying functions, with markdown narration between steps. Already executed once, so you can read it without rerunning anything. |
| `requirements.txt` | Python dependencies. |

## Package Requirements

- **`apbs-binary`** - a self-contained pip-installable build of APBS
  (Adaptive Poisson-Boltzmann Solver) and its companion command-line
  tools (`dxmath`, `multivalue`, `smooth`, ...). No separate system
  install of APBS is needed. Each tool is exposed as a Python function,
  e.g. `apbs_binary.run_apbs(...)`.
- **`pymol-open-source`** - the open-source PyMOL build, importable as a
  normal package. Its `pymol2` submodule drives a full PyMOL session
  headlessly (no GUI/display needed), used for the optional rendering
  step.
- `numpy`, `pandas` - array/table handling (`pandas` is only used for
  the summary table in the notebook).
- `ipykernel` - only needed if you want to *execute* the notebook
  yourself rather than just read it.

Tested with `apbs-binary==3.4.1.1` and `pymol-open-source==3.2.0a0`.

## Quick start

```bash
python pka_workflow.py \
    --mn2-prot complex_MNP2.pqr \
    --mn2-deprot complex_MNP2_deprot.pqr \
    --mn3-prot complex_MNP3.pqr \
    --mn3-deprot complex_MNP3_deprot.pqr \
    --pdie 2.0 \
    --outdir results \
    --run-pymol
```

Total runtime is roughly 2-4 minutes (four APBS solves plus, if
`--run-pymol` is given, one headless PyMOL render).

## The four input states

The script expects four PQR files representing every combination of
{reduced / oxidized} x {protonated / deprotonated}:

| CLI flag | Meaning | Example | Ligand net charge |
|---|---|---|---|
| `--mn2-prot`   | reduced state, protonated  | Mn(II), N-H present  | +1 |
| `--mn2-deprot` | reduced state, deprotonated | Mn(II), N-H removed  |  0 |
| `--mn3-prot`   | oxidized state, protonated | Mn(III), N-H present | +2 |
| `--mn3-deprot` | oxidized state, deprotonated | Mn(III), N-H removed | +1 |

**All four files must share one fixed atomic geometry** (identical
coordinates and radii for every atom they have in common - one file may
simply be missing the leaving proton). Only the ligand's partial charges,
and the presence/absence of that one proton, should differ between them.
This is what lets the double-difference cancellation in the pKa
calculation, and the point-by-point subtraction in the difference map,
work cleanly. The script checks this automatically at startup
(`check_pqr_coordinates_consistent`) and prints a warning - not a hard
failure - if it finds a mismatch.

In the reference system this was built for, the protonated-state ligand
charges came from a DFT/MK (Merz-Kollman) electrostatic potential fit at
each oxidation state (B3LYP/Def2-TZVP, high-spin), and the
deprotonated-state charges from the same kind of fit with the proton
removed; protein atoms carried one consistent force field (AMBER ff14SB)
throughout.

## What it does

1. **Runs APBS** on each of the four states (linearized Poisson-Boltzmann
   equation, `mg-auto` focusing multigrid, molecular/Connolly dielectric
   surface, cubic-spline charge discretization, 150 mM 1:1 ionic
   strength by default), producing one electrostatic potential grid
   (`.dx`) and one total electrostatic energy per state.

   All four states are solved on **one shared grid box**, auto-sized
   from the union of all four structures' atomic extents plus padding
   (or pinned explicitly with `--cglen`/`--fglen`). This matters: if
   each state's box were sized off its own atoms independently, removing
   the leaving proton could shift the box just enough to misalign the
   four grids, which would break both the difference-map subtraction and
   the energy double-difference below.

2. **Computes the pKa shift** via the thermodynamic square:

   ```
      Mn(II), protonated  --- oxidation --->  Mn(III), protonated
             |                                        |
      dG_deprot(Mn2+)                          dG_deprot(Mn3+)
             |                                        |
             v                                        v
      Mn(II), deprotonated -- oxidation --->  Mn(III), deprotonated
   ```

   ```
   dG_deprot(Mn^n+) = G(Mn^n+, deprot) - G(Mn^n+, prot)
   ddG_coupling      = dG_deprot(Mn3+) - dG_deprot(Mn2+)
   dpKa              = ddG_coupling / (RT ln10)
   ```

   Because this is a double difference, the intrinsic (bulk-solvent /
   model-compound) deprotonation free energy of the titratable group and
   the proton's own reference chemical potential are identical in both
   oxidation-state legs and cancel out exactly - no separate small-molecule
   reference calculation is needed, as long as it's the *shift* you're
   after and not the absolute pKa. This follows the classic
   Bashford & Karplus (*Biochemistry* **29**, 10219-10225, 1990)
   electrostatic pKa/redox-coupling formalism.

3. **Computes the difference map** (reduced minus oxidized, protonated
   states) with APBS's `dxmath` tool - a direct point-by-point
   subtraction of the two `.dx` grids.

4. **Writes a PyMOL script** (`.pml`) that loads the structure and the
   difference map and renders it as +/- isosurfaces. The isosurface
   level defaults to **3x the difference map's own standard deviation**
   (override with `--isovalue`), computed at runtime rather than
   hardcoded - a fixed level like +/-1 kT/e can look fine at one
   `pdie`/grid combination and swallow the entire structure in a solid
   blob at another, because raw (unsmoothed) APBS maps carry large,
   sharp values right at each atom's own position (a self-energy
   artifact of representing a point charge on a finite grid, not a
   real feature of the field).

5. **(Optional, `--run-pymol`)** drives PyMOL headlessly via `pymol2` to
   execute that script and save a rendered PNG.

## Outputs

Everything lands in `--outdir` (default `pka_workflow_output/`):

```
<label>.in, <label>.log          APBS input/log for each of the 4 states
pot_<label>.dx                   electrostatic potential map for each state
pot_diff_mn2_minus_mn3.dx        difference map (reduced minus oxidized)
diff_map_viz.pml                 standalone PyMOL visualization script
diff_map_view.png                rendered image (only with --run-pymol)
summary.json                     all energies, ddG, and the predicted pKa shift
```

`diff_map_viz.pml` is fully standalone (absolute paths baked in) - you
can run it directly against any PyMOL install later, independent of this
script:

```bash
pymol diff_map_viz.pml
```

## CLI options

| Flag | Default | Meaning |
|---|---|---|
| `--mn2-prot`, `--mn2-deprot`, `--mn3-prot`, `--mn3-deprot` | *(required)* | the four PQR files |
| `--pdie` | `2.0` | interior (solute) dielectric constant |
| `--sdie` | `78.54` | solvent dielectric constant |
| `--ionic-strength` | `0.150` | 1:1 ionic strength, mol/L |
| `--temp` | `298.15` | temperature, K |
| `--dime` | `161 161 161` | grid points per axis (coarse and fine) |
| `--cglen`, `--fglen` | auto | coarse/fine box lengths, Angstrom (auto-computed from the four structures' union extent if omitted) |
| `--isovalue` | auto | +/- isosurface level, kT/e (auto = 3 x std. dev. of the difference map if omitted) |
| `--outdir` | `pka_workflow_output` | output directory |
| `--run-pymol` | off | also render and save a PNG |

Run `python pka_workflow.py --help` for the same information from the
script itself.

## Example results

From a real run on the reference Mn porphyrin system:

| pdie | ddG_coupling (kJ/mol) | predicted dpKa |
|---|---|---|
| 2.0 | -50.17 | -8.79 |
| 4.0 | -26.07 | -4.57 |

(auto-sized grid box: coarse 119.5 x 123.1 x 119.5 Angstrom, fine
69.5 x 73.1 x 69.5 Angstrom, 161^3 grid points, ~40-55 s per APBS solve.)

A negative shift means oxidation makes the group *more* acidic (lower
pKa) - consistent with the extra positive charge on the oxidized metal
destabilizing the protonated (more positively charged) form of the
titratable group more than the neutral, deprotonated form.

The coupling energy (and therefore the predicted shift) scales roughly
as `1/pdie` over this range, which is expected for an electrostatic
interaction dominated by direct charge-charge terms in a low-polarity,
largely rigid interior, and is a useful internal-consistency check when
you re-run at a different `pdie`.

## The notebook

`pka_workflow.ipynb` walks through the exact same steps as the script,
one at a time, with markdown explanations in between - importing and
calling the actual functions from `pka_workflow.py` rather than
duplicating logic, so there's one source of truth. It has already been
executed once (`jupyter nbconvert --execute`), so all cell outputs -
including the printed energies, the summary table, and the rendered
PNG - are baked into the file; you can read it top to bottom without
running anything.

To re-execute it yourself (e.g. with different input files):

```bash
pip install -r requirements.txt
python -m ipykernel install --user --name python3 --display-name "Python 3"
jupyter notebook pka_workflow.ipynb
```

`pka_workflow.py` must be in the same directory as the notebook, since
it's imported directly (`from pka_workflow import ...`).

## Caveats and limitations

- **Frozen geometry.** All four states share one fixed atomic geometry;
  no relaxation of the ligand or protein upon protonation or oxidation
  is modeled.
- **Dielectric constant.** The predicted shift's magnitude is sensitive
  to `--pdie`, roughly inversely proportional to it over a reasonable
  range (see table above). The right value for a given system isn't
  known a priori and ideally should be validated against experiment.
- **Linearized Poisson-Boltzmann.** May be less accurate very close to a
  densely-charged, buried region (small `pdie`, large local potentials).
- **Single static structure.** No conformational sampling is performed.
  Treat the predicted shift as a semi-quantitative estimate of sign and
  rough magnitude, not a rigorous free energy.

## Troubleshooting

- **`.dx` file naming.** Depending on how the underlying APBS binary was
  built, the written potential map is named either `pot_<label>.dx` or
  `pot_<label>-PE0.dx` (the `-PE0` suffix indicates a parallel-aware
  build labeling rank 0's output). `run_apbs_calculation` checks for
  both automatically - if you're inspecting `--outdir` by hand, be aware
  either name is possible.
- **APBS input file paths.** If you're extending the script and calling
  `apbs_binary.run_apbs` yourself, remember that the `cwd` argument
  changes the subprocess's working directory - pass the input filename
  alone (not `outdir/filename`) or you'll get a "cannot open input file"
  error.
- **`NoSuchKernel: python3`** when executing the notebook: run
  `python -m ipykernel install --user --name python3` first (see above).
