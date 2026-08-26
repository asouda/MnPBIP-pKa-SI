#!/usr/bin/env python3
"""
pka_workflow.py
================

End-to-end workflow to estimate a redox-linked pKa shift from four
protein-ligand PQR structures (protonated / deprotonated, at two metal
oxidation states) using APBS continuum electrostatics, combined in an
isodesmic-style double-difference thermodynamic cycle. Also produces a
PyMOL visualization of the resulting electrostatic potential difference
map.

Pipeline
--------
1. Run APBS (linearized Poisson-Boltzmann, focusing multigrid) on each of
   the four supplied PQR structures, producing one electrostatic potential
   grid (.dx) and one total electrostatic energy per state.
2. Combine the four energies in a thermodynamic square (isodesmic cycle)
   to obtain the electrostatic coupling free energy between the redox
   site and the protonation site, and convert it to a pKa shift.
3. Compute the electrostatic potential DIFFERENCE grid between the two
   protonated states (reduced minus oxidized) with APBS's own `dxmath`
   tool.
4. Write a ready-to-run PyMOL (.pml) script that visualizes that
   difference map as a pair of +/- isosurfaces on the protein/ligand
   structure.
5. Optionally, drive PyMOL itself (headless, via the `pymol2` API) to
   execute that script and save a rendered PNG snapshot.

Requires
--------
    pip install apbs-binary pymol-open-source

Usage
-----
    python pka_workflow.py \\
        --mn2-prot complex_MNP2.pqr \\
        --mn2-deprot complex_MNP2_deprot.pqr \\
        --mn3-prot complex_MNP3.pqr \\
        --mn3-deprot complex_MNP3_deprot.pqr \\
        --pdie 2.0 \\
        --outdir results \\
        --run-pymol

Run `python pka_workflow.py --help` for the full list of options.
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import math
import os
import re
import sys
from pathlib import Path
from turtle import title
from typing import Optional, Sequence

try:
    import apbs_binary
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        "This script requires the 'apbs-binary' package.\n"
        "Install it with:  pip install apbs-binary"
    ) from exc

import numpy as np


# ----------------------------------------------------------------------
# Physical constants
# ----------------------------------------------------------------------
GAS_CONSTANT_KJ = 8.314462618e-3  # kJ / (mol K)
KJ_TO_KCAL = 0.23900573613766727  # kJ / (mol K) * (1 kcal / 4.184 kJ)


# ----------------------------------------------------------------------
# Data structures
# ----------------------------------------------------------------------
@dataclasses.dataclass
class StateResult:
    """Result of one APBS electrostatics calculation for a single charge state."""

    label: str
    pqr_file: Path
    in_file: Path
    log_file: Path
    dx_file: Path
    total_energy_kjmol: float


@dataclasses.dataclass
class PKaShiftResult:
    """Result of the isodesmic double-difference pKa-shift calculation."""

    dG_deprot_reduced: float
    dG_deprot_oxidized: float
    ddG_coupling: float
    dpka: float
    rt_ln10: float


# ----------------------------------------------------------------------
# PQR helpers
# ----------------------------------------------------------------------
def read_pqr_extent(pqr_file: Path):
    """Return (xmin, xmax, ymin, ymax, zmin, zmax) over all ATOM/HETATM records.

    Coordinates are read from the 3 fields preceding charge and radius,
    counted from the end of the line -- this is robust to PQR files with
    or without an explicit chain-ID column.
    """
    xs, ys, zs = [], [], []
    with open(pqr_file) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                parts = line.split()
                x, y, z = (float(v) for v in parts[-5:-2])
                xs.append(x)
                ys.append(y)
                zs.append(z)
    if not xs:
        raise ValueError(f"No ATOM/HETATM records found in {pqr_file}")
    return min(xs), max(xs), min(ys), max(ys), min(zs), max(zs)


def compute_box_dimensions(pqr_file: Path, coarse_padding: float, fine_padding: float):
    """Auto-size APBS coarse/fine grid boxes from the molecular extent plus padding."""
    xmin, xmax, ymin, ymax, zmin, zmax = read_pqr_extent(pqr_file)
    extent = (xmax - xmin, ymax - ymin, zmax - zmin)
    cglen = tuple(e + 2 * coarse_padding for e in extent)
    fglen = tuple(e + 2 * fine_padding for e in extent)
    return cglen, fglen


def compute_common_box_dimensions(pqr_files: Sequence[Path], coarse_padding: float, fine_padding: float):
    """Auto-size ONE shared coarse/fine grid box that comfortably contains all
    of the given PQR structures (union of their extents).

    Using a single shared box (rather than sizing each state's grid off its
    own, possibly slightly different, atom count/extent) guarantees that all
    four charge states are solved on identical grids -- which is what lets
    their electrostatic energies be combined in a clean double difference,
    and what lets their .dx maps be subtracted point-for-point.
    """
    xmins, xmaxs, ymins, ymaxs, zmins, zmaxs = [], [], [], [], [], []
    for f in pqr_files:
        xmin, xmax, ymin, ymax, zmin, zmax = read_pqr_extent(f)
        xmins.append(xmin); xmaxs.append(xmax)
        ymins.append(ymin); ymaxs.append(ymax)
        zmins.append(zmin); zmaxs.append(zmax)
    extent = (max(xmaxs) - min(xmins), max(ymaxs) - min(ymins), max(zmaxs) - min(zmins))
    cglen = tuple(e + 2 * coarse_padding for e in extent)
    fglen = tuple(e + 2 * fine_padding for e in extent)
    return cglen, fglen


def check_pqr_coordinates_consistent(pqr_a: Path, pqr_b: Path, atol: float = 1e-3) -> bool:
    """Sanity check that two PQR files share one protein-embedded geometry.

    Compares the coordinates of every (atom name, residue number) pair the
    two files have in common (one file may be missing atoms present in the
    other -- e.g. a deprotonated state missing the leaving proton). This is
    important because the double-difference cancellation this workflow
    relies on assumes all four states sit on the *same* fixed geometry.
    """

    def read_coords(path: Path):
        coords = {}
        with open(path) as fh:
            for line in fh:
                if line.startswith(("ATOM", "HETATM")):
                    parts = line.split()
                    key = (parts[2], parts[4])  # (atom name, residue number)
                    coords[key] = tuple(float(v) for v in parts[-5:-2])
        return coords

    coords_a = read_coords(pqr_a)
    coords_b = read_coords(pqr_b)
    for key in set(coords_a) & set(coords_b):
        ax, ay, az = coords_a[key]
        bx, by, bz = coords_b[key]
        if abs(ax - bx) > atol or abs(ay - by) > atol or abs(az - bz) > atol:
            return False
    return True


# ----------------------------------------------------------------------
# APBS input generation and execution
# ----------------------------------------------------------------------
APBS_INPUT_TEMPLATE = """\
read
    mol pqr {pqr_file}
end

elec name {elec_name}
    mg-auto
    dime {dimex} {dimey} {dimez}
    cglen {cglenx:.3f} {cgleny:.3f} {cglenz:.3f}
    fglen {fglenx:.3f} {fgleny:.3f} {fglenz:.3f}
    cgcent mol 1
    fgcent mol 1
    mol 1
    lpbe
    bcfl mdh
    pdie {pdie:.4f}
    sdie {sdie:.4f}
    srfm mol
    chgm spl2
    srad 1.4
    swin 0.3
    sdens 10.0
    temp {temp:.2f}
    ion charge 1 conc {ionic_strength:.4f} radius 2.0
    ion charge -1 conc {ionic_strength:.4f} radius 1.8
    calcenergy total
    calcforce no
    write pot dx {dx_stem}
end

quit
"""


def build_apbs_input(pqr_file: Path, elec_name: str, dx_stem: str, pdie: float,
                      sdie: float, ionic_strength: float, temp: float,
                      dime: Sequence[int], cglen: Sequence[float], fglen: Sequence[float]) -> str:
    return APBS_INPUT_TEMPLATE.format(
        pqr_file=str(pqr_file),
        elec_name=elec_name,
        dimex=dime[0], dimey=dime[1], dimez=dime[2],
        cglenx=cglen[0], cgleny=cglen[1], cglenz=cglen[2],
        fglenx=fglen[0], fgleny=fglen[1], fglenz=fglen[2],
        pdie=pdie, sdie=sdie, temp=temp, ionic_strength=ionic_strength,
        dx_stem=dx_stem,
    )


def parse_fine_grid_energy(log_file: Path) -> float:
    """Extract the converged (focused/fine-grid) total electrostatic energy.

    APBS's `mg-auto` focusing prints 'Total electrostatic energy' twice per
    calculation: once for the coarse grid and once (final) for the focused
    fine grid. The fine-grid value -- the last one printed -- is the
    converged result we want.
    """
    pattern = re.compile(r"Total electrostatic energy\s*=\s*([-\d.eE+]+)\s*kJ/mol")
    values = [float(m.group(1)) for m in pattern.finditer(log_file.read_text())]
    if not values:
        raise ValueError(f"No 'Total electrostatic energy' found in {log_file}")
    return values[-1]


def run_apbs_calculation(pqr_file: Path, label: str, pdie: float, outdir: Path, *,
                          sdie: float = 78.54, ionic_strength: float = 0.150,
                          temp: float = 298.15, dime: Sequence[int] = (161, 161, 161),
                          cglen: Sequence[float], fglen: Sequence[float]) -> StateResult:
    """Write an APBS input file for one charge state, run APBS, and parse the result.

    `cglen`/`fglen` must be supplied by the caller (see
    `compute_common_box_dimensions`) so that every state in a given workflow
    run is solved on an identical grid.
    """
    pqr_file = Path(pqr_file).resolve()

    dx_stem = f"pot_{label}"
    in_text = build_apbs_input(pqr_file, label, dx_stem, pdie, sdie, ionic_strength,
                                temp, dime, cglen, fglen)

    in_file = outdir / f"{label}.in"
    log_file = outdir / f"{label}.log"
    in_file.write_text(in_text)

    print(f"[apbs] running '{label}'  (pdie={pdie:g}) ...")
    # NOTE: cwd is set to outdir below, so the input-file argument must be
    # given relative to outdir (just the filename), not as outdir/filename.
    result = apbs_binary.run_apbs(in_file.name, cwd=str(outdir), capture_output=True, text=True)
    log_file.write_text((result.stdout or "") + (result.stderr or ""))
    if result.returncode != 0:
        raise RuntimeError(f"APBS failed for '{label}' (see {log_file})")

    energy = parse_fine_grid_energy(log_file)

    # Depending on how APBS was built (serial vs. rank-aware), the written
    # map is named either "pot_<label>.dx" or "pot_<label>-PE0.dx". Accept
    # either.
    candidates = [outdir / f"{dx_stem}.dx", outdir / f"{dx_stem}-PE0.dx"]
    dx_file = next((c for c in candidates if c.exists()), None)
    if dx_file is None:
        raise FileNotFoundError(
            f"Expected APBS output not found (tried: {', '.join(str(c) for c in candidates)})"
        )

    print(f"[apbs] '{label}': total electrostatic energy = {energy*KJ_TO_KCAL:.3f} kcal/mol")
    return StateResult(label, pqr_file, in_file, log_file, dx_file, energy)


# ----------------------------------------------------------------------
# pKa shift (isodesmic double-difference cycle)
# ----------------------------------------------------------------------
def compute_pka_shift(G_reduced_prot: float, G_reduced_deprot: float,
                       G_oxidized_prot: float, G_oxidized_deprot: float,
                       temp: float = 298.15) -> PKaShiftResult:
    """Combine the four state energies into a pKa shift.

    Because only the *shift* in pKa (not its absolute value) is sought, the
    intrinsic/model-compound deprotonation free energy and the proton's own
    reference chemical potential are identical in both oxidation states and
    cancel exactly in the double difference below -- no bulk-solvent
    reference-compound calculation is required.
    """
    dG_reduced = G_reduced_deprot - G_reduced_prot
    dG_oxidized = G_oxidized_deprot - G_oxidized_prot
    ddG = dG_oxidized - dG_reduced
    rt_ln10 = GAS_CONSTANT_KJ * temp * math.log(10)
    dpka = ddG / rt_ln10
    return PKaShiftResult(dG_reduced, dG_oxidized, ddG, dpka, rt_ln10)


# ----------------------------------------------------------------------
# Difference map (dxmath)
# ----------------------------------------------------------------------
def compute_difference_map(dx_a: Path, dx_b: Path, out_name: str, outdir: Path) -> Path:
    """Compute dx_a - dx_b with APBS's dxmath tool; returns the output .dx path."""
    mie_script = outdir / "diffmap.mie"
    mie_script.write_text(f"{dx_a.name}\n{dx_b.name} -\n{out_name} =\n")

    print(f"[dxmath] computing {out_name} = {dx_a.name} - {dx_b.name}")
    result = apbs_binary.run_dxmath(mie_script.name, cwd=str(outdir), capture_output=True, text=True)
    out_dx = outdir / out_name
    if result.returncode != 0 or not out_dx.exists():
        raise RuntimeError(f"dxmath failed:\n{result.stdout}\n{result.stderr}")
    return out_dx


def read_dx_values(dx_file: Path) -> np.ndarray:
    """Load the scalar grid values out of an OpenDX file as a flat numpy array."""
    lines = Path(dx_file).read_text().splitlines()
    start = next(i for i, l in enumerate(lines) if "data follows" in l) + 1
    values = []
    for line in lines[start:]:
        if line.startswith("attribute"):
            break
        values.extend(float(v) for v in line.split())
    return np.array(values)


def suggest_isovalue(dx_file: Path, sigma_multiple: float = 3.0) -> float:
    """Suggest a +/- isosurface level from a map's own value distribution.

    Raw (unsmoothed) APBS potential/difference maps carry large near-atom
    self-energy spikes (see the accompanying notebook for a full discussion);
    a fixed isovalue that looks reasonable at one pdie/grid combination can
    engulf the entire structure at another. Scaling the isovalue to a
    multiple of the map's own standard deviation keeps the isosurface
    informative (showing structured signal, not background noise or a
    single solid blob) across different runs.
    """
    values = read_dx_values(dx_file)
    return float(sigma_multiple * values.std())


# ----------------------------------------------------------------------
# PyMOL visualization
# ----------------------------------------------------------------------
PML_TEMPLATE = """\
# {title}
bg_color white
load {structure_pqr}, complex

hide everything
show surface, complex
set transparency, 0.15, complex
show cartoon, complex
show sticks, resn MNP

load {diff_dx}, diffmap
ramp_new diff_ramp, diffmap, [{minvalue}, {midvalue}, {maxvalue}], [red, white, blue]
set surface_color, diff_ramp, complex

select hn6_atom, resn MNP and name HN6
set sphere_scale, 0.25, hn6_atom
show spheres, hn6_atom
color yellow, hn6_atom
deselect

select mn_atom, resn MNP and name MN
show spheres, mn_atom
deselect

set_view (\\
     0.313743979,   -0.820566177,    0.477735579,\\
    -0.942059696,   -0.331898630,    0.048604969,\\
     0.118676759,   -0.465306908,   -0.877151370,\\
     0.000000000,    0.000000000, -170.399078369,\\
    36.760910034,   29.815954208,   27.141380310,\\
   134.343933105,  206.454223633,  -20.000000000 )

zoom resn MNP, 15
"""


def write_pml_script(structure_pqr: Path, diff_dx: Path, out_pml: Path, *,
                      minvalue: float = -2.0,
                      maxvalue: float = 2.0,
                      title: str = "Electrostatic potential difference map") -> Path:
    """Write a standalone PyMOL script that visualizes the difference map.

    Absolute paths are used throughout so the script can be run with
    `pymol out_pml` from any working directory.
    """
    text = PML_TEMPLATE.format(
        title=title,
        structure_pqr=str(Path(structure_pqr).resolve()),
        diff_dx=str(Path(diff_dx).resolve()),
        minvalue=minvalue,
        maxvalue=maxvalue,
        midvalue=minvalue+(maxvalue-minvalue)/2,
    )
    out_pml.write_text(text)
    print(f"[pymol] wrote visualization script: {out_pml}")
    return out_pml


def run_pymol_visualization(pml_file: Path, png_out: Optional[Path] = None,
                             width: int = 1200, height: int = 900, dpi: int = 150) -> None:
    """Execute a .pml script headlessly via the pymol2 API and optionally save a PNG."""
    try:
        import pymol2
    except ImportError as exc:  # pragma: no cover
        raise SystemExit(
            "PyMOL visualization requires the 'pymol-open-source' package.\n"
            "Install it with:  pip install pymol-open-source"
        ) from exc

    print("[pymol] launching headless PyMOL session ...")
    session = pymol2.PyMOL()
    session.start()
    try:
        session.cmd.do(f"@{pml_file}")
        if png_out is not None:
            session.cmd.ray(width, height)
            session.cmd.png(str(png_out), dpi=dpi)
            print(f"[pymol] saved image: {png_out}")
    finally:
        session.stop()


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------
def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="APBS + PyMOL workflow: redox-linked pKa shift from four PQR charge states.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--mn2-prot", required=True, type=Path,
                    help="PQR file: reduced state, protonated (e.g. Mn(II), N-H present)")
    p.add_argument("--mn2-deprot", required=True, type=Path,
                    help="PQR file: reduced state, deprotonated")
    p.add_argument("--mn3-prot", required=True, type=Path,
                    help="PQR file: oxidized state, protonated (e.g. Mn(III), N-H present)")
    p.add_argument("--mn3-deprot", required=True, type=Path,
                    help="PQR file: oxidized state, deprotonated")

    p.add_argument("--pdie", type=float, default=2.0,
                    help="Solute (interior) dielectric constant")
    p.add_argument("--sdie", type=float, default=78.54, help="Solvent dielectric constant")
    p.add_argument("--ionic-strength", type=float, default=0.150,
                    help="1:1 (monovalent) ionic strength, in mol/L")
    p.add_argument("--temp", type=float, default=298.15, help="Temperature, in K")

    p.add_argument("--dime", type=int, nargs=3, default=(161, 161, 161), metavar=("NX", "NY", "NZ"),
                    help="APBS grid point counts (same for coarse and fine grids)")
    p.add_argument("--cglen", type=float, nargs=3, default=None, metavar=("X", "Y", "Z"),
                    help="Coarse grid box lengths in Angstrom [default: auto from molecular extent + padding]")
    p.add_argument("--fglen", type=float, nargs=3, default=None, metavar=("X", "Y", "Z"),
                    help="Fine grid box lengths in Angstrom [default: auto from molecular extent + padding]")

    p.add_argument("--minvalue", type=float, default=-2.0,
                    help="minimum value of electrostatic potential (kT/e) for the PyMOL difference-map visualization "
                         "[default: -2.0]")
    p.add_argument("--maxvalue", type=float, default=2.0,
                    help="maximum value of electrostatic potential (kT/e) for the PyMOL difference-map visualization "
                         "[default: +2.0]")
    p.add_argument("--outdir", type=Path, default=Path("pka_workflow_output"),
                    help="Directory for all generated files")
    p.add_argument("--run-pymol", action="store_true",
                    help="Also launch headless PyMOL, render, and save a PNG of the difference map")

    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> dict:

    args = parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("Redox-linked pKa shift workflow  (APBS + PyMOL)")
    print("=" * 70)
    print(f"pdie = {args.pdie:g}   outdir = {args.outdir}")

    # ---- 1. sanity-check that all four states share one geometry --------
    checks = [
        (args.mn2_prot, args.mn2_deprot, "Mn2+ prot vs deprot"),
        (args.mn3_prot, args.mn3_deprot, "Mn3+ prot vs deprot"),
        (args.mn2_prot, args.mn3_prot, "Mn2+ vs Mn3+ (protonated)"),
    ]
    for a, b, tag in checks:
        ok = check_pqr_coordinates_consistent(a, b)
        print(f"[check] shared-atom coordinate consistency ({tag}): {'OK' if ok else 'MISMATCH'}")
        if not ok:
            print(f"  WARNING: {a} and {b} disagree on the coordinates of one or more shared atoms.")
            print("  This workflow's double-difference cancellation assumes one shared geometry")
            print("  across all four states -- treat the resulting pKa shift with caution.")

    # ---- 2. run APBS on each of the four states --------------------------
    # One shared grid box for all four states (see compute_common_box_dimensions
    # docstring for why this matters), computed from the union of their extents
    # unless the user pinned explicit --cglen/--fglen values.
    all_pqrs = [args.mn2_prot, args.mn2_deprot, args.mn3_prot, args.mn3_deprot]
    if args.cglen and args.fglen:
        cglen, fglen = tuple(args.cglen), tuple(args.fglen)
    else:
        auto_cglen, auto_fglen = compute_common_box_dimensions(all_pqrs, coarse_padding=40.0, fine_padding=15.0)
        cglen = tuple(args.cglen) if args.cglen else auto_cglen
        fglen = tuple(args.fglen) if args.fglen else auto_fglen
    print(f"[grid] cglen = {cglen[0]:.2f} x {cglen[1]:.2f} x {cglen[2]:.2f} A "
          f"(coarse), fglen = {fglen[0]:.2f} x {fglen[1]:.2f} x {fglen[2]:.2f} A (fine)")

    common_grid_kwargs = dict(
        sdie=args.sdie, ionic_strength=args.ionic_strength, temp=args.temp,
        dime=tuple(args.dime), cglen=cglen, fglen=fglen,
    )
    states = {}
    for label, pqr in [("mn2_prot", args.mn2_prot), ("mn2_deprot", args.mn2_deprot),
                        ("mn3_prot", args.mn3_prot), ("mn3_deprot", args.mn3_deprot)]:
        states[label] = run_apbs_calculation(pqr, label, args.pdie, args.outdir, **common_grid_kwargs)

    # ---- 3. isodesmic double-difference pKa shift -------------------------
    pka = compute_pka_shift(
        states["mn2_prot"].total_energy_kjmol, states["mn2_deprot"].total_energy_kjmol,
        states["mn3_prot"].total_energy_kjmol, states["mn3_deprot"].total_energy_kjmol,
        temp=args.temp,
    )
    print("-" * 70)
    print(f"dG_deprot(reduced)   = {pka.dG_deprot_reduced*KJ_TO_KCAL:9.3f} kcal/mol")
    print(f"dG_deprot(oxidized)  = {pka.dG_deprot_oxidized*KJ_TO_KCAL:9.3f} kcal/mol")
    print(f"ddG_coupling         = {pka.ddG_coupling*KJ_TO_KCAL:9.3f} kcal/mol")
    print(f"RT ln10              = {pka.rt_ln10*KJ_TO_KCAL:7.4f} kcal/mol")
    print(f"Predicted pKa shift (oxidized - reduced) = {pka.dpka:+.2f}")
    print("-" * 70)

    # ---- 4. electrostatic potential difference map -----------------------
    diff_dx = compute_difference_map(
        states["mn2_prot"].dx_file, states["mn3_prot"].dx_file,
        "pot_diff_mn2_minus_mn3.dx", args.outdir,
    )

    # ---- 5. PyMOL visualization script ------------------------------------
    minvalue = args.minvalue
    maxvalue = args.maxvalue
    print(f"[pymol] selected range of electrostatic potential (kT/e) for the PyMOL difference-map visualization: {minvalue:.3f} to {maxvalue:.3f}")
    pml_file = write_pml_script(
        states["mn2_prot"].pqr_file, diff_dx, args.outdir / "diff_map_viz.pml",
        minvalue=minvalue,
        maxvalue=maxvalue,
    )

    # ---- 6. summary --------------------------------------------------------
    summary = {
        "pdie": args.pdie,
        "energies_kcal_per_mol": {k: v.total_energy_kjmol*KJ_TO_KCAL for k, v in states.items()},
        "dG_deprot_reduced_kcal_per_mol": pka.dG_deprot_reduced*KJ_TO_KCAL,
        "dG_deprot_oxidized_kcal_per_mol": pka.dG_deprot_oxidized*KJ_TO_KCAL,
        "ddG_coupling_kcal_per_mol": pka.ddG_coupling*KJ_TO_KCAL,
        "RT_ln10_kcal_per_mol": pka.rt_ln10*KJ_TO_KCAL,
        "delta_pKa": pka.dpka,
        "difference_map": str(diff_dx),
        "pml_script": str(pml_file),
        "minvalue": minvalue,
        "maxvalue": maxvalue,
    }
    summary_file = args.outdir / "summary.json"
    summary_file.write_text(json.dumps(summary, indent=2))
    print(f"[summary] wrote {summary_file}")

    # ---- 7. optional: render with PyMOL -------------------------------------
    if args.run_pymol:
        png_out = args.outdir / "diff_map_view.png"
        run_pymol_visualization(pml_file, png_out=png_out)

    print("Done.")
    return summary


if __name__ == "__main__":
    main()
