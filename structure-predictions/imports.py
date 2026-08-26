"""Centralized imports and shared helpers for the protein_mn_ligand Chai-1 pipeline.

Mirrors R/imports.py with naming logic adapted for single-ligand smuprot.json entries
(no _with_..._and_PMP pattern — the FASTA stem IS the directory name).
"""

# Environment tweaks
import os
os.environ['DISABLE_PANDERA_IMPORT_WARNING'] = 'True'

# Standard library
import sys
import re
import glob
import json
import time
import argparse
import inspect
import subprocess
import shutil
from pathlib import Path
from datetime import datetime

# Third-party
try:
    import numpy as np
except Exception:
    np = None  # type: ignore

try:
    from tqdm import tqdm
except Exception:
    def tqdm(iterable, **kwargs):  # type: ignore
        for item in iterable:
            yield item


def get_chai_run_inference(engine: str):
    """Return the appropriate Chai run_inference callable for the selected engine."""
    engine_lower = (engine or "chai1").lower()
    if engine_lower == "chai2":
        try:
            from chai_lab.chai2 import run_inference  # type: ignore
            return run_inference
        except Exception as exc:
            raise ImportError(f"Chai-2 engine requested but unavailable: {exc}")
    else:
        from chai_lab.chai1 import run_inference  # type: ignore
        return run_inference


def parse_fasta_stem_to_dirname(fasta_stem: str) -> str:
    """Convert a FASTA stem to its output directory name.

    For the original R/ proteins the stem was 'ArRmut11_with_1a_and_PMP' → 'ArRmut11_1a'.
    For smuprot entries the stem is already the complex name (e.g. 'protein_mn_ligand'),
    so it is returned unchanged.
    """
    # Legacy pattern kept for compatibility if ever reused with original FASTA naming
    if "_with_" in fasta_stem and "_and_" in fasta_stem:
        parts = fasta_stem.split("_with_")
        if len(parts) == 2:
            protein = parts[0]
            substrate_part = re.sub(r"_and_\w+$", "", parts[1])
            return f"{protein}_{substrate_part}"
    return fasta_stem
