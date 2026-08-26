#!/bin/bash
# =============================================================================
# run.sh  —  Chai-1 structure prediction pipeline for R_smuprot
#
# Converts every smuprot_N.json in the current directory to a Chai-1 FASTA,
# then runs ULTRA-HIGH quality prediction (500 steps, 5 samples) on each one.
#
# Usage
# -----
#   ./run.sh                        # process all smuprot_*.json files
#   ./run.sh smuprot_0              # process smuprot_0.json only
#   ./run.sh smuprot_0 chai2        # use Chai-2 engine
#
# Outputs
# -------
#   smuprot_0/smuprot_0.fasta       generated FASTA (protein + ligand SMILES)
#   smuprot_0/<complex_name>.cif    best predicted structure (CIF format)
#   smuprot_0/<complex_name>.pdb    same structure in PDB format
#   smuprot_0/summary_*.txt         score statistics across the 5 samples
# =============================================================================

echo "======================================================================"
echo "  Chai-1 Prediction Pipeline  —  R_smuprot"
echo "======================================================================"

# Load Python module on HPC / NERSC
source /usr/share/lmod/lmod/init/bash 2>/dev/null
module load python/3.11 2>/dev/null
echo "Python: $(python3 --version)"

TARGET="${1:-all}"    # e.g. smuprot_0, smuprot_1, smuprot_2, or 'all'
ENGINE="${2:-chai1}"

echo "Target : $TARGET"
echo "Engine : $ENGINE"
echo ""

# ---------------------------------------------------------------------------
# Helper: run one JSON → FASTA → Chai-1
# ---------------------------------------------------------------------------
run_one() {
    local json_file="$1"         # e.g. smuprot_0.json
    local stem="${json_file%.json}"  # e.g. smuprot_0

    echo "----------------------------------------------------------------------"
    echo "  Processing $json_file"
    echo "----------------------------------------------------------------------"

    # Step 1: JSON → FASTA  (Mn is preserved; UFF warning is suppressed)
    echo "[1/2] Converting JSON to FASTA ..."
    python3 json_to_fasta.py --json "$json_file"
    if [ $? -ne 0 ]; then
        echo "  ERROR: json_to_fasta.py failed for $json_file"
        return 1
    fi

    # Step 2: Chai-1 batch prediction on the generated folder
    echo "[2/2] Running Chai-1 ($ENGINE) on folder: $stem ..."
    python3 run_batch.py --protein "$stem" --engine "$ENGINE"
    if [ $? -ne 0 ]; then
        echo "  ERROR: run_batch.py failed for $stem"
        return 1
    fi

    echo "  Done: $stem"
}

# ---------------------------------------------------------------------------
# Main dispatch
# ---------------------------------------------------------------------------
if [ "$TARGET" = "all" ]; then
    found=0
    for json in smuprot_*.json; do
        [ -f "$json" ] && found=1 && break
    done
    if [ "$found" -eq 0 ]; then
        echo "Error: No smuprot_*.json files found in $(pwd)"
        exit 1
    fi
    echo "Processing all smuprot_*.json files ..."
    echo ""
    for json in smuprot_*.json; do
        [ -f "$json" ] && run_one "$json"
    done
else
    json_file="${TARGET}.json"
    if [ ! -f "$json_file" ]; then
        echo "Error: $json_file not found in $(pwd)"
        exit 1
    fi
    run_one "$json_file"
fi

echo ""
echo "======================================================================"
echo "  All predictions complete."
echo "  CIF/PDB files are saved inside each smuprot_N/ folder."
echo "======================================================================"
