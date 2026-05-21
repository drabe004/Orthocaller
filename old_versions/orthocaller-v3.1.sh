#!/bin/bash -l
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8g
#SBATCH --tmp=8g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=email@blank.edu
#SBATCH --job-name=orthocaller
#SBATCH --output=logs/orthocaller_%A_%a.out ###CHANGE THIS LOGS DIR TO RUN WITHOUT OVERWRITING LOGS
#SBATCH --error=logs/orthocaller_%A_%a.err  ###CHANGE THIS LOGS DIR TO RUN WITHOUT OVERWRITING LOGS



# ---------------------------
# Set working directory
# ---------------------------
cd path/to/your/Orthocaller3.1
mkdir -p logs   ###CHANGE THIS LOGS DIR TO RUN WITHOUT OVERWRITING LOGS



# ---------------------------
# Load environment
# ---------------------------
# Load environment
module load compatibility/agate-centos7
module load conda
source activate orthocaller

# ---------------------------
# Define variables
# ---------------------------
OG="${SLURM_ARRAY_TASK_ID}_generax"

RECON_DIR="/path/to/GeneRax_Run_Results/"
MAP_DIR="/path/to/maps"
FOCAL_DIR="path/to/focal_species_files"
TREE_DIR="path/to/pruned_trees"
OUT_DIR="path/to/output/dir"
CAVEFISH="path/to/Cavefish_List.txt"
BACKGROUND="path/to/BackgroundFish_List.txt"

FOCAL_FILE="${FOCAL_DIR}/${OG}.focal.txt"
SPECIES_TREE="${TREE_DIR}/${OG}.tre"

# ---------------------------
# Sanity checks
# ---------------------------
if [[ ! -f "$SPECIES_TREE" ]]; then
    echo "[ERROR] Missing species tree: $SPECIES_TREE"
    exit 1
fi

if [[ ! -f "$FOCAL_FILE" ]]; then
    echo "[ERROR] Missing focal species file: $FOCAL_FILE"
    exit 1
fi

if [[ ! -f "$CAVEFISH" || ! -f "$BACKGROUND" ]]; then
    echo "[ERROR] Cavefish or background list not found."
    exit 1
fi

echo "[INFO] Running orthocaller on $OG"
echo "[INFO] Species tree: $SPECIES_TREE"
echo "[INFO] Focal species file: $FOCAL_FILE"
echo "[INFO] Output directory: ${OUT_DIR}/${OG}"

# ---------------------------
# Run OrthoCaller
# ---------------------------
python orthocaller-v3.1.py \
    -og "$OG" \
    -st "$SPECIES_TREE" \
    -rd "$RECON_DIR" \
    -md "$MAP_DIR" \
    -od "$OUT_DIR" \
    --cavefish_list "$CAVEFISH" \
    --background_list "$BACKGROUND" \
    --inparalog_strategy farthest \
    --focal_species_file "$FOCAL_FILE"
