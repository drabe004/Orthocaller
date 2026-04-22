#!/bin/bash -l
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8g
#SBATCH --tmp=8g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=you@university.edu
#SBATCH --job-name=orthocaller
#SBATCH --output=logs/orthocaller_%A_%a.out
#SBATCH --error=logsv/orthocaller_%A_%a.err
#SBATCH --array=1-4000

# ---------------------------
# Set working directory
# ---------------------------
cd /projects/your/projects

mkdir -p logs

# ---------------------------
# Load environment
# ---------------------------
# Load environment
module load compatibility/agate-centos7
module load conda
source activate orthocaller
echo "[DEBUG] SLURM_ARRAY_TASK_ID='${SLURM_ARRAY_TASK_ID}'"
# ---------------------------
# Define variables
# ---------------------------
OG="${SLURM_ARRAY_TASK_ID}_generax"

RECON_DIR="Path/to/your/GeneRax_Run_Results"
MAP_DIR="Path/to/your/maps_Directory"
FOCAL_DIR="Path/to/your//focal_species_files_Directory"
TREE_DIR="Path/to/your/pruned_trees_Directory"
OUT_DIR="Path/to/your/Orthocaller_output"
CAVEFISH="Path/to/your/Foreground_Species_List.txt"
BACKGROUND="Path/to/your/Background_Species_List.txt"

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
python ortho_caller_Version8_local_standard.py \
    -og "$OG" \
    -st "$SPECIES_TREE" \
    -rd "$RECON_DIR" \
    -md "$MAP_DIR" \
    -od "$OUT_DIR" \
    --cavefish_list "$CAVEFISH" \
    --background_list "$BACKGROUND" \
    --inparalog_strategy shortest \
    --focal_species_file "$FOCAL_FILE" \
    --long_branch_threshold 5 \  ##Reassign node as a duplication if branch is > 5x as long as the parent branch
    --short_branch_threshold 0.1 ##Reassign node as a speciation if branch is < .1 as long as the parent branch 
