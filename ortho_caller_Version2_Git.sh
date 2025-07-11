#!/bin/bash -l
#SBATCH --time=400:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8g
#SBATCH --tmp=8g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --job-name=orthocaller
#SBATCH -p privatenode
#SBATCH --output=logs/orthocaller_%A_%a.out
#SBATCH --error=logs/orthocaller_%A_%a.err
#SBATCH --array=1-3000%120  # Run up to 120 jobs in parallel

# Set working directory and create logs dir
cd your/paths/to/home/dir
mkdir -p logs

# Load environment
module load compatibility/agate-centos7
module load conda
source activate orthocaller

# Variables
OG=${SLURM_ARRAY_TASK_ID}_generax
SPECIES_TREE="your/paths/to/home/dir/SpeciesTree.tre"
RECON_DIR="your/paths/to/home/dir/GeneRax"
OUT_DIR="your/paths/to/home/dir/5.1_OG_classifications"
CAVEFISH="your/paths/to/home/dir/Species_List1.txt"
BACKGROUND="your/paths/to/home/dir/Specieslist2.txt"
MAP_DIR="your/paths/to/home/dir/maps"

# Run the script
python ortho_caller_Version2_Official.py \
    -og "$OG" \
    -st "$SPECIES_TREE" \
    -rd "$RECON_DIR" \
    -md "$MAP_DIR" \
    -od "$OUT_DIR" \
    --cavefish_list "$CAVEFISH" \
    --background_list "$BACKGROUND"

# Log failed orthogroups (if this job errored)
if [ $? -ne 0 ]; then
    echo "$OG" >> failed_orthogroups_runtime.txt
fi

# Final summary on the last array task
if [ "$SLURM_ARRAY_TASK_ID" -eq 3000 ]; then
    echo "Running failure summary scan from task $SLURM_ARRAY_TASK_ID..."
    find logs -name "orthocaller_*.err" -size +0 -exec basename {} .err \; | \
        sed 's/orthocaller_\([0-9]*\)_[0-9]*$/\1_generax/' | \
        sort -n > failed_jobs_from_logs.txt
    echo "Finished writing failed_jobs_from_logs.txt"
fi
