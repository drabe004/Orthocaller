#!/bin/bash -l        
#SBATCH --time=96:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8g
#SBATCH --tmp=10g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=


c
module load conda 
source activate orthocaller

# Set paths
MASTER_TREE="/mastertree.tre"
FOCAL_DIR="path/to/focal_species_files"
OUTPUT_DIR="pruned_trees"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Run the pruning script (positional arguments only!)
python writeTrees.py "$MASTER_TREE" "$FOCAL_DIR" "$OUTPUT_DIR"
