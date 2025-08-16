#!/bin/bash -l
#SBATCH --time=96:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --mem=64g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=you@example.org
#SBATCH --job-name=IQTREEArrays
#SBATCH -o iqtree_%j_%a.out
#SBATCH -e iqtree_%j_%a.err
#SBATCH -p YOUR_PARTITION
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --array=1-1000  # Adjust this to your actual list length

# Set working directory
cd /path/to/workdir/IQTREE
module purge
module load iqtree2/2.1.2   # adjust to your env/module name

# Input & Output paths
alignment_dir="/path/to/alignments"
output_dir="IQTREE_Results_array"
mkdir -p "$output_dir"

# File list
ALIST="${alignment_dir}/list.txt"
ALIGNMENT_FILE="$(sed "${SLURM_ARRAY_TASK_ID}q;d" "${ALIST}")"

# Output file prefix (no .iqtree extension)
output_file_name="$(basename "$ALIGNMENT_FILE" .fa)"
output_prefix="${output_dir}/${output_file_name}"

# Run IQ-TREE2 with ultrafast bootstrap and SH-aLRT support
iqtree -s "$ALIGNMENT_FILE" \
       -st AA \
       -m MFP \
       -bb 1000 \
       -alrt 1000 \
       -T AUTO \
       -pre "$output_prefix"

echo "IQ-TREE finished for $ALIGNMENT_FILE"
