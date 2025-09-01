#!/bin/bash -l
#SBATCH --time=96:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8g
#SBATCH --tmp=8g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --job-name=gen_focalspeciesfile
#SBATCH -p 

module load conda 
source activate orthocaller




python generate_focal_species_from_maps.py -m maps -o focal_species_files
