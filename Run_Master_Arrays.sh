#!/bin/bash -l
#SBATCH --time=96:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=1g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=email@.edu
#SBATCH --job-name=masterarray
#SBATCH --account=mcgaughs
#SBATCH --output=orthocaller_masterarray.out
#SBATCH --error=orthocaller_masterarray.err


cd /path/to/orthocaller/directory

##Runs iterative arrays 3000 at a time (SLURM SYSTEMS USUALLY SET LIMITS FOR ARRAY JOBS ON SHARED HPC CLUSTERS) %120 uses 120 cores (an entire node), so will run 120 orthogroups at a time. Adjust as needed based on computing resources available. 

sbatch --array=1-3000%120 --wait orthocaller-v3.1.sh
sbatch --array=3001-6000%120 --wait orthocaller-v3.1.sh
sbatch --array=6001-10000%120 --wait orthocaller-v3.1.sh


###etc... add as many lines as needed 