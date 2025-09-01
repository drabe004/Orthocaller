#!/bin/bash -l
#SBATCH --time=96:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8g
#SBATCH --tmp=8g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --job-name=gen_map




module load compatibility/agate-centos7
module load conda 
source activate orthocaller



species_list="Species_List.txt"
mapdir="maps"

mkdir -p "$mapdir"

for newick in /GeneRax_Run_Results/*_generax/reconciliations/*.newick
do
    orthogroup=$(basename $(dirname $(dirname $newick)))   # e.g. 1_generax

    outmap="${mapdir}/${orthogroup}.map"

    echo "Generating $outmap from $newick"

    python generate_mapsv2.py \
       -n "$newick" \
       -s "$species_list" \
       -o "$outmap"
done