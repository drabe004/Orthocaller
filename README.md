Orthocaller-v2.0
Orthocaller v2.0 — A Reconciliation-Based Pipeline for Evolutionary Orthology Assignment
This is a modified version of the orthocaller.py pipeline from https://github.com/masonaj157/Comparative_genomics_snakes adapted for broad use on large orthogroup data. 

Overview
This script classifies orthogroups using reconciled gene trees (e.g., from GeneRax), a rooted species tree, and a gene-to-species map. It identifies evolutionary events such as conserved orthologs, species-specific duplications, gene losses, and paralogy, and outputs ortholog group tables, class summaries, and NHX/NWK tree files for downstream comparative genomic and functional analyses.

Designed for evolutionary biologists, genome scientists, and comparative physiologists, this tool is optimized for large-scale, high-throughput orthogroup classification, particularly in cases where deep divergence, gene duplication, or gene loss might obscure one-to-one orthology relationships.

Use case: This pipeline was developed to study convergent and lineage-specific patterns of gene evolution across cavefish species and background lineages—distinguishing true orthologs from confounding paralogs or loss events using tree topology and evolutionary events inferred from reconciliation.

Why Reconciliation?
Orthology inference based solely on clustering or sequence similarity (e.g., OrthoFinder or BLAST) often misses the mark in non-model or deeply diverged species. Reconciliation-based methods incorporate both gene trees and species trees, enabling better resolution of:

In-paralogs (species-specific duplications)

False orthology (where topology or branch lengths suggest a misclassification)

Gene loss or gain events along specific branches

This script wraps around reconciled trees from GeneRax and builds downstream annotations and ortholog tables.

Input Requirements
Orthogroup name (-og): Typically 1_generax, 2_generax, etc. This will be used to locate the reconciled tree directory.

Species tree (-st): A rooted Newick-format species tree with consistent naming.

Gene tree directory (-rd): Directory with GeneRax output. Assumes each orthogroup has a subdirectory with reconciled trees in it.

Map directory (-md): Folder with .map files that relate gene sequence names to species names (e.g., Species:Seq1;Seq2).

Output directory (-od): Where the classification results for each orthogroup are written.

Cavefish species list (--cavefish_list): Plain-text file of cavefish species names (one per line).

Background species list (--background_list): Plain-text file of comparison species.

Outputs
Each orthogroup will produce the following in its own output folder:

File	Description
*.nhx	NHX-format gene tree annotated with species and evolutionary event types (S, D, etc.)
*_decostar.nwk	Pruned Newick tree formatted for DeCoSTAR
*_codeml.nwk	Codeml-compatible Newick tree
*_expr.nhx	Expression-ready tree for downstream annotation
*_pruned_species.nwk	Pruned species tree containing only focal species
*_orthogroup.csv	Final table of orthologous genes by group
*_genes.csv	Key of each gene to its assigned orthogroup
*_classes.csv	Classification table showing conserved, loss, or duplication events for each paralog group
summary.txt	High-level description of species coverage, event types, and counts for this orthogroup

Running the Script
You can run the script manually:
python ortho_caller_Version2_Official.py \
    -og 1_generax \
    -st SpeciesTree_rooted.tre \
    -rd GeneRax_Output_Dir \
    -md maps/ \
    -od 5.1_OG_classifications/ \
    --cavefish_list Cavefish_List.txt \
    --background_list BackgroundFish_List.txt

Or at scale using the included SLURM script (ortho_caller_Version2_Official.sh), which parallelizes across 3000 orthogroups using SLURM job arrays:
sbatch ortho_caller_Version2_Official.sh

This SLURM script:
Submits 3000 jobs (--array=1-3000)
Logs failed orthogroups for reruns
Scans logs to generate a list of failed jobs (failed_jobs_from_logs.txt)

Evolutionary Reasoning: 
This pipeline aims to improve orthogroup classification by integrating information from both gene and species trees in a biologically informed way. Gene trees alone can sometimes be misleading—topological artifacts, missing data, or uneven evolutionary rates can obscure which genes are truly orthologous. By examining branch lengths and tree structure, we flag cases where one gene copy appears to evolve much faster than its sibling, which may indicate misassignment or unresolved duplication. We also recognize that not every inferred duplication event reflects true evolutionary history; some may be artifacts of limited phylogenetic resolution, so we apply conservative filters to reduce false positives. The fundamental assumption is that true orthologs will appear once per focal species. When this pattern is broken—due to missing species or extra copies—it suggests gene loss or duplication, which we classify accordingly. In doing so, this pipeline provides a more nuanced and testable framework for understanding

Dependencies
Python ≥ 3.8
ete3
pandas


Install dependencies via conda:
conda create -n orthocaller python=3.10 ete3 pandas
conda activate orthocaller

Citation
Drabeck & Mason, 2024
Orthocaller Version 2: A branch-length-aware orthogroup classification pipeline
Part of the Comparative Genomics Snake Toolkit

Orthocaller Version 2: A branch-length-aware orthogroup classification pipeline
Part of the Comparative Genomics Snake Toolkit
Available at: https://github.com/masonaj157/Comparative_genomics_snakes

Orthocaller v2 was developed by Andrew J. Mason (cit) and adapted by Danielle H Drabeck to address key challenges in distinguishing orthology, paralogy, and gene loss in large comparative datasets. The pipeline integrates reconciled gene trees, a reference species tree, and curated species mappings to classify orthogroups with sensitivity to topological uncertainty and branch length asymmetry.


Questions? Bugs?
Open an issue or contact:
drabe004@umn.edu


