# 🧬 Orthocaller v2.0  
**A Reconciliation-Based Pipeline for Evolutionary Orthology Assignment**  
_Adapted from [Comparative Genomics Snakes Toolkit](https://github.com/masonaj157/Comparative_genomics_snakes)_

---
# 🧬 Orthocaller v3.0

**Orthocaller v3** is a gene tree-based orthogroup classification tool for comparative genomics. It uses reconciled gene trees (e.g., from **GeneRax**) and a rooted species tree to identify conserved, duplicated, and lost orthologs across species.

---

## 🚨 What's New in v3.0?

This version introduces **flexible in-paralog resolution**, with three user-selectable strategies for deciding which in-paralog copy to retain:

| Strategy            | Description                                                                 |
|---------------------|-----------------------------------------------------------------------------|
| `farthest` (default) | Retains the most diverged in-paralog (longest branch from duplication node). |
| `shortest`          | Retains the least diverged in-paralog (shortest branch).                     |
| `average_divergence`| Retains the copy whose branch length ratio is closest to the duplication node — a proxy for balanced divergence. |



---
# 🧬 Orthocaller v3.1
---

## 🚨 What's New in v3.1?
1) Improved branch length handling
The branch length options remain functional, but very short branches occasionally caused problems when re-rooting. Version 3.1 introduces a safeguard: if the ancestral node is already the root, the algorithm now skips re-rooting and grafts the subtree directly. This prevents errors associated with near-zero branch lengths.

2) Updated classification logic for in-paralog pruning
After pruning in-paralogs at SD nodes, the grafted node is now reclassified from a duplication (“D”) to a speciation (“S”) event. This change reflects the biological reality that once redundant in-paralogs are removed, the remaining structure represents speciation. Retaining the “D” label in this context was misleading and inconsistent with the model.

3) Safeguard against reclassification of artificial short branches
Re-grafted subtrees introduced by in-paralog pruning are often extremely short by design. These branches are now excluded from further branch-length–based reclassification to avoid artifacts, ensuring that only biologically meaningful nodes are considered in downstream analyses.

Together, these updates provide more consistent handling of tree re-rooting and classification, especially in cases involving short branches and in-paralog pruning.

Why this matters
These updates ensure that tree handling and classification better reflect the underlying biology. By treating pruned in-paralogs as true speciation events and preventing artifacts from very short branches, Version 3.1 improves the accuracy and interpretability of downstream analyses such as selection tests. In short: results are now both more robust and more biologically meaningful.

---

---
# 🧬 Orthocaller v3.2
---

Long branch thresholds and Short branch thresholds are passable now as arguments 

---

## 📖 Overview
**Orthocaller v2.0 +** classifies orthogroups using reconciled gene trees (e.g., from **GeneRax**), a rooted species tree, and a gene-to-species map.  
It identifies evolutionary events like:
- ✅ Conserved orthologs  
- 🔁 Species-specific duplications  
- ❌ Gene losses  
- 🧩 Paralogy  

It outputs ortholog group tables, annotated trees, and classification summaries—ideal for comparative genomics and evolutionary analyses.

## Description
This script classifies orthogroups using reconciled gene trees (e.g., from GeneRax), a rooted species tree, and a gene-to-species map. It identifies evolutionary events such as conserved orthologs, species-specific duplications, gene losses, and paralogy, and outputs ortholog group tables, class summaries, and NHX/NWK tree files for downstream comparative genomic and functional analyses.Designed for evolutionary biologists, genome scientists, and comparative physiologists, this tool is optimized for large-scale, high-throughput orthogroup classification, particularly in cases where deep divergence, gene duplication, or gene loss might obscure one-to-one orthology relationships.Use case: This pipeline was developed to study convergent and lineage-specific patterns of gene evolution across cavefish species and background lineages—distinguishing true orthologs from confounding paralogs or loss events using tree topology and evolutionary events inferred from reconciliation.

## Why Reconciliation?
Orthology inference based solely on clustering or sequence similarity (e.g., OrthoFinder or BLAST) often misses the mark in non-model or deeply diverged species. Reconciliation-based methods incorporate both gene trees and species trees, enabling better resolution of:
-In-paralogs (species-specific duplications)
-False orthology (where topology or branch lengths suggest a misclassification)
-Gene loss or gain events along specific branches
-This script wraps around reconciled trees from GeneRax and builds downstream annotations and ortholog tables.

## 📥 Input Requirements

| Argument | Description |
|----------|-------------|
| `-og` | Orthogroup name (e.g., `1_generax`) |
| `-st` | Rooted species tree in Newick format |
| `-rd` | Directory with GeneRax reconciled output |
| `-md` | Folder with gene-to-species `.map` files (e.g., Species:Seq1;Seq2).|
| `-od` | Output directory |
| `--cavefish_list` | Text file with cavefish species (aka foreground species of interest) (1 per line) |
| `--background_list` | Text file with background species (aka background species of interest)  |
| '--focal_species_file' | a text file listing all the species in the orthogroup |
| '--inparalog_strategy' | How do you want the script to choose between paralogs at the tree tips? Choices--- 1) shortest (chooses the paralog with the shortest branch from the D node) 2) average_distance (chooses the paralog with a branch length that is the clostest to the branch length of the parent branch leading to the D node) 3) farthest (chooses the paralog with the longest branch length from the parent D node) |


## 📤 Output Files per Orthogroup

| File | Purpose |
|------|---------|
| `*.nhx` | Annotated gene tree with species + evolutionary events |
| `*_decostar.nwk` | Pruned tree for DeCoSTAR input |
| `*_codeml.nwk` | Codeml-compatible tree |
| `*_expr.nhx` | Tree for expression overlay |
| `*_pruned_species.nwk` | Species tree with focal species only |
| `*_orthogroup.csv` | Table of orthologous genes |
| `*_genes.csv` | Gene-to-orthogroup key |
| `*_classes.csv` | Classifications per gene group |
| `summary.txt` | High-level stats & summary |

### Running the Script
You can run the script manually:
python orthocaller-v3.1.py \
    -og 1_generax \
    -st SpeciesTree_rooted.tre \
    -rd GeneRax_Output_Dir \
    -md maps/ \
    -od 5.1_OG_classifications/ \
    --cavefish_list Cavefish_List.txt \
    --background_list BackgroundFish_List.txt
    --focal_species_file focalfile.txt \
    --inparalog_strategy shortest 

Or at scale using the included SLURM script (ortho_caller_Version2_Official.sh), which parallelizes across 3000 orthogroups using SLURM job arrays:
sbatch orthocaller-v3.1.sh 

This SLURM script:
Submits 3000 jobs (--array=1-3000)
Logs failed orthogroups for reruns
Scans logs to generate a list of failed jobs (failed_jobs_from_logs.txt)

## ⚙️ Dependencies

- Python ≥ 3.8  
- `ete3`  
- `pandas`  

Create an environment:

```bash
conda create -n orthocaller python=3.10 ete3 pandas
conda activate orthocaller
```

## Evolutionary Reasoning
This pipeline aims to improve orthogroup classification by integrating information from both gene and species trees in a biologically informed way. Gene trees alone can sometimes be misleading—topological artifacts, missing data, or uneven evolutionary rates can obscure which genes are truly orthologous. By examining branch lengths and tree structure, we flag cases where one gene copy appears to evolve much faster than its sibling, which may indicate misassignment or unresolved duplication. We also recognize that not every inferred duplication event reflects true evolutionary history; some may be artifacts of limited phylogenetic resolution, so we apply conservative filters to reduce false positives. The fundamental assumption is that true orthologs will appear once per focal species. When this pattern is broken—due to missing species or extra copies—it suggests gene loss or duplication, which we classify accordingly. In doing so, this pipeline provides a more nuanced and testable framework for understanding

## Dependencies
Python ≥ 3.8
ete3
pandas


## Install dependencies via conda:
conda create -n orthocaller python=3.10 ete3 pandas
conda activate orthocaller

## Citation
Drabeck & Mason, 2024
Orthocaller Version 2: A branch-length-aware orthogroup classification pipeline
Part of the Comparative Genomics Snake Toolkit

## Orthocaller Version 2: A branch-length-aware orthogroup classification pipeline
Adapted from the Comparative Genomics Snake Toolkit
Available at: https://github.com/masonaj157/Comparative_genomics_snakes

Orthocaller v2 was developed by Andrew J. Mason (cit) and adapted by Danielle H Drabeck to address key challenges in distinguishing orthology, paralogy, and gene loss in large comparative datasets. The pipeline integrates reconciled gene trees, a reference species tree, and curated species mappings to classify orthogroups with sensitivity to topological uncertainty and branch length asymmetry.


Questions? Bugs?
Open an issue or contact:
drabe004@umn.edu


## 🛠 Accessory Scripts

The following two accessory scripts assist with filtering and extracting gene-specific sequence alignments based on species representation in orthogroups. These are especially useful for downstream analysis when focusing on genes with broad taxonomic sampling.


## Preparing input from Generax 
# Orthocaller – Accessory Pipeline

This accessory pipeline uses **GeneRax output** to prepare all the necessary inputs for **Orthocaller**.  
It automates the generation of species maps, focal species files, and pruned species trees derived from GeneRax reconciliations.  

---

## Workflow Overview

1. **Generate species maps** – Parse GeneRax reconciled trees to create mapping files linking sequence IDs to species names.  
2. **Generate focal species files** – Use species maps to define per-orthogroup focal species sets.  
3. **Generate pruned species trees** – Build reduced species trees from the focal species lists, suitable for Orthocaller.  

The outputs of this pipeline are used directly by Orthocaller.

---

## Step 1: Generate Species Maps

Run the loop script to write species maps from GeneRax reconciliations:

1) Make map files 
generate_map_loopz2.sh

Step 2: Generate Focal Species Files from map files 

generate_focal_species_from_maps.py
generate_focal_species_from_maps.sh


Step 3: Generate Pruned Species Trees
Prune the master species tree to reflect only focal species, producing per-orthogroup trees for Orthocaller:

writeTrees.py
WriteTrees.sh

Notes
These scripts assume completed GeneRax runs with reconciled gene trees.

Outputs include:

*.map files (species ↔ sequence mappings)
focal species files (plain text)
pruned species trees (Newick format)


## Running Orthocaller arrays on HPC clusters

### 1. 'Run_Master_Arrays.sh'
--allows a user to run thousands of orthogroups at once, iterative waiting to acomodate shared HPC cluster array job limits (ours is 3k). 

## Processing output from orthocaller
### 1. `FilterSummaryFiles.py`

**Purpose**  
Scans the `summary.txt` file for each orthogroup and extracts only those gene entries that meet user-defined thresholds for the number of cavefish and background species. This is useful for identifying well-sampled gene trees for comparative or statistical analysis.

**Inputs**
- A directory containing multiple `summary.txt` files (e.g., `1_generax/summary.txt`, `2_generax/summary.txt`, ...).
- Thresholds for minimum number of cavefish and background species (default: ≥8 cavefish, ≥31 background).

**Output**
- A single master summary file containing only gene lines that pass the filtering criteria (e.g., `master_summary_file_7_OG_classifications.txt`).

**Example**
```bash
python FilterSummaryFiles.py \
  -i /path/to/orthogroup_directories/ \
  -o master_summary_file_7_OG_classifications.txt \
  --min-cavefish 8 \
  --min-background 31

### 2. `GetAlignmentsFromOrthocallerOutput.py`

**Purpose**  
Takes a master summary file (e.g., generated by `FilterSummaryFiles.py`) and extracts the corresponding gene-specific sequences from precomputed alignment FASTA files. It maps orthogroups to alignments using a `GeneRaxKey.txt` file and uses `.csv` metadata to identify the relevant sequences per gene.

**Inputs**
- Master summary file with filtered gene entries
- Directory of orthogroup classification `.csv` files
- Directory of alignment FASTA files (GeneRax-formatted)
- `GeneRaxKey.txt` file mapping orthogroup names to alignment filenames

**Output**
- Individual FASTA files for each gene, written to the specified output directory  
  Each file contains only the sequences corresponding to that gene.

**Example**
```bash
python GetAlignmentsFromOrthocallerOutput.py \
  -s master_summary_file_7_OG_classifications.txt \
  -o ./output_alignments/ \
  -b ./7_OG_classifications/ \
  -a ./GeneraxFormattedAlns2/ \
  -k GeneRaxKey.txt
