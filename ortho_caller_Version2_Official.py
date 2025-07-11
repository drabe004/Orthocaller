#!/usr/bin/env python

# Additional software necessary to run this:
# -*- coding: utf-8 -*-

# (1) biopython
# (2) pandas
# (3) dfply
# (4) mafft
# (5) iqtree
# (6) generax v.XXX

import glob
import argparse
import copy
from re import sub
import threading
from ete3 import Tree
from ete3 import PhyloTree
from itertools import combinations
import subprocess as sp
import pandas as pd
import csv
#from ete3 import Phyloxml, phyloxml

from ctypes import alignment
import sys, os, shutil
import datetime as dt
#import numpy as np
from dfply import *


########################################
############### ARGUMENTS ##############
########################################

# ---------------------------
# Argument Parser
# ---------------------------
# This script classifies orthogroups based on reconciled gene trees (e.g., from GeneRax),
# using a species tree and a mapping of sequences to species.
# It also uses cavefish vs. background species lists to track evolutionary events.
# ---------------------------


parser = argparse.ArgumentParser(description='')
parser.add_argument("-og", "--orthogroup", type=str, required=True, help="Classify orthogroups using reconciled gene trees, species tree, and sequence maps.")

parser.add_argument("-st","--species_tree",
					type=str,
                    default="Classify orthogroups using reconciled gene trees, species tree, and sequence maps.",
					help=(
        "Species tree in Newick format, used as the guide/reference tree for reconciliation.\n"
        "Tip: This must be rooted and use species names consistent with your map files and sequence headers."
    )
parser.add_argument("-rd","--reconciled_dir",
					type=str,
                    default="3.5_OG_reconciliations",
					help="Directory containing generax input and output files")
parser.add_argument("-md","--map_dir",
					type=str,
                    default="3.1_Species_maps",
					help= help=(
        "Directory containing map files (Species:Species_geneID1; Species_geneID2) 1_generax.map with columns:\n"
        "    Sequence_ID, Species_Name\n"
        "Used to map gene tree tips to species names."
    )
parser.add_argument("-od","--output_dir",
					type=str,
                    default="3.6_OG_classifications",
					help="Directory of where to write a folder for this orthogroup and subsequent output files")
parser.add_argument("--cavefish_list",
                    type=str,
                    default="Cavefish_List.txt",
                    help="File with cavefish species names")
parser.add_argument("--background_list",
                    type=str,
                    default="BackgroundFish_List.txt",
                    help="File with background species names")
args = parser.parse_args()

########################################
################# SETUP ################
########################################
# Load cavefish and background species lists
with open(args.cavefish_list) as f:
    cavefish_list = [line.strip() for line in f if line.strip()]

with open(args.background_list) as f:
    background_list = [line.strip() for line in f if line.strip()]


species_map_df = None  # Global variable for the species map

orthogroup = args.orthogroup
reconciled_dir = args.reconciled_dir
map_dir = args.map_dir
species_tree = args.species_tree
output_dir = args.output_dir


########################################
############## FUNCTIONS ###############
########################################

# Returns a list of species names by stripping file extensions from all filenames in a directory, replaced this with the global_species_lookup function
#def find_species(genomes_dir) :
#	species_list = list(os.listdir(genomes_dir))
#	species_list = [x.split('.')[0] for x in species_list]
#	return(species_list)


###################Danielle edited some of this to be able to take a wildcard file name for looped runs
# Looks up the species name for a leaf node based on its name using a prefix match to the species map

def global_species_lookup(node):
    name = node.name

    if not node.is_leaf():
        return None  # Don't try to assign species to internal nodes

    if name is None:
        return None

    # Extract species name
    if "_" in name:
        species_prefix = "_".join(name.split("_")[:-1])
    else:
        species_prefix = name

    match = species_map_df[species_map_df["New_name"] == species_prefix]
    if match.empty:
        raise ValueError(f"No species match in map for: {name}")
    return species_prefix

#################################THIS FUNCTION IS SUPER IMPORTANT AND IF YOUR HAVING NAMING ERRORS ITS PROBABLY HERE####################
# Loads the reconciled gene tree for an orthogroup from the GeneRax output, verifies uniqueness of gene names, and prints basic diagnostics
###############################################################################################
def read_and_reconcile_tree(orthogroup, reconciled_dir):
    import glob
    tree_dir = os.path.join(reconciled_dir, orthogroup, "reconciliations")
    tree_candidates = glob.glob(os.path.join(tree_dir, "*events.newick"))

    if not tree_candidates:
        raise FileNotFoundError(f"No events.newick file found in {tree_dir}")

    tree_path = tree_candidates[0]
    tree = PhyloTree(tree_path, format=1)

    # DEBUG: Print first 10 leaf names
    print("[DEBUG] First 10 leaf names from gene tree:")
    for leaf in list(tree.iter_leaves())[:10]:
        print(leaf.name)

    # DEBUG: Print map head
    print("[DEBUG] First 5 entries in map file:")
    print(species_map_df.head())

    # Check for duplicate gene names
    leaf_names = [leaf.name for leaf in tree.iter_leaves()]
    duplicates = set(name for name in leaf_names if leaf_names.count(name) > 1)
    if duplicates:
        print(f"[ERROR] Duplicate gene leaf names in tree: {duplicates}")
        raise ValueError("Non-unique gene names in tree.")

    # Use the global top-level species naming function (avoids pickling errors)
    tree.set_species_naming_function(global_species_lookup)

    # Trigger species parsing for sanity check
    for leaf in tree.iter_leaves():
        try:
            _ = leaf.species
        except Exception as e:
            print(f"[WARN] Failed to assign species for leaf: {leaf.name} - {e}")

    # Reconcile
    recon_tree, events = tree.reconcile(st)
    return recon_tree
##########################################End major DD edits #######################################
######################################################################################################

# Applies classification steps to the gene tree: resolves topology issues, defines paralog groups, and filters to focal species
def modify_classifications(tree, st, focal_species):
    tree = modify_in_paralog(tree)
    identify_and_resolve_topology_problems(tree,st)
    identify_putative_false_orthology(tree)
    final_paralogs_list,nhx_tree = fix_in_paralogs_define_groups(tree, focal_species)
    final_paralogs_list = filter_nonfocal_species(final_paralogs_list, focal_species)
    return(final_paralogs_list, nhx_tree)

# Identifies child speciation nodes (S) beneath duplication nodes (D) that have more than one species and unequal branch lengths, and reclassifies them as duplications (D).
# Corrects false orthology assignments by detecting asymmetric tree topologies where an 'S' node has two children with large branch length imbalance, inconsistent with true speciation.
def identify_putative_false_orthology(tree):
    for node in tree.traverse():
        if (hasattr(node,'evoltype') and node.evoltype == 'D'):
            for child in node.get_children():
                if hasattr(child,'evoltype') and child.evoltype == 'S' and len(list(child.get_species())) != 1:
                    ref = child.dist
                    print(ref)
                    ## get_children means we are only looking at two nodes. Will also speed things up.
                    for descendent in child.get_children():
                        if long_br(descendent,ref) and descendent.up.evoltype != 'SD':
                            print('')
                            print("ref ", ref)
                            print("made a change based on ")
                            print(descendent)
                            descendent.up.evoltype = 'D'
                            print(descendent.up)
                            #print('')

# Calculates the ratio of a descendent’s branch length to its sibling’s; returns True if the descendent’s branch is >5× longer, excluding zero-lengths via pseudovalues.
# This flags extreme branch asymmetry at a speciation (S) node under a duplication (D), where one child likely represents a fast-evolving or misassigned paralog.
# Reclassifies misidentified speciation (S) nodes as duplications (D) if one child branch is >5× longer than the other,
# suggesting the pair are likely paralogs due to asymmetric divergence, not true orthologs.
def long_br(descendent,ref):
    if ref == 0.0:
        ref = 0.000001
    test = descendent.dist
    if test == 0.0:
        test = 0.000001
    print("test", test)
    val = test /ref
    print("val", val)
    if val > 5 and ref != 0.000001:
        return True
    else:
        return False

# Looks for duplication (D) nodes directly above speciation (S) nodes where the branch is very short and species overlap,
# suggesting incorrect topology — then reroots or rearranges the tree to fix likely misplacements.
def identify_and_resolve_topology_problems(tree,st):
    for node in tree.traverse():
        if hasattr(node, 'evoltype') and node.evoltype == 'D' :
            for child in node.get_children():
                if hasattr(child,'evoltype') and child.evoltype == 'S' :
                    #print("running test")
                    if short_br(node,child) and species_rep_test(node):
                        #node_to_fix = node.detach()
                        #print(node)
                        fix_node(tree,node,st)

###################################################################################################################################################
 ###################################### Determine whether the branch leading to a child node is unusually short relative to its parent.
#####################################################################################################################################################
    # Args:
    #     node (TreeNode): The parent node in the tree.
    #     child (TreeNode): The child node to compare.
    #
    # Returns:
    #     bool: True if the child branch is less than 10% the length of the parent's branch, False otherwise.
    #
    # Notes:
    #     - Very short branches can indicate poor phylogenetic resolution or redundancy.
    #     - To avoid division by zero errors, extremely small values are substituted for zero-length branches.

def short_br(node,child):
    ref = node.dist
    if ref == 0.0:
        ref = 0.00000000001
    test = node.get_distance(child)
    if test == 0.0:
        test = test = 0.00000000001
    val = test /ref
    #print("parent is ", str(ref))
    #print("child is ",str(test))
    #print("test ratio ", str(val))
    if val < 0.1:
        return True
    else:
        return False

################################################################################################################################################################
# ============================================
# Step: Validate that each species is uniquely represented
#        (excluding in-paralogs from duplication events)
# ============================================
# Tests whether a node contains no more than one representative per species,
# accounting for species with in-paralogs (i.e., duplication events labeled as "SD").
#
# Returns:
#     True if, after excluding in-paralog descendants, each species appears only once.
#     False if any species appears more than once.
#
# Args:
#     node (TreeNode): A node in a gene tree with annotated species and duplication events.
#
# Notes:
#     - Assumes that in-paralogs are identified by internal nodes with evoltype == "SD".
#     - Traverses all descendants of the node to remove species associated with in-paralog duplications.def species_rep_test(node):
################################################################################################################################################################
def species_rep_test(node):
    species = [x.species for x in node.get_leaves()]
    sd_species = [x.species for x in node.traverse() if hasattr(node,"evoltype") and node.evoltype == "SD"]
    if len(sd_species) != 0:
        for i in sd_species:
            species.remove(i)
    if any([species.count(x) > 1 for x in species]):
        return False
    else:
        return True
##############################################################################################################################################################
# Replaces a problematic subtree (bad_node) in a gene tree with a re-rooted and reconciled version,
# using the species tree to find an appropriate outgroup and reassign evolutionary events.
#
# Args:
#     tree (PhyloTree): The full gene tree (not directly modified here but useful for context if needed).
#     bad_node (TreeNode): The problematic subtree suspected of incorrect rooting or structure.
#     st (PhyloTree): The species tree used for reconciliation and identifying outgroups.
#
# Steps:
#     1. Deep copy the problematic node and species tree to preserve original structures.
#     2. Prune the species tree to only include species represented in the subtree.
#     3. Identify a suitable outgroup using `find_reference_outgroup()`.
#     4. If multiple instances of the outgroup species are present (e.g., paralogs), re-root using the ancestor of those sequences.
#        Otherwise, re-root using a single leaf.
#     5. Reconcile the re-rooted subtree with the pruned species tree to update duplication/speciation events.
#     6. Replace the old subtree in the original tree with the corrected version.
##################################################################################################################################################################
def fix_node(tree, bad_node, st):
    node = copy.deepcopy(bad_node)
    spare_st = copy.deepcopy(st)

    # Get species represented in the subtree
    represented_species = [x for x in spare_st.get_leaf_names() if x in node.get_species()]
    spare_st.prune(represented_species)

    # Identify appropriate outgroup
    outgroup_species = find_reference_outgroup(spare_st)
    species_instances = [x.species for x in node.get_leaves()]

    # Re-root by common ancestor if outgroup is duplicated; otherwise by leaf
    if len(outgroup_species) > 1 or species_instances.count(outgroup_species[0]) > 1:
        species_leaves = search_for_species(node, outgroup_species)
        out_anc = node.get_common_ancestor(species_leaves)
        node.set_outgroup(out_anc)
    else:
        out_leaf = search_for_species(node, outgroup_species)[0]
        node.set_outgroup(out_leaf)

    # Reconcile to assign duplication/speciation events
    rec_node, node_events = node.reconcile(spare_st)

    # Replace old node with fixed one
    bad_node.up.add_child(node)
    bad_node.detach()
    

# Searches the gene tree for all leaf nodes whose assigned species match entries in the input list;
# used to locate where focal species occur in the gene tree for downstream classification or filtering.
def search_for_species(tree, species):
    node_list = []
    for entry in species:
        for node in tree.traverse():
            if hasattr(node,"species") and node.species == entry:
                node_list.append(node)
    return(node_list)

# Identifies a reference outgroup from the species tree by selecting the smallest clade (usually the most basal lineage).
def find_reference_outgroup(st):
    outgroup = list(st.get_leaf_names())
    for child in st.get_children():
        if child.is_leaf() == True:
            outgroup = [child.name]
            break
        elif len(list(child.get_leaf_names())) < len(outgroup):
            outgroup = list(child.get_leaf_names())
    return(outgroup)

# Walks through the gene tree.
# For each internal node labeled as a duplication event (D), it checks if the node only includes genes from one species.
# If yes, that means the duplication likely happened within that species (i.e., in-paralog), 
# so it changes the event type to "SD" (same-species duplication).
# It also explicitly records the species on the node.
def modify_in_paralog(tree):
    for node in tree.traverse():
        if hasattr(node, "evoltype") and node.evoltype == "D":
            if len(list(node.get_species())) == 1:
                node.evoltype = "SD"
                sd_species = list(node.get_species())[0]
                node.species = sd_species
    return(tree)


# This function processes the gene tree to define groups of orthologs and in-paralogs.
# First, it creates a copy of the tree (nhx_tree) and resets all "SD" nodes (same-species duplications) to "D" (duplication).
# Then, it splits the original gene tree into subtrees at duplication nodes.
# For each resulting subtree:
#   - If it contains no SD events, it's considered a clean ortholog group and is classified directly.
#   - If it contains SD events, these are handled using SD_graft(), which extracts and reclassifies the in-paralog groups.
# All final groups (orthologs or in-paralogs) are collected into final_paralog_list.
# Returns the list of classified groups and the cleaned nhx_tree for downstream annotation.
def fix_in_paralogs_define_groups(tree,focal_species):
    #print(tree)
    nhx_tree = copy.deepcopy(tree)
    for node in nhx_tree.traverse():
        if hasattr(node, "evoltype") and node.evoltype == "SD":
            node.evoltype = "D"
    final_paralog_list = []
    for subtree in tree.split_by_dups(autodetect_duplications=False):
        SD_events = count_SD_events(subtree)
        print(SD_events)
        print(subtree)
        if SD_events == 0:
            add_evolevent_feature(subtree,focal_species)
            final_paralog_list.append(subtree)
        else:
            in_paralogs = SD_graft(subtree)
            for paralog in in_paralogs:
                final_paralog_list.append(paralog)
    return(final_paralog_list , nhx_tree)

# Traverses a gene subtree and counts how many internal nodes are labeled as "SD" (same-species duplication).
# Used to determine whether a subtree contains in-paralog events that need special handling.
def count_SD_events(subtree):
    counter = 0
    for node in subtree.traverse():
        if hasattr(node, "evoltype") and node.evoltype == "SD":
            counter += 1
    return(counter)
        



#######################################################################################################################################
# Resolves subtrees containing in-paralogs (nodes labeled as "SD") to enforce single-copy representation per species.
#
# Behavior:
# - If the subtree contains only one species:
#     - All "SD" nodes are reclassified as true duplications ("D").
#     - The subtree is split into separate paralogs using .split_by_dups().
#
# - If multiple species are present:
#     - Iteratively identifies "SD" nodes where duplicated species are present.
#     - Selects one representative leaf (the farthest leaf) from the duplicated species.
#     - Makes a deep copy of this leaf and **grafts it back** onto the parent node.
#       ? This step ensures that **one representative is retained** and prevents over-pruning of that species.
#     - The rest of the SD subtree is detached and split as paralogs, with appropriate annotations.
#
# Annotations:
# - Subtrees are labeled with:
#     - "conserved" if all focal species are retained.
#     - "loss" with a list of missing species.
#     - "duplication" for each paralog group split from the original.
#
# Returns a list of cleaned and annotated subtrees, including one representative per focal species.
#####################################################################################################################################

def SD_graft(subtree):
    node_list = []
    if len(subtree.get_species()) == 1:
        for node in subtree.traverse():
            if hasattr(node, "evoltype") and node.evoltype == "SD":
                node.evoltype = "D"
        ##print("test")
        paralogs = subtree.detach()
        ##print(paralogs)
        for paralog in paralogs.split_by_dups(autodetect_duplications=False):
            ##print(9)
            ##print(paralog) 
            species = [x.species for x in paralog.get_leaves()]
            paralog.add_feature("evolevent","duplication")
            paralog.add_feature("duplicated_taxon",species)
            node_list.append(paralog)
    else:
        nnodespecies = len(subtree.get_species())
        while len([x.species for x in subtree.get_leaves()]) != nnodespecies:
            for node in subtree.traverse():
                ##print(1)
                #while hasattr(node, "evoltype") and node.evoltype == "SD" :
                if hasattr(node, "evoltype") and node.evoltype == "SD" :
                    nnodespecies = len(list(set([x.species for x in node.get_leaves()])))
                    farthest, test = node.get_farthest_leaf()
                    for leaf in node.get_leaves():
                        ##print(2)
                        dist = node.get_distance(leaf)
                        ##print(dist)
                        ##print(test)
                        if dist <= test:
                            dist = test
                            chosen_leaf = leaf
                            chosen_leaf_copy = copy.deepcopy(chosen_leaf)
                            chosen_leaf_copy.dist = dist
                        ## now we have the leaf to graft to node
                    print(3)
                    node.up.add_child(chosen_leaf_copy)
                    ### Here is our break.
                    paralogs = node.detach()
                    nspecies = [x.species for x in subtree.get_leaves()]
                    print(4)
                    if all([nspecies.count(x) == 1 for x in focal_species]):
                        print(5)
                        subtree.add_feature("evolevent","conserved")
                    elif any([nspecies.count(x) < 1 for x in focal_species]):
                        print(6)
                        lost_species_list = [x for x in focal_species if nspecies.count(x) < 1]
                        subtree.add_feature("evolevent","loss")
                        subtree.add_feature("lost_taxa",lost_species_list)
                    print(7)
                    break
            for node in paralogs.traverse():
                print(8)
                if hasattr(node, "evoltype") :
                    node.evoltype ="D"
            print(paralogs.evoltype)
            for paralog in paralogs.split_by_dups(autodetect_duplications=False):
                print(9)
                print(paralog.get_leaf_names())
                print(list(chosen_leaf_copy.get_leaf_names())[0])
                if list(paralog.get_leaf_names())[0] != list(chosen_leaf_copy.get_leaf_names())[0]:
                    species = [x.species for x in paralog.get_leaves()]
                    paralog.add_feature("evolevent","duplication")
                    paralog.add_feature("duplicated_taxon",species)
                    node_list.append(paralog)
                print(10)
            nnodespecies = len(subtree.get_species())
        node_list.append(subtree)
        print(subtree)
    return(node_list)


####################################################################################################################
# Assigns an evolutionary event label ("evolevent") to a node based on representation of focal species.
#
# - Labels as "conserved" if exactly one copy of each focal species is present.
# - Labels as "loss" if one or more focal species is missing, and records the missing taxa.
# - Labels internal "SD" nodes as "duplication" if any species appears more than once.
# - Prints a warning if the node doesn't match any expected pattern.
########################################################################################################################

def add_evolevent_feature(node,focal_species):
    #allspecies = list(tree.get_species())
    nspecies = [x.species for x in node.get_leaves()]
    if all([nspecies.count(x) == 1 for x in focal_species]):
        node.add_feature("evolevent","conserved")
    elif any([nspecies.count(x) < 1 for x in focal_species]):
        lost_species_list = [x for x in focal_species if nspecies.count(x) < 1]
        node.add_feature("evolevent","loss")
        node.add_feature("lost_taxa",lost_species_list)
    elif any([nspecies.count(x) > 1 for x in focal_species]):
        for child in node.traverse():
            if hasattr(child,"evolevent") and child.evolevent == "SD":
                child.add_feature("evolevent","duplication")
    else:
        print("uhhh")



# Filters a list of gene tree subtrees to retain only those containing focal species.
# Subtrees missing all focal species are removed; subtrees with partial focal coverage are pruned to retain only focal species tips.
def filter_nonfocal_species(subtree_list, focal_species):
    for subtree in subtree_list:
        species = list(subtree.get_species())
        if not any(x in species for x in focal_species):
            subtree_list.remove(subtree)
        elif not all(x in species for x in focal_species):
            keep_leaves = []
            for node in subtree:
                if hasattr(node, "species") and node.species in focal_species:
                    keep_leaves.append(node.name)
            subtree.prune(keep_leaves)
    return(subtree_list)
                    

# Builds a list of ortholog group outputs from gene tree subtrees.
# Each group is given a unique ID and its member gene names are concatenated into a space-separated string.
def build_orthologs_output(paralogs_list):
    orthologs_out = []
    counter = 0
    for paralog in paralogs_list:
        counter += 1
        gene_name = orthogroup + "-Gene-" + str(counter)
        string = ''
        for leaf_name in paralog.get_leaf_names():
            string = string + leaf_name
            if leaf_name != paralog.get_leaf_names()[-1]:
                string = string + ' '
        orthologs_out.append([gene_name,string])
    return(orthologs_out)

# Builds a key mapping each gene (leaf) to its assigned ortholog group.
# Used to track which genes belong to which orthogroup subcluster.
def build_ortholog_keys(paralogs_list):
    ortholog_keys = []
    counter = 0
    for paralog in paralogs_list:
        counter += 1
        gene_name = orthogroup + "-Gene-" + str(counter)
        for leaf_name in paralog.get_leaf_names():
            entry = [leaf_name,gene_name]
            ortholog_keys.append(entry)
    return(ortholog_keys)

# Builds a table summarizing the evolutionary event (e.g., conserved, loss, duplication)
# between every pairwise comparison of focal species for each ortholog group.
def build_classes_table(paralogs_list, focal_species):
    classes = []
    counter = 0
    species_comparisons = [comb for comb in combinations(focal_species, 2)]
    header = ['paralog']
    for comp in species_comparisons:
        comp_string = comp[0] + '-' + comp[1]
        header.append(comp_string)
    classes.append(header)    
    for paralog in paralogs_list:
        counter += 1
        gene_name = orthogroup + "-Gene-" + str(counter)
        entry = [gene_name]
        if paralog.evolevent == "conserved":
            n_conserved = ["conserved"] * len(species_comparisons)
            for i in n_conserved:
                entry.append(i)
        elif paralog.evolevent == "loss":
            for comp in species_comparisons:
                if comp[0] not in paralog.lost_taxa and comp[1] not in paralog.lost_taxa:
                    entry.append("conserved")
                elif comp[0] in paralog.lost_taxa and comp[1] in paralog.lost_taxa:
                    entry.append("NA")
                else:
                    entry.append("loss")
        elif paralog.evolevent == "duplication":
            for comp in species_comparisons:
                if comp[0] in paralog.duplicated_taxon or comp[1] in paralog.duplicated_taxon:
                    entry.append("duplication")
                else:
                    entry.append("NA")
        classes.append(entry)
    return(classes)

# Prints the unique leaf names in the NHX-format gene tree and species map for debugging name mismatches.
def fix_nhx_names(nhx_tree, species_map_df):
    print("[DEBUG] Unique NHX leaf names:")
    print(set(leaf.name for leaf in nhx_tree.iter_leaves()))
    print("[DEBUG] Map names:")
    print(set(species_map_df["New_name"]))

    # Add assertion here to catch unmatched species names early
    leaf_species_names = set(leaf.name.split("_")[0] for leaf in nhx_tree.iter_leaves())
    map_species_names = set(species_map_df["New_name"])

    missing = leaf_species_names - map_species_names
    assert leaf_species_names <= map_species_names, \
        f"[ASSERTION FAILED] NHX tree contains species not present in the mapping file: {missing}"

    # Now rename leaf nodes
    for leaf in nhx_tree.iter_leaves():
        try:
            species_name = leaf.name.split("_")[0]  # Only take the species part
            match = species_map_df[species_map_df["New_name"] == species_name]
            if match.empty:
                raise ValueError(f"[ERROR] No match in map for NHX leaf species: {species_name}")
            new_name = match["Sequence"].iloc[0]
            leaf.name = new_name
        except Exception as e:
            print(f"[ERROR] While fixing NHX name for {leaf.name}: {e}")
            raise

    return nhx_tree

# Renames leaf nodes in each paralog tree using the species map and filters out trees with unmapped species;
# ensures final trees contain only species present in the provided map and uses correct sequence names.
def fix_paralog_names(paralog_list, species_map_df):
    valid_species = set(species_map_df["New_name"])

    filtered_paralog_list = []
    skipped_species = set()

    for paralog in paralog_list:
        tree_species = {leaf.name.split("_")[0] for leaf in paralog.get_leaves()}
        if not tree_species <= valid_species:
            missing = tree_species - valid_species
            skipped_species.update(missing)
            continue  # Skip this paralog group entirely
        else:
            for leaf in paralog.get_leaves():
                try:
                    species_name = leaf.name.split("_")[0]
                    match = species_map_df[species_map_df["New_name"] == species_name]
                    if match.empty:
                        raise ValueError(f"[ERROR] No match in map for paralog leaf: {leaf.name}")
                    new_name = match["Sequence"].iloc[0]
                    leaf.name = new_name
                except Exception as e:
                    print(f"[ERROR] While fixing paralog name for {leaf.name}: {e}")
                    raise
            filtered_paralog_list.append(paralog)

    if skipped_species:
        print(f"[WARNING] Skipped {len(skipped_species)} species not in map: {skipped_species}")

    return filtered_paralog_list

# Creates a version of the NHX gene tree containing only focal species present in that tree;
# preserves topology and branch lengths by pruning based on actual leaf objects, not names.
def modify_nhx(nhx, focal_species):
    new_tree = copy.deepcopy(nhx)

    # Get all species present in this tree
    available_species = {leaf.species for leaf in new_tree.get_leaves()}

    # Dynamically filter the focal_species to only those found in this tree
    focal_present = available_species & set(focal_species)

    # Now keep leaf *objects* (not names) where species is in focal
    keep_leaves = [leaf for leaf in new_tree.get_leaves() if leaf.species in focal_present]

    # This avoids ambiguity and preserves branch lengths
    new_tree.prune(keep_leaves, preserve_branch_length=True)

    return new_tree


# Prunes the species tree to include only the focal species, used for consistency with pruned gene trees in downstream reconciliation or visualization.
def decostar_sp_tree(species_tree,focal_species):
    keep_leaves = [leaf.name for leaf in species_tree.get_leaves() if leaf.name in focal_species]
    edited_species_tree = copy.deepcopy(species_tree)
    edited_species_tree.prune(keep_leaves)
    return(edited_species_tree)

# Adds ancestral node labels to the species tree: assigns unique names to internal unlabeled nodes and full species names to terminal nodes for clarity in reconciliation.
def add_anc_names(st):
    possible_species = st.get_leaf_names()
    counter = 0
    for node in st.traverse():
        if node.species == '':
            sp_name = "anc" + str(counter)
            node.add_feature("S",sp_name)
            counter += 1
        else:
            putative_species = [x for x in possible_species if node.species in x]
            if len(putative_species) == 1:
                node.add_feature("S", putative_species[0])
    return(st)

# Converts a gene tree to Decostar-compatible NHX format by identifying the shared ancestor of focal species in the species tree and annotating the root accordingly.
def make_decostar_nhx(tree, focal_st):
    from ete3 import Tree

    decostar_tree = tree.copy()
    observed_species = set()

    # Collect species names from gene tree
    for leaf in decostar_tree.iter_leaves():
        try:
            observed_species.add(leaf.species)
        except:
            print(f"[WARN] Could not extract species from leaf {leaf.name}")

    # DEBUG: Print all species observed in this orthogroup
    print(f"[DEBUG] Observed species (raw): {observed_species}")

    # Convert observed species to match focal_st node names (remove underscores)
    formatted_species = [sp.replace("_", "") for sp in observed_species]

    # DEBUG: Print translated species names
    print(f"[DEBUG] Formatted species for common ancestor lookup: {formatted_species}")

    try:
        anc = focal_st.get_common_ancestor(formatted_species).S
    except Exception as e:
        print(f"[ERROR] Could not get common ancestor for: {formatted_species}")
        raise e

    # Annotate ancestor node in NHX format
    decostar_tree.add_features(S=anc)
    return decostar_tree

# BEGIN PATCH: write_summary function--- writes a summary file to output with the two species lists--- can be customized or skipped

def write_summary(orthogroup_name, orthologs_table, focal_species, cavefish_list, background_list, output_dir):
    summary_lines = []
    summary_lines.append(f"Orthogroup: {orthogroup_name}")
    summary_lines.append(f"Total species in input: {len(focal_species)}")

    # collect all unique species across all genes
    unique_species = set()
    for _, members in orthologs_table:
        for m in members.split():
            species = m.split("_")[0]  # Adjust if species parsing differs
            unique_species.add(species)
    summary_lines.append(f"Unique species in orthogroup: {len(unique_species)}")

    cavefish_present = [s for s in unique_species if s in cavefish_list]
    background_present = [s for s in unique_species if s in background_list]

    summary_lines.append(f"Cavefish species present: {len(cavefish_present)}")
    summary_lines.append(f"Background species present: {len(background_present)}\n")

    summary_lines.append("Breakdown per gene:")
    for gene, members in orthologs_table:
        these_species = set()
        for m in members.split():
            sp = m.split("_")[0]
            these_species.add(sp)
        n_cave = len([s for s in these_species if s in cavefish_list])
        n_bg = len([s for s in these_species if s in background_list])
        summary_lines.append(f"  {gene}: {len(these_species)} species ({n_cave} cavefish, {n_bg} background)")

    summary_file = os.path.join(output_dir, orthogroup_name, "summary.txt")
    with open(summary_file, "w") as f:
        f.write("\n".join(summary_lines))

# END PATCH
########################################
################# CODE #################
########################################
# Construct the full path to the species map file for the given orthogroup
map_file = os.path.join(map_dir, orthogroup + '.map')

# Load the species map file into a DataFrame with columns 'New_name' (species) and 'Sequence' (semicolon-delimited sequence IDs)
species_map_df = pd.read_csv(map_file, sep=":", header=None, names=["New_name", "Sequence"])

# Split multi-sequence entries into individual rows so each gene is mapped to a single species (1:1 format)
species_map_df = species_map_df.assign(Sequence=species_map_df["Sequence"].str.split(";")).explode("Sequence").reset_index(drop=True)

# Load the species tree from Newick format using ETE3 PhyloTree (used to guide reconciliation)
st = PhyloTree(species_tree)

# Dynamically determine which species are represented in this orthogroup’s map file
focal_species = list(species_map_df["New_name"].unique())

# Copy the full species tree so we can trim it down to only include species in this orthogroup
focal_st = copy.deepcopy(st)

# Prune the species tree to only contain leaves corresponding to the focal species for this orthogroup
focal_st.prune(focal_species)

# Add internal node names (ancestral species labels) to the pruned focal species tree
focal_st = add_anc_names(focal_st)

# Read in the reconciled gene tree for this orthogroup and check for consistency or duplication issues
tree = read_and_reconcile_tree(orthogroup, reconciled_dir)

# Run the core classification logic on the gene tree: define paralogs, identify duplications/losses, and clean up bad nodes
paralogs_list, nhx_tree = modify_classifications(tree, st, focal_species)

# Prune the NHX-format gene tree to include only focal species, preserving branch lengths for downstream tools
nwk_tree = modify_nhx(nhx_tree, focal_species)

# Convert the cleaned-up gene tree to Decostar-compatible NHX format using ancestral node naming from the focal species tree
decostar_nhx = make_decostar_nhx(nwk_tree, focal_st)

# Sanity check: verify the new tree is rooted properly before writing out
print(decostar_nhx.is_root())

# Fix leaf node names in the NHX tree to match the full sequence IDs from the species map (required for Decostar)
nhx_tree = fix_nhx_names(nwk_tree, species_map_df)

# Apply same name fix to all individual paralog subtrees, skipping any that don't match map species
paralogs_list = fix_paralog_names(paralogs_list, species_map_df)

# Build output table mapping orthogroup gene numbers to full sequence IDs
orthologs_table = build_orthologs_output(paralogs_list)

# Create a key file mapping every leaf gene ID to its assigned ortholog group name
ortholog_keys = build_ortholog_keys(paralogs_list)

# Generate the classification matrix summarizing conserved/loss/duplication status across species pairs
classes_table = build_classes_table(paralogs_list, focal_species)

# Prune the species tree again to include only species from the current orthogroup, for Decostar input
edited_species_tree = decostar_sp_tree(st, focal_species)

########################################
############### Outputs ################
########################################
########################################
############### Outputs ################
########################################

# Create the output directory for this orthogroup if it doesn't already exist
os.makedirs(os.path.join(output_dir, orthogroup), exist_ok=True)

# Write the full NHX gene tree with 'species' and 'evoltype' tags to file (used for tracing duplication/loss events)
out_tree = output_dir + "/" + orthogroup + "/" + orthogroup + ".nhx"
nhx_tree.write(features=["species", "evoltype"], outfile=out_tree, format_root_node=True)

# Write the pruned NHX tree in Newick format (format 9) for input to Decostar
out_nwk = output_dir + "/" + orthogroup + "/" + orthogroup + "_decostar.nwk"
nwk_tree.write(format=9, outfile=out_nwk, format_root_node=True)

# Write a version of the pruned Newick tree for codeml (phylogenetic model testing); same tree, just a different label
out_codeml = output_dir + "/" + orthogroup + "/" + orthogroup + "_codeml.nwk"
nwk_tree.write(outfile=out_codeml, format_root_node=True)

# Write Decostar NHX tree with internal node species assignments (S) and duplication flags (D)
out_deco = output_dir + "/" + orthogroup + "/" + orthogroup + "_decostar.nhx"
decostar_nhx.write(format=9, features=["S", "D"], outfile=out_deco, format_root_node=True)

# Write a simplified version of the Decostar tree for visualization or expression annotation (format 0 = plain NHX)
out_expr = output_dir + "/" + orthogroup + "/" + orthogroup + "_expr.nhx"
decostar_nhx.write(format=0, features=["S", "D"], outfile=out_expr, format_root_node=True)

# Write the pruned species tree (only focal species) in Newick format for Decostar or visual summaries
out_species = output_dir + "/" + orthogroup + "/" + orthogroup + "_pruned_species.nwk"
edited_species_tree.write(format=9, outfile=out_species, format_root_node=True)

# Write a CSV file listing each ortholog group with all gene names (space-separated)
orthologs_outfile = output_dir + '/' + orthogroup + "/" + orthogroup + '_orthogroup.csv'
with open(orthologs_outfile, 'w') as csv_file:
    csv_writer = csv.writer(csv_file, delimiter=',')
    for row in orthologs_table:
        csv_writer.writerow(row)
csv_file.close()

# Write a CSV key file linking each individual gene to its orthogroup ID
ortho_keys_outfile = output_dir + '/' + orthogroup + "/" + orthogroup + '_genes.csv'
with open(ortho_keys_outfile, 'w') as csv_file:
    csv_writer = csv.writer(csv_file, delimiter=',')
    csv_writer.writerow(['gene', 'orthogroup'])  # Header row
    for row in ortholog_keys:
        csv_writer.writerow(row)
csv_file.close()

# Write the classification matrix CSV (e.g., conserved/loss/duplication status for each species-pair comparison)
classes_table_outfile = output_dir + '/' + orthogroup + "/" + orthogroup + '_classes.csv'
with open(classes_table_outfile, 'w') as csv_file:
    csv_writer = csv.writer(csv_file, delimiter=',')
    for row in classes_table:
        csv_writer.writerow(row)
csv_file.close()

# Write a summary .txt file of orthogroup-level statistics for this gene family
write_summary(
    orthogroup_name=orthogroup,
    orthologs_table=orthologs_table,
    focal_species=focal_species,
    cavefish_list=cavefish_list,
    background_list=background_list,
    output_dir=output_dir
)
print(f"Summary file written to {os.path.join(output_dir, orthogroup, 'summary.txt')}")
