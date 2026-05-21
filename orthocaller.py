#!/usr/bin/env python

# Additional software necessary to run this:

# (1) biopython
# (2) pandas
# (3) dfply
# (4) mafft
# (5) iqtree
# (6) generax v.XXX

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
import glob ##Need this for wildcard handling

#from ete3 import Phyloxml, phyloxml

from ctypes import alignment
import sys, os, shutil
import datetime as dt
#import numpy as np
from dfply import *


########################################
############### ARGUMENTS ##############
########################################

parser = argparse.ArgumentParser(description='')
parser.add_argument("-og","--orthogroup",
					type=str,
					help="Orthogroup to process.")
parser.add_argument("-st","--species_tree",
					type=str,
                    default="0_Inputs/species_tree.nwk",
					help="Species tree file for orthofinder. Newick format")
parser.add_argument("-rd","--reconciled_dir",
					type=str,
                    default="3.5_OG_reconciliations",
					help="Directory containing generax input and output files")
parser.add_argument("-md","--map_dir",
					type=str,
                    default="3.1_Species_maps",
					help="Directory containing sequence")
parser.add_argument("-fs", "--focal_species_file",
                    type=str,
                    required=True,
                    help="Path to a file containing a list of focal species, one per line.")
parser.add_argument("-od","--output_dir",
					type=str,
                    default="3.6_OG_classifications",
					help="Directory of where to write a folder for this orthogroup and subsequent output files")
parser.add_argument(
    "--inparalog_strategy",
    choices=["farthest", "shortest", "average_divergence"],
    default="farthest",
    help=(
        "Strategy for choosing which in-paralog to retain:\n"
        "  farthest - the most diverged copy (default)\n"
        "  shortest - the least diverged copy\n"
        "  average_divergence - the copy whose branch length ratio is closest to the duplication node's branch"
    )
)
parser.add_argument(
    "--cavefish_list",
    type=str,
    default="Cavefish_List.txt",
    help="Text file with newline-separated cavefish species names."
)

parser.add_argument(
    "--background_list",
    type=str,
    default="BackgroundFish_List.txt",
    help="Text file with newline-separated background (non-cavefish) species names."
)

parser.add_argument(
    "--long_branch_threshold",
    type=float,
    default=100000,
    help="Threshold for long branches"
)

parser.add_argument(
    "--short_branch_threshold",
    type=float,
    default=0,
    help="Threshold for short branches"
)

args = parser.parse_args()

########################################
################# SETUP ################
########################################
# Load cavefish and background species lists
with open(args.cavefish_list) as f:
    cavefish_list = [line.strip() for line in f if line.strip()]

with open(args.background_list) as f:
    background_list = [line.strip() for line in f if line.strip()]




orthogroup = args.orthogroup
reconciled_dir = args.reconciled_dir
map_dir = args.map_dir
species_tree = args.species_tree
focal_file = args.focal_species_file
output_dir = args.output_dir
LONG_BRANCH_THRESHOLD = args.long_branch_threshold
SHORT_BRANCH_THRESHOLD = args.short_branch_threshold

########################################
############## FUNCTIONS ###############
########################################

def find_species(focal_file):
    """
    Read a newline-delimited focal species file.

    Parameters
    ----------
    focal_file : str
        Path to a text file containing one species name per line.

    Returns
    -------
    list
        Non-empty species names stripped of leading and trailing whitespace.
    """
    with open(focal_file, 'r') as f:
        return [line.strip() for line in f if line.strip()]



def parse_sp_name(node_name):
    """
    Extract the species identifier from a renamed sequence header.

    Parameters
    ----------
    node_name : str
        Sequence name formatted with the species identifier before the first underscore.

    Returns
    -------
    str
        Species identifier parsed from the sequence name.
    """
    return node_name.split("_")[0]


def read_and_reconcile_tree(orthogroup, reconciled_dir, map):
    """
    Load the GeneRax-reconciled gene tree and reconcile it against the species tree.

    Parameters
    ----------
    orthogroup : str
        Orthogroup identifier being processed.
    reconciled_dir : str
        Directory containing GeneRax reconciliation output.
    map : pandas.DataFrame
        Mapping table with original sequence names and renamed species-aware identifiers.

    Returns
    -------
    ete3.PhyloTree
        Reconciled gene tree with renamed leaves and species parsing enabled.
    """
    tree_string = glob.glob(reconciled_dir + "/" + orthogroup + "/reconciliations/*.newick")[0]
    ##Danielle had to change this because Phylotree cannot handle wildcards like *.newick
    #tree_string = reconciled_dir + "/" + orthogroup + "/" + "reconciliations/*.newick"
    tree = PhyloTree(tree_string, format=1)
    ## Renames tree leaves to have Species
    for leaf in tree.get_tree_root():
        new_name = map[map["Sequence"] == leaf.name]["New_name"].iloc[0]
        leaf.name = new_name
    tree.set_species_naming_function(parse_sp_name)
    #genetree = PhyloTree(tree)
    #st = PhyloTree(st)
    recon_tree, events = tree.reconcile(st)
    return(tree)

def reclassify_postprune_duplications(tree):
    """
    Reclassify duplication nodes that no longer contain duplicated species.

    Parameters
    ----------
    tree : ete3.Tree
        Reconciled gene tree containing evoltype annotations.

    Returns
    -------
    None
        The input tree is modified in place.
    """
    for node in tree.traverse():
        if hasattr(node, "evoltype") and node.evoltype == "D":
            species = [x.species for x in node.get_leaves()]
            if all(species.count(sp) == 1 for sp in species):
                print(f"[INFO] Reclassifying node {node.name} from D ? S (unique species)")
                node.evoltype = "S"



def modify_classifications(tree, st, focal_species):
    """
    Run the Orthocaller classification workflow on a reconciled gene tree.

    This applies post-pruning duplication reclassification, in-paralog labeling,
    topology correction, putative false-orthology detection, paralog splitting,
    and non-focal species filtering.

    Parameters
    ----------
    tree : ete3.Tree
        Reconciled gene tree.
    st : ete3.Tree
        Species tree used for reconciliation and topology correction.
    focal_species : list
        Species retained for downstream orthogroup classification.

    Returns
    -------
    tuple
        Final paralog subtree list and an NHX tree retaining duplication annotations.
    """
    reclassify_postprune_duplications(tree)
    tree = modify_in_paralog(tree)
    identify_and_resolve_topology_problems(tree,st)
    identify_putative_false_orthology(tree)
    final_paralogs_list,nhx_tree = fix_in_paralogs_define_groups(tree, focal_species)
    final_paralogs_list = filter_nonfocal_species(final_paralogs_list, focal_species)
    return(final_paralogs_list, nhx_tree)


def identify_putative_false_orthology(tree):
    """
    Identify putative false orthology caused by long branch patterns.

    Parameters
    ----------
    tree : ete3.Tree
        Reconciled gene tree containing speciation and duplication annotations.

    Returns
    -------
    None
        Nodes may be reclassified in place when long-branch criteria are met.
    """
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


def long_br(descendent,ref):
    """
    Test whether a descendant branch is long relative to a reference branch.

    Parameters
    ----------
    descendent : ete3.TreeNode
        Descendant node whose branch length is evaluated.
    ref : float
        Reference branch length used as the denominator for the ratio.

    Returns
    -------
    bool
        True if the descendant/reference ratio exceeds LONG_BRANCH_THRESHOLD.
    """
    if ref == 0.0:
        ref = 0.000001
    test = descendent.dist
    if test == 0.0:
        test = 0.000001
    print("test", test)
    val = test /ref
    print("val", val)
    if val > LONG_BRANCH_THRESHOLD and ref != 0.000001:
        return True
    else:
        return False


def identify_and_resolve_topology_problems(tree,st):
    """
    Detect and repair topology issues near duplication/speciation boundaries.

    Parameters
    ----------
    tree : ete3.Tree
        Reconciled gene tree to inspect and modify.
    st : ete3.Tree
        Species tree used to re-root and re-reconcile problematic subtrees.

    Returns
    -------
    None
        Problematic nodes are modified in place when criteria are met.
    """
    for node in tree.traverse():
        if hasattr(node, 'evoltype') and node.evoltype == 'D' :
            for child in node.get_children():
                if hasattr(child,'evoltype') and child.evoltype == 'S' :
                    #print("running test")
                    if short_br(node,child) and species_rep_test(node):
                        #node_to_fix = node.detach()
                        #print(node)
                        fix_node(tree,node,st)


def short_br(node,child):
    """
    Test whether a child branch is short relative to its parent duplication branch.

    Parameters
    ----------
    node : ete3.TreeNode
        Parent node used as the reference branch.
    child : ete3.TreeNode
        Child node whose distance from the parent is evaluated.

    Returns
    -------
    bool
        True if the child/reference ratio is below SHORT_BRANCH_THRESHOLD.
    """
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
    if val < SHORT_BRANCH_THRESHOLD:
        return True
    else:
        return False


## Tests number of species representatives below a given node accounting for in paralogs.
## Passes (True) if only 1 per species. Fails (False) if more than 1 after account for in-paralogs.
def species_rep_test(node):
    """
    Test whether a node contains at most one representative per species.

    Species-specific duplication events are accounted for before evaluating
    whether any species appears more than once below the node.

    Parameters
    ----------
    node : ete3.TreeNode
        Tree node whose descendant species representation is evaluated.

    Returns
    -------
    bool
        True if species representation is compatible with one copy per species.
    """
    species = [x.species for x in node.get_leaves()]
    sd_species = [x.species for x in node.traverse() if hasattr(node,"evoltype") and node.evoltype == "SD"]
    if len(sd_species) != 0:
        for i in sd_species:
            species.remove(i)
    if any([species.count(x) > 1 for x in species]):
        return False
    else:
        return True


def fix_node(tree, bad_node, st):
    """
    Re-root and re-reconcile a problematic subtree against a pruned species tree.

    Parameters
    ----------
    tree : ete3.Tree
        Full gene tree containing the problematic node.
    bad_node : ete3.TreeNode
        Node to replace after topology correction.
    st : ete3.Tree
        Species tree used to select an outgroup and reconcile the corrected subtree.

    Returns
    -------
    None
        The corrected subtree is grafted back into the original tree in place.
    """
    node = copy.deepcopy(bad_node)
    spare_st = copy.deepcopy(st)
    represented_species = [x for x in spare_st.get_leaf_names() if x in node.get_species()]
    spare_st.prune(represented_species)
    outgroup_species = find_reference_outgroup(spare_st)
    species_instances = [x.species for x in node.get_leaves()]

    if len(outgroup_species) > 1 or species_instances.count(outgroup_species[0]) > 1:
        species_leaves = search_for_species(node, outgroup_species)
        if not species_leaves:
            print(f"[WARNING] No species leaves found for outgroup {outgroup_species} - skipping re-rooting.")
        else:
            out_anc = node.get_common_ancestor(species_leaves)
            if node == out_anc:
                print(f"[INFO] Skipping re-rooting: node {node.name} is its own outgroup (symmetric or shallow).")
            else:
                node.set_outgroup(out_anc)
    else:
        species_leaves = search_for_species(node, outgroup_species)
        if not species_leaves:
            print(f"[WARNING] No species leaves found for fallback outgroup {outgroup_species} - skipping.")
        else:
            out_leaf = species_leaves[0]
            if node == out_leaf:
                print(f"[INFO] Skipping re-rooting: node {node.name} is its own fallback outgroup leaf.")
            else:
                node.set_outgroup(out_leaf)

    rec_node, node_events = node.reconcile(spare_st)
    bad_node.up.add_child(node)
    bad_node.detach()
    

def search_for_species(tree, species):
    """
    Find all nodes in a tree matching one or more species names.

    Parameters
    ----------
    tree : ete3.Tree
        Tree to search.
    species : list
        Species names to locate.

    Returns
    -------
    list
        Nodes whose species attribute matches one of the requested species.
    """
    node_list = []
    for entry in species:
        for node in tree.traverse():
            if hasattr(node,"species") and node.species == entry:
                node_list.append(node)
    return(node_list)


def find_reference_outgroup(st):
    """
    Select a reference outgroup from a species tree.

    Parameters
    ----------
    st : ete3.Tree
        Species tree pruned to the species represented in a subtree.

    Returns
    -------
    list
        Species name or species names used as the reference outgroup.
    """
    outgroup = list(st.get_leaf_names())
    for child in st.get_children():
        if child.is_leaf() == True:
            outgroup = [child.name]
            break
        elif len(list(child.get_leaf_names())) < len(outgroup):
            outgroup = list(child.get_leaf_names())
    return(outgroup)


def modify_in_paralog(tree):
    """
    Label species-specific duplication nodes as SD events.

    Parameters
    ----------
    tree : ete3.Tree
        Reconciled gene tree containing duplication annotations.

    Returns
    -------
    ete3.Tree
        Tree with single-species duplication nodes relabeled as SD.
    """
    for node in tree.traverse():
        if hasattr(node, "evoltype") and node.evoltype == "D":
            if len(list(node.get_species())) == 1:
                node.evoltype = "SD"
                sd_species = list(node.get_species())[0]
                node.species = sd_species
    return(tree)



def fix_in_paralogs_define_groups(tree,focal_species):
    """
    Split duplication-defined subtrees and resolve in-paralog groups.

    Parameters
    ----------
    tree : ete3.Tree
        Reconciled gene tree with SD events labeled.
    focal_species : list
        Species used to classify conservation, loss, or duplication.

    Returns
    -------
    tuple
        Final paralog subtree list and a copy of the NHX tree with SD nodes restored to D.
    """
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


def count_SD_events(subtree):
    """
    Count species-specific duplication events within a subtree.

    Parameters
    ----------
    subtree : ete3.Tree
        Subtree to scan for nodes labeled with evoltype == "SD".

    Returns
    -------
    int
        Number of SD events detected.
    """
    counter = 0
    for node in subtree.traverse():
        if hasattr(node, "evoltype") and node.evoltype == "SD":
            counter += 1
    return(counter)
        

#  PATCHED SD_graft with strategy block inserted
def SD_graft(subtree):
    """
    Resolve species-specific duplication events by retaining one in-paralog representative.

    For each SD node, this function selects one retained copy according to
    args.inparalog_strategy and separates remaining copies into duplication-derived
    paralog groups.

    Parameters
    ----------
    subtree : ete3.Tree
        Subtree containing one or more SD events.

    Returns
    -------
    list
        Resolved ortholog and paralog subtrees after SD processing.
    """
    node_list = []
    if len(subtree.get_species()) == 1:
        for node in subtree.traverse():
            if hasattr(node, "evoltype") and node.evoltype == "SD":
                node.evoltype = "D"
        paralogs = subtree.detach()
        for paralog in paralogs.split_by_dups(autodetect_duplications=False):
            species = [x.species for x in paralog.get_leaves()]
            paralog.add_feature("evolevent","duplication")
            paralog.add_feature("duplicated_taxon",species)
            node_list.append(paralog)
    else:
        nnodespecies = len(subtree.get_species())
        while len([x.species for x in subtree.get_leaves()]) != nnodespecies:
            for node in subtree.traverse():
                if hasattr(node, "evoltype") and node.evoltype == "SD":
                    nnodespecies = len(list(set([x.species for x in node.get_leaves()])))

                    # === In-paralog Options Code Patch ===
                    leaves = node.get_leaves()
                    distances = {leaf: node.get_distance(leaf) for leaf in leaves}

                    if args.inparalog_strategy == "farthest":
                        chosen_leaf = max(distances, key=distances.get)

                    elif args.inparalog_strategy == "shortest":
                        chosen_leaf = min(distances, key=distances.get)

                    elif args.inparalog_strategy == "average_divergence":
                        ref = node.dist
                        if ref == 0.0:
                            ref = 0.000001
                        chosen_leaf = min(
                            leaves,
                            key=lambda leaf: abs((node.get_distance(leaf) / ref) - 1.0)
                        )

                    chosen_leaf_copy = copy.deepcopy(chosen_leaf)
                    chosen_leaf_copy.dist = distances[chosen_leaf]
                    node.up.add_child(chosen_leaf_copy)
                    print(f"[INFO] Retained {chosen_leaf.name} using strategy: {args.inparalog_strategy}")
                    # === End Patch ===

                    print(3)
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
                if hasattr(node, "evoltype"):
                    node.evoltype = "D"
            print(paralogs.evoltype)
            for paralog in paralogs.split_by_dups(autodetect_duplications=False):
                print(9)
                print(paralog.get_leaf_names())
                print(list(chosen_leaf_copy.get_leaf_names())[0])
                if list(paralog.get_leaf_names())[0] != list(chosen_leaf_copy.get_leaf_names())[0]:
                    species = [x.species for x in paralog.get_leaves()]
                    paralog.add_feature("evolevent", "duplication")
                    paralog.add_feature("duplicated_taxon", species)
                    node_list.append(paralog)
                print(10)
            nnodespecies = len(subtree.get_species())
        node_list.append(subtree)
        print(subtree)
    return(node_list)


def add_evolevent_feature(node,focal_species):
    """
    Assign conserved, loss, or duplication classification to a subtree.

    Parameters
    ----------
    node : ete3.TreeNode
        Subtree root to classify.
    focal_species : list
        Species expected in the orthogroup.

    Returns
    -------
    None
        Classification features are added directly to the node.
    """
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




def filter_nonfocal_species(subtree_list, focal_species):
    """
    Remove or prune subtrees that lack complete focal-species representation.

    Parameters
    ----------
    subtree_list : list
        Candidate ortholog/paralog subtrees.
    focal_species : list
        Species to retain for downstream analysis.

    Returns
    -------
    list
        Filtered subtree list containing focal species only.
    """
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
                    

def build_orthologs_output(paralogs_list):
    """
    Build the orthogroup membership output table.

    Parameters
    ----------
    paralogs_list : list
        Final ortholog/paralog subtrees.

    Returns
    -------
    list
        Rows containing generated gene names and space-delimited member sequence names.
    """
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


def build_ortholog_keys(paralogs_list):
    """
    Build a sequence-to-orthogroup key table.

    Parameters
    ----------
    paralogs_list : list
        Final ortholog/paralog subtrees.

    Returns
    -------
    list
        Rows mapping each sequence name to its generated Orthocaller gene name.
    """
    ortholog_keys = []
    counter = 0
    for paralog in paralogs_list:
        counter += 1
        gene_name = orthogroup + "-Gene-" + str(counter)
        for leaf_name in paralog.get_leaf_names():
            entry = [leaf_name,gene_name]
            ortholog_keys.append(entry)
    return(ortholog_keys)


def build_classes_table(paralogs_list, focal_species):
    """
    Build pairwise focal-species classification table for each paralog group.

    Parameters
    ----------
    paralogs_list : list
        Final classified ortholog/paralog subtrees.
    focal_species : list
        Species used to construct pairwise comparison columns.

    Returns
    -------
    list
        Classification table with conserved, loss, duplication, or NA entries.
    """
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


def fix_nhx_names(nhx_tree,map):
    """
    Convert NHX tree leaf names back to original sequence identifiers.

    Parameters
    ----------
    nhx_tree : ete3.Tree
        Annotated NHX tree with renamed leaves.
    map : pandas.DataFrame
        Mapping table linking renamed leaves to original sequence identifiers.

    Returns
    -------
    ete3.Tree
        NHX tree with original sequence names restored.
    """
    for leaf in nhx_tree.get_leaves():
        new_name = map[map["New_name"] == leaf.name]["Sequence"].iloc[0]
        leaf.name = new_name
    return(nhx_tree)


def fix_paralog_names(paralog_list, map):
    """
    Convert paralog subtree leaf names back to original sequence identifiers.

    Parameters
    ----------
    paralog_list : list
        Final paralog subtrees with renamed leaves.
    map : pandas.DataFrame
        Mapping table linking renamed leaves to original sequence identifiers.

    Returns
    -------
    list
        Paralog subtrees with original sequence names restored.
    """
    for paralog in paralog_list:
        for leaf in paralog.get_leaves():
            new_name = map[map["New_name"] == leaf.name]["Sequence"].iloc[0]
            leaf.name = new_name
    return(paralog_list)


def modify_nhx(nhx,focal_species):
    """
    Prune an NHX tree to focal species.

    Parameters
    ----------
    nhx : ete3.Tree
        Annotated tree to copy and prune.
    focal_species : list
        Species retained in the output tree.

    Returns
    -------
    ete3.Tree
        Pruned copy of the input tree.
    """
    new_tree = copy.deepcopy(nhx)
    keep_leaves = [leaf.name for leaf in new_tree.get_leaves() if leaf.species in focal_species]
    new_tree.prune(keep_leaves)
    return(new_tree)


def decostar_sp_tree(species_tree,focal_species):
    """
    Prune the species tree for DeCoSTAR output.

    Parameters
    ----------
    species_tree : ete3.Tree
        Full species tree.
    focal_species : list
        Species to retain.

    Returns
    -------
    ete3.Tree
        Copy of the species tree pruned to focal species.
    """
    keep_leaves = [leaf.name for leaf in species_tree.get_leaves() if leaf.name in focal_species]
    edited_species_tree = copy.deepcopy(species_tree)
    edited_species_tree.prune(keep_leaves)
    return(edited_species_tree)


def add_anc_names(st):
    """
    Add ancestral species labels required for downstream NHX output.

    Parameters
    ----------
    st : ete3.Tree
        Species tree to annotate.

    Returns
    -------
    ete3.Tree
        Species tree with S features added to ancestral and terminal nodes.
    """
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



#def make_anc_dict(focal_st):
#    anc_dict = {}
#    for node in focal_st.traverse():
#        species = tuple(node.get_leaf_names())
#        anc_dict[species] = node.S
#    return(anc_dict)



def make_decostar_nhx(nwk_tree, focal_st):
    """
    Create a DeCoSTAR-compatible NHX tree with duplication and species annotations.

    Parameters
    ----------
    nwk_tree : ete3.Tree
        Pruned gene tree containing evoltype annotations.
    focal_st : ete3.Tree
        Focal species tree annotated with ancestral species names.

    Returns
    -------
    ete3.Tree
        NHX tree annotated with D, S, and geneName features.
    """
    decostar_nhx = copy.deepcopy(nwk_tree)
    for node in decostar_nhx.traverse():
        if hasattr(node, 'evoltype') and node.evoltype == "D":
            node.add_feature("D",'Y')
        else:
            node.add_feature("D","N")
    for node in decostar_nhx.traverse():
        all_observed_species = [x.split('_')[0] for x in node.get_leaf_names()]
        uniq_observed_species = list(set(all_observed_species))
        if node.is_leaf():
            anc = node.get_leaf_names()[0].split('_')[0]
            node.add_feature("S",anc)
            node.add_feature("geneName",node.get_leaf_names()[0])
        elif node.evoltype == "S":
            anc = focal_st.get_common_ancestor(uniq_observed_species).S
            node.add_feature("S",anc)
        elif node.evoltype == "D" and len(uniq_observed_species) == 1:
            anc = uniq_observed_species[0]
            node.add_feature("S",anc)
        elif node.evoltype == "D" :
            anc = focal_st.get_common_ancestor(uniq_observed_species).S
            node.add_feature("S",anc)
    return(decostar_nhx)
# BEGIN PATCH: write_summary function--- writes a summary file to output with the two species lists--- can be customized or skipped
  
def write_summary(orthogroup_name, orthologs_table, focal_species, cavefish_list, background_list, output_dir):
    """
    Write an orthogroup-level summary file.

    The summary includes total focal species, unique species observed,
    cavefish/background counts, and per-gene species breakdowns.

    Parameters
    ----------
    orthogroup_name : str
        Orthogroup identifier.
    orthologs_table : list
        Ortholog membership table produced by build_orthologs_output.
    focal_species : list
        Focal species used in the analysis.
    cavefish_list : list
        Species classified as cavefish.
    background_list : list
        Species classified as background fish.
    output_dir : str
        Base output directory for Orthocaller results.

    Returns
    -------
    None
        Writes summary.txt into the orthogroup output directory.
    """
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

st = PhyloTree(species_tree)

focal_species = find_species(focal_file)
if '' in focal_species:
    focal_species.remove('')
if '.dir' in focal_species:
    focal_species.remove('.dir')


focal_st = copy.deepcopy(st)
focal_st.prune(focal_species)
focal_st = add_anc_names(focal_st)

map_file = map_dir + "/" + orthogroup + '.map'
map = pd.read_csv(map_file)


tree = read_and_reconcile_tree(orthogroup, reconciled_dir, map)

paralogs_list,nhx_tree = modify_classifications(tree, st, focal_species)

nwk_tree = modify_nhx(nhx_tree,focal_species)

decostar_nhx = make_decostar_nhx(nwk_tree, focal_st)
print(decostar_nhx.is_root())


nhx_tree = fix_nhx_names(nhx_tree, map)

paralogs_list = fix_paralog_names(paralogs_list, map)

orthologs_table = build_orthologs_output(paralogs_list)

ortholog_keys = build_ortholog_keys(paralogs_list)

classes_table = build_classes_table(paralogs_list, focal_species)

edited_species_tree = decostar_sp_tree(st,focal_species)

########################################
############### Outputs ################
########################################

#mk_output = "mkdir " + output_dir + "/" + orthogroup
#sp.call(mk_output, shell=True)
##This threw dir errors when used in loop
os.makedirs(os.path.join(output_dir, orthogroup), exist_ok=True)



out_tree = output_dir + "/" + orthogroup + "/" + orthogroup + ".nhx"
nhx_tree.write(features = ["species","evoltype"], outfile=out_tree, format_root_node=True)

out_nwk = output_dir + "/" + orthogroup + "/" + orthogroup + "_decostar.nwk"
nwk_tree.write(format=9, outfile=out_nwk, format_root_node=True)
out_codeml = output_dir + "/" + orthogroup + "/" + orthogroup + "_codeml.nwk"
nwk_tree.write(outfile=out_codeml, format_root_node=True)


out_deco = output_dir + "/" + orthogroup + "/" + orthogroup + "_decostar.nhx"
decostar_nhx.write(format=9,features = ["S","D"], outfile=out_deco, format_root_node=True)

out_expr = output_dir + "/" + orthogroup + "/" + orthogroup + "_expr.nhx"
decostar_nhx.write(format=0,features = ["S","D"], outfile=out_expr, format_root_node=True)


out_species = output_dir + "/" + orthogroup + "/" + orthogroup + "_pruned_species.nwk"
edited_species_tree.write(format=9,outfile=out_species, format_root_node=True)


orthologs_outfile = output_dir + '/' + orthogroup + "/" + orthogroup + '_orthogroup.csv'
with open(orthologs_outfile, 'w') as csv_file:
    csv_writer = csv.writer(csv_file, delimiter = ',')
    for row in orthologs_table:
        csv_writer.writerow(row)		
csv_file.close()


ortho_keys_outfile = output_dir + '/' + orthogroup + "/" + orthogroup + '_genes.csv'
with open(ortho_keys_outfile, 'w') as csv_file:
    csv_writer = csv.writer(csv_file, delimiter = ',')
    csv_writer.writerow(['gene','orthogroup'])
    for row in ortholog_keys:
        csv_writer.writerow(row)		
csv_file.close()


classes_table_outfile = output_dir + '/' + orthogroup + "/" + orthogroup + '_classes.csv'
with open(classes_table_outfile, 'w') as csv_file:
    csv_writer = csv.writer(csv_file, delimiter = ',')
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



#tree_string = reconciled_dir + "/" + orthogroup + "/" + orthogroup + "_output/results/family_1/geneTree.newick"
#    tree = PhyloTree(tree_string)
