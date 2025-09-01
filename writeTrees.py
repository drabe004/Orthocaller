
#!/usr/bin/env python3

from Bio import Phylo
import os
import argparse
import copy
import subprocess
import re


def prune_tree(tree, species_list):
    for leaf in tree.get_terminals():
        if leaf.name not in species_list:
            tree.prune(leaf)

def main():
    parser = argparse.ArgumentParser(description="Generate pruned trees from focal species files.")
    parser.add_argument("master_tree", help="Path to the master tree file (in Newick format)")
    parser.add_argument("focal_dir", help="Directory containing *.focal.txt files")
    parser.add_argument("output_dir", help="Directory to write pruned trees")
    args = parser.parse_args()

    master_tree_file = args.master_tree
    focal_dir = args.focal_dir
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)

    master_tree = Phylo.read(master_tree_file, "newick")

    for file in os.listdir(focal_dir):
        if not file.endswith(".focal.txt"):
            continue

        focal_path = os.path.join(focal_dir, file)
        with open(focal_path) as f:
            species_list = [re.sub(r'_\d+$', '', line.strip()) for line in f if line.strip()]

        pruned_tree = copy.deepcopy(master_tree)
        prune_tree(pruned_tree, species_list)

        output_path = os.path.join(output_dir, file.replace(".focal.txt", ".tre"))
        Phylo.write(pruned_tree, output_path, "newick")

        subprocess.run(["sed", "-i", "-E", "s/:[0-9]+\.?[0-9]*//g", output_path])

if __name__ == "__main__":
    main()
