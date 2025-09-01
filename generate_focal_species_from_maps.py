#!/usr/bin/env python3

import os
import argparse
import pandas as pd
import glob
import re

def extract_focal_species(map_file, output_file):
    try:
        df = pd.read_csv(map_file)
        if "New_name" not in df.columns:
            raise ValueError("Missing 'New_name' column in map file.")

        # Strip any _number suffix from New_name to match species tree
        df["Clean_name"] = df["New_name"].apply(lambda x: re.sub(r"_\d+$", "", x))
        focal_species = sorted(df["Clean_name"].drop_duplicates())

        with open(output_file, "w") as out:
            for sp in focal_species:
                out.write(sp + "\n")

        print(f"? Wrote {len(focal_species)} species to {output_file}")
    except Exception as e:
        print(f"? Failed to process {map_file}: {e}")

def main():
    parser = argparse.ArgumentParser(description="Generate focal species files from map files.")
    parser.add_argument("-m", "--map_dir", required=True, help="Directory containing .map files")
    parser.add_argument("-o", "--output_dir", required=True, help="Directory to write .focal.txt files")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    map_files = glob.glob(os.path.join(args.map_dir, "*.map"))

    if not map_files:
        print("No .map files found in the specified directory.")
        return

    for map_file in map_files:
        base = os.path.splitext(os.path.basename(map_file))[0]
        output_file = os.path.join(args.output_dir, f"{base}.focal.txt")
        extract_focal_species(map_file, output_file)

if __name__ == "__main__":
    main()
