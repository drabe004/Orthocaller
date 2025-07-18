#!/usr/bin/env python3

import os
import re
import argparse

def parse_args():
    parser = argparse.ArgumentParser(
        description="Summarize gene lines with cavefish and background species thresholds."
    )
    parser.add_argument(
        "-i", "--input_dir",
        required=True,
        help="Base input directory containing <n>_generax/summary.txt folders"
    )
    parser.add_argument(
        "-o", "--output_file",
        required=True,
        help="Path to write master summary file"
    )
    parser.add_argument(
        "-n", "--num_orthogroups",
        type=int,
        default=12171,
        help="Number of orthogroups to search (default: 12171)"
    )
    parser.add_argument(
        "--min-cavefish",
        type=int,
        default=8,
        help="Minimum number of cavefish species (default: 8)"
    )
    parser.add_argument(
        "--min-background",
        type=int,
        default=31,
        help="Minimum number of background species (default: 31)"
    )
    return parser.parse_args()

def main():
    args = parse_args()
    input_dir = args.input_dir
    output_path = args.output_file
    num_ogs = args.num_orthogroups
    min_cavefish = args.min_cavefish
    min_background = args.min_background

    matched_lines = 0
    opened_files = 0

    with open(output_path, "w") as out:
        for i in range(1, num_ogs + 1):
            folder_name = f"{i}_generax"
            summary_path = os.path.join(input_dir, folder_name, "summary.txt")

            if not os.path.exists(summary_path):
                continue
            opened_files += 1

            with open(summary_path, "r") as f:
                lines = f.readlines()
                for line in lines:
                    line = line.replace('\r', '').strip()  # Handle Windows-style carriage returns

                    if "cavefish" in line and "Gene" in line:
                        # Match gene lines and extract cavefish/background counts
                        match = re.search(r"Gene-\d+: \d+ species \((\d+) cavefish, (\d+) background\)", line)
                        if match:
                            cavefish = int(match.group(1))
                            background = int(match.group(2))
                            if cavefish >= min_cavefish and background >= min_background:
                                out.write(line + "\n")
                                matched_lines += 1

    print(f"? Processed {opened_files} summary files.")
    print(f"? Found {matched_lines} matching gene lines.")
    print(f"? Output written to: {output_path}")

if __name__ == "__main__":
    main()
