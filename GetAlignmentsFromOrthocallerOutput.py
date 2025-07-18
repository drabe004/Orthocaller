# Debugging: Add more informative print statements
# and enhance record ID matching logic for robustness

import os
import re
from pathlib import Path
from Bio import SeqIO
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Extract gene-specific alignments from orthogroup data using GeneRaxKey.")
    parser.add_argument("-s", "--summary_file", required=True, help="Path to master summary file.")
    parser.add_argument("-o", "--output_dir", required=True, help="Output directory to write alignments.")
    parser.add_argument("-b", "--base_dir", required=True, help="Base directory for orthogroup classifications.")
    parser.add_argument("-a", "--alignment_dir", required=True, help="Directory containing all alignment FASTA files.")
    parser.add_argument("-k", "--key_file", required=True, help="GeneRaxKey file mapping alignment file to orthogroup.")
    parser.add_argument("--dry-run", action="store_true", help="Only print what would be done, do not write any output files.")
    parser.add_argument("--limit", type=int, default=None, help="Only process the first N genes from the summary file.")
    return parser.parse_args()

def read_generax_key(key_path):
    mapping = {}
    with open(key_path) as f:
        for line in f:
            aln, og = line.strip().split(",")
            mapping[og.strip()] = aln.strip()
    return mapping

def extract_sequences_to_fasta(seq_names, fasta_path, out_path, dry_run=False):
    matched_records = []

    # Flatten space-separated entries if they were read as a single string
    if len(seq_names) == 1 and ' ' in seq_names[0]:
        seq_names = seq_names[0].split()

    target_ids = set(name.strip() for name in seq_names)
    print(f"  [DEBUG] Target sequence IDs: {list(target_ids)[:5]} ... (total {len(target_ids)})")

    for record in SeqIO.parse(fasta_path, "fasta"):
        # Try exact match
        if record.id in target_ids:
            matched_records.append(record)
            continue
        # Try fuzzy match (substring)
        for target in target_ids:
            if target in record.id:
                matched_records.append(record)
                break

    print(f"  [DEBUG] Matched {len(matched_records)} of {len(target_ids)} from {fasta_path.name}")

    if dry_run:
        return len(matched_records)

    with open(out_path, "w") as out_f:
        SeqIO.write(matched_records, out_f, "fasta")
    return len(matched_records)

def main():
    args = parse_args()

    summary_path = Path(args.summary_file)
    output_dir = Path(args.output_dir)
    base_dir = Path(args.base_dir)
    alignment_dir = Path(args.alignment_dir)
    key_path = Path(args.key_file)

    output_dir.mkdir(parents=True, exist_ok=True)
    generax_map = read_generax_key(key_path)

    processed = 0
    log_lines = []

    with open(summary_path) as summary_file:
        for line in summary_file:
            if args.limit is not None and processed >= args.limit:
                break
            if not line.strip():
                continue

            match = re.match(r"(\d+_generax)-Gene-(\d+):", line)
            if not match:
                continue

            og_name, gene_num = match.groups()
            gene_label = f"{og_name}-Gene-{gene_num}"
            gene_short = f"Gene{gene_num}"

            csv_path = base_dir / og_name / f"{og_name}_orthogroup.csv"
            if not csv_path.exists():
                log_lines.append(f"??  Missing .csv for {gene_label}")
                continue

            seq_names = []
            with open(csv_path) as csv_f:
                for row in csv_f:
                    if row.startswith(gene_label):
                        row_parts = row.strip().split(",")
                        seq_names = row_parts[1:]
                        break

            if not seq_names:
                log_lines.append(f"??  No line found for {gene_label} in {csv_path}")
                continue

            if og_name not in generax_map:
                log_lines.append(f"??  No mapping in GeneRaxKey for {og_name}")
                continue

            alignment_file = generax_map[og_name]
            alignment_path = alignment_dir / alignment_file
            if not alignment_path.exists():
                log_lines.append(f"??  Alignment file not found: {alignment_path}")
                continue

            gene_symbol = alignment_file.split("__")[0] if "__" in alignment_file else "UNKNOWN"
            out_name = f"{gene_symbol}__{gene_short}_{alignment_file}_orthocaller.fasta"
            out_path = output_dir / out_name

            n_written = extract_sequences_to_fasta(seq_names, alignment_path, out_path, dry_run=args.dry_run)

            log_lines.append(f"DRY-RUN ? {gene_label}: would extract {n_written} sequences")
            log_lines.append(f"         ? alignment: {alignment_path}")
            log_lines.append(f"         ? write to:  {out_path}")

            processed += 1

    with open("dryrun.txt", "w") as logf:
        for l in log_lines:
            logf.write(l + "\n")

if __name__ == "__main__":
    main()
