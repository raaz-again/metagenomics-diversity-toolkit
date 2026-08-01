#!/usr/bin/env python3
"""
diversity_metrics.py

Compute alpha diversity (Shannon, Simpson) for each sample and pairwise
beta diversity (Bray-Curtis dissimilarity) between samples, from an
OTU/ASV abundance table.

Expected input (tab-separated):
    OTU_ID    sample_1    sample_2    ...    sample_N

Usage:
    python diversity_metrics.py --table otu_table.tsv --alpha-output alpha_diversity.tsv --beta-output beta_diversity.tsv
"""

import argparse
import csv
import math
import sys
from itertools import combinations
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute alpha and beta diversity metrics from an OTU/ASV abundance table."
    )
    parser.add_argument("--table", required=True, help="TSV abundance table: OTU_ID plus one column per sample.")
    parser.add_argument("--alpha-output", required=True, help="Output TSV file with per-sample alpha diversity metrics.")
    parser.add_argument("--beta-output", required=True, help="Output TSV file with pairwise Bray-Curtis dissimilarities.")
    return parser.parse_args()


def load_table(path):
    """Load an OTU table into {sample_name: [counts...]} and a list of OTU ids."""
    with open(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        sample_names = header[1:]
        samples = {name: [] for name in sample_names}
        otu_ids = []
        for row in reader:
            if not row:
                continue
            otu_ids.append(row[0])
            for name, value in zip(sample_names, row[1:]):
                samples[name].append(float(value))
    return otu_ids, samples


def shannon_index(counts):
    total = sum(counts)
    if total == 0:
        return 0.0
    proportions = [c / total for c in counts if c > 0]
    return -sum(p * math.log(p) for p in proportions)


def simpson_index(counts):
    total = sum(counts)
    if total == 0:
        return 0.0
    proportions = [c / total for c in counts if c > 0]
    return 1 - sum(p ** 2 for p in proportions)


def bray_curtis(counts_a, counts_b):
    numerator = sum(abs(a - b) for a, b in zip(counts_a, counts_b))
    denominator = sum(a + b for a, b in zip(counts_a, counts_b))
    if denominator == 0:
        return 0.0
    return numerator / denominator


def write_alpha(samples, output_path):
    with open(output_path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "shannon", "simpson"])
        for name, counts in samples.items():
            writer.writerow([name, round(shannon_index(counts), 4), round(simpson_index(counts), 4)])


def write_beta(samples, output_path):
    names = list(samples.keys())
    with open(output_path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample_a", "sample_b", "bray_curtis"])
        for name_a, name_b in combinations(names, 2):
            dissimilarity = bray_curtis(samples[name_a], samples[name_b])
            writer.writerow([name_a, name_b, round(dissimilarity, 4)])


def main():
    args = parse_args()

    if not Path(args.table).exists():
        sys.exit(f"Abundance table not found: {args.table}")

    otu_ids, samples = load_table(args.table)
    write_alpha(samples, args.alpha_output)
    write_beta(samples, args.beta_output)

    print(f"Computed alpha diversity for {len(samples)} samples ({len(otu_ids)} OTUs). Written to {args.alpha_output}")
    print(f"Computed pairwise Bray-Curtis beta diversity. Written to {args.beta_output}")


if __name__ == "__main__":
    main()
