#!/usr/bin/env python3

import csv
from collections import defaultdict
from pathlib import Path


SOURCE_TABLE = Path("/scratch/p312334/project/10--HGT_isolates/data/HGT_table_2000bp_isolates_within_plasmid_func.csv")
TARGET_TABLE = Path("/scratch/p312334/project/10--HGT_isolates/data/all_within_individual_cross_species_genome_pairs_2000bp_hgt_flag.csv")
OUTPUT_TABLE = Path("/scratch/p312334/project/10--HGT_isolates/data/all_within_individual_cross_species_genome_pairs_2000bp_hgt_flag_func.csv")

SOURCE_COLUMNS = [
    "COG24_FUNCTION",
    "COG24_CATEGORY",
    "Pfam",
    "COG24_PATHWAY",
    "Resfams_core_v1.2",
    "CAZyme",
    "VFDB_setB_20260327",
]


def canonical_pair(a: str, b: str):
    return tuple(sorted((a, b)))


def append_unique_tokens(token_string: str, seen: set, out: list):
    if not token_string:
        return
    for token in token_string.split("|"):
        token = token.strip()
        if not token or token in seen:
            continue
        seen.add(token)
        out.append(token)


def build_pair_annotations():
    aggregated = defaultdict(lambda: {col: [] for col in SOURCE_COLUMNS})
    seen_tokens = defaultdict(lambda: {col: set() for col in SOURCE_COLUMNS})

    with SOURCE_TABLE.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            pair = canonical_pair(row["MAG_1"], row["MAG_2"])
            for col in SOURCE_COLUMNS:
                append_unique_tokens(
                    row.get(col, ""),
                    seen_tokens[pair][col],
                    aggregated[pair][col],
                )

    return aggregated


def main():
    aggregated = build_pair_annotations()

    with TARGET_TABLE.open(newline="") as in_handle, OUTPUT_TABLE.open("w", newline="") as out_handle:
        reader = csv.DictReader(in_handle)
        fieldnames = list(reader.fieldnames) + SOURCE_COLUMNS
        writer = csv.DictWriter(out_handle, fieldnames=fieldnames)
        writer.writeheader()

        for row in reader:
            pair = canonical_pair(row["genome_1"], row["genome_2"])
            annotations = aggregated.get(pair, {})
            for col in SOURCE_COLUMNS:
                row[col] = "|".join(annotations.get(col, []))
            writer.writerow(row)

    print(f"output={OUTPUT_TABLE}")
    print(f"genome_pairs_with_annotations={len(aggregated)}")


if __name__ == "__main__":
    main()
