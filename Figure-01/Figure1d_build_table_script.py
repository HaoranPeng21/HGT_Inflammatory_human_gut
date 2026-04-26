#!/usr/bin/env python3

import csv
import os
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt


INPUT_TABLE = Path(os.getenv("INPUT_TABLE", "/scratch/p312334/project/10--HGT_isolates/data/all_within_individual_cross_species_genome_pairs_2000bp_hgt_flag_plasmid_func.csv"))
OUT_DIR = Path(os.getenv("OUT_DIR", "/scratch/p312334/project/10--HGT_isolates/1--General_infomation/function/result"))
MATRIX_OUT = OUT_DIR / "species_pair_by_COG24_CATEGORY_presence.csv"
FREQ_OUT = OUT_DIR / "COG24_CATEGORY_species_pair_frequency.csv"
PLOT_OUT = OUT_DIR / "COG24_CATEGORY_species_pair_frequency_barplot.png"
MIN_GENOME_PAIRS = 10


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    species_pair_letters = defaultdict(set)
    all_species_pairs = set()
    all_letters = set()
    species_pair_genome_pairs = defaultdict(set)

    with INPUT_TABLE.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            species_pair = row["species_pair"]
            all_species_pairs.add(species_pair)
            species_pair_genome_pairs[species_pair].add(row["genome_pair"])
            category_string = row.get("COG24_CATEGORY", "").strip()
            if not category_string:
                continue
            for letter in category_string.split("|"):
                letter = letter.strip()
                if not letter:
                    continue
                species_pair_letters[species_pair].add(letter)
                all_letters.add(letter)

    eligible_species_pairs = sorted(
        sp for sp in all_species_pairs if len(species_pair_genome_pairs[sp]) >= MIN_GENOME_PAIRS
    )
    ordered_letters = sorted(all_letters)

    with MATRIX_OUT.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["species_pair"] + ordered_letters)
        for species_pair in eligible_species_pairs:
            present = species_pair_letters.get(species_pair, set())
            writer.writerow([species_pair] + [1 if letter in present else 0 for letter in ordered_letters])

    total_species_pairs = len(eligible_species_pairs)
    frequencies = []
    for letter in ordered_letters:
        count = sum(1 for species_pair in eligible_species_pairs if letter in species_pair_letters.get(species_pair, set()))
        freq = count / total_species_pairs if total_species_pairs else 0
        frequencies.append((letter, count, total_species_pairs, freq))

    frequencies.sort(key=lambda x: (-x[3], x[0]))

    with FREQ_OUT.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["COG24_CATEGORY", "species_pair_count_with_1", "species_pair_total", "frequency"])
        writer.writerows(frequencies)

    letters = [x[0] for x in frequencies]
    freqs = [x[3] for x in frequencies]

    plt.figure(figsize=(12, 5))
    bars = plt.bar(letters, freqs, color="#4C7A9F", edgecolor="black", linewidth=0.5)
    plt.xlabel("COG24_CATEGORY")
    plt.ylabel("Frequency")
    plt.title("Frequency of COG24 categories across species pairs")
    plt.xticks(rotation=0)
    plt.tight_layout()
    plt.savefig(PLOT_OUT, dpi=300)
    plt.close()

    print(f"matrix={MATRIX_OUT}")
    print(f"frequency={FREQ_OUT}")
    print(f"plot={PLOT_OUT}")
    print(f"species_pairs_total={len(all_species_pairs)}")
    print(f"species_pairs_with_at_least_{MIN_GENOME_PAIRS}_genome_pairs={total_species_pairs}")
    print(f"categories={len(ordered_letters)}")


if __name__ == "__main__":
    main()
