#!/usr/bin/env python3
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator, StrMethodFormatter


BASE = Path("/scratch/p312334/project/10--HGT_isolates/1--General_infomation/conserve_cutoff")
REPRESENTATIVES = BASE / "secondary_cluster_best_genomes.tsv"
SNP_LONG = BASE / "secondary_cluster_best.res" / "ribosomal42_snp" / "ribosomal42_pairwise_snp_long.tsv"
V3 = Path("/scratch/p312334/project/10--HGT_isolates/data/__binstoget_isolates_summary_v3.csv")
OUTDIR = BASE / "result"
PREFIX = "ALL_speciespair_ribosomal42_snp"


def canonical_pair(a: str, b: str) -> tuple[str, str]:
    return (a, b) if a <= b else (b, a)


def read_v3_classification() -> dict[str, str]:
    out = {}
    with V3.open(newline="") as handle:
        for row in csv.DictReader(handle):
            genome = row["Genome_file"].strip()
            classification = row.get("classification", "").strip()
            if genome and classification and genome not in out:
                out[genome] = classification
    return out


def normalize_species(species: str, classification: str) -> str:
    species = species.strip()
    if species and species != "s__":
        return species
    parts = [p.strip() for p in classification.split(";") if p.strip()]
    for part in reversed(parts):
        if part.startswith("s__") and part != "s__":
            return part
    return species or "s__unknown"


def read_representatives(v3_classification: dict[str, str]) -> dict[str, dict[str, str]]:
    reps = {}
    with REPRESENTATIVES.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            genome = row["Genome_file"].strip()
            genome_with_suffix = f"{genome}.unicycler"
            species = normalize_species(row.get("Species", ""), v3_classification.get(genome, ""))
            reps[genome_with_suffix] = {
                "secondary_cluster": row["secondary_cluster"].strip(),
                "genome": genome,
                "genome_with_suffix": genome_with_suffix,
                "species": species,
            }
    return reps


def save_integer_zoom_histogram(
    values: np.ndarray,
    out_prefix: Path,
    ylabel: str,
    title: str,
    max_value: int,
) -> None:
    zoom_values = values[(values >= 0) & (values <= max_value)]
    zoom_bins = np.arange(-0.5, max_value + 1.5, 1)

    fig, ax = plt.subplots(figsize=(6.2, 5.5))
    ax.hist(zoom_values, bins=zoom_bins, color="#8A8A8A", edgecolor="#6A6A6A", linewidth=0.8, alpha=1.0)
    ax.set_xlim(-0.5, max_value + 0.5)
    ax.set_xticks(range(0, max_value + 1))
    ax.xaxis.set_major_formatter(StrMethodFormatter("{x:.0f}"))
    ax.set_axisbelow(True)
    ax.grid(axis="y", color="#E3E3E3", linewidth=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#777777")
    ax.spines["bottom"].set_color("#777777")
    ax.tick_params(axis="both", colors="#4F4F4F", labelsize=10)
    ax.set_xlabel("Pairwise SNP per 2000 bp")
    ax.set_ylabel(ylabel)
    ax.set_title(f"{title} (0-{max_value})")
    fig.tight_layout()
    fig.savefig(out_prefix.with_suffix(".pdf"))
    fig.savefig(out_prefix.with_suffix(".png"), dpi=200)
    plt.close(fig)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    v3_classification = read_v3_classification()
    reps = read_representatives(v3_classification)

    detailed_rows: list[dict[str, object]] = []
    skipped_nonrep = 0

    with SNP_LONG.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            genome_1 = row["genome_1"].strip()
            genome_2 = row["genome_2"].strip()
            rep1 = reps.get(genome_1)
            rep2 = reps.get(genome_2)
            if rep1 is None or rep2 is None:
                skipped_nonrep += 1
                continue

            comparable = int(row["comparable_sites"])
            snps = int(row["snp_count"])
            snp_rate = float(row["snp_rate"]) if row["snp_rate"] else np.nan
            snp_per_2000bp = (
                float(row["snp_per_2000bp"])
                if row.get("snp_per_2000bp")
                else (snps / comparable * 2000 if comparable else np.nan)
            )

            species_1 = rep1["species"]
            species_2 = rep2["species"]
            sp1, sp2 = canonical_pair(species_1, species_2)
            cl1, cl2 = canonical_pair(rep1["secondary_cluster"], rep2["secondary_cluster"])
            detailed_rows.append(
                {
                    "secondary_cluster_pair": f"{cl1}-{cl2}",
                    "secondary_cluster_1": cl1,
                    "secondary_cluster_2": cl2,
                    "species_1": species_1,
                    "species_2": species_2,
                    "species_pair": f"{sp1} | {sp2}",
                    "genome_1": genome_1,
                    "genome_2": genome_2,
                    "snp_count": snps,
                    "comparable_sites": comparable,
                    "snp_rate": snp_rate,
                    "snp_per_2000bp": snp_per_2000bp,
                }
            )

    detailed_path = OUTDIR / f"{PREFIX}.tsv"
    with detailed_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "secondary_cluster_pair",
                "secondary_cluster_1",
                "secondary_cluster_2",
                "species_1",
                "species_2",
                "species_pair",
                "genome_1",
                "genome_2",
                "snp_count",
                "comparable_sites",
                "snp_rate",
                "snp_per_2000bp",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(detailed_rows)

    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in detailed_rows:
        grouped[str(row["species_pair"])].append(row)

    summary_rows = []
    for species_pair, rows in sorted(grouped.items(), key=lambda item: (-len(item[1]), item[0])):
        snps = np.array([int(r["snp_count"]) for r in rows], dtype=float)
        rates = np.array([float(r["snp_rate"]) for r in rows], dtype=float)
        per2000 = np.array([float(r["snp_per_2000bp"]) for r in rows], dtype=float)
        summary_rows.append(
            {
                "species_pair": species_pair,
                "n_cluster_pairs": len(rows),
                "median_snp_count": float(np.median(snps)),
                "mean_snp_count": float(np.mean(snps)),
                "min_snp_count": int(np.min(snps)),
                "max_snp_count": int(np.max(snps)),
                "median_snp_rate": float(np.nanmedian(rates)),
                "mean_snp_rate": float(np.nanmean(rates)),
                "median_snp_per_2000bp": float(np.nanmedian(per2000)),
                "mean_snp_per_2000bp": float(np.nanmean(per2000)),
                "min_snp_per_2000bp": float(np.nanmin(per2000)),
                "max_snp_per_2000bp": float(np.nanmax(per2000)),
            }
        )

    summary_path = OUTDIR / f"{PREFIX}_summary.tsv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)

    stats_path = OUTDIR / f"{PREFIX}_stats.txt"
    snp_values = np.array([int(r["snp_count"]) for r in detailed_rows], dtype=float)
    rate_values = np.array([float(r["snp_rate"]) for r in detailed_rows], dtype=float)
    per2000_values = np.array([float(r["snp_per_2000bp"]) for r in detailed_rows], dtype=float)
    with stats_path.open("w") as handle:
        handle.write(f"pairwise_genome_pairs_in_input\t{sum(1 for _ in SNP_LONG.open()) - 1}\n")
        handle.write(f"pairwise_genome_pairs_used\t{len(detailed_rows)}\n")
        handle.write(f"unique_species_pairs\t{len(grouped)}\n")
        handle.write(f"skipped_nonrepresentative_pairs\t{skipped_nonrep}\n")
        handle.write(f"snp_count_median\t{np.median(snp_values):.2f}\n")
        handle.write(f"snp_count_mean\t{np.mean(snp_values):.2f}\n")
        handle.write(f"snp_count_min\t{int(np.min(snp_values))}\n")
        handle.write(f"snp_count_max\t{int(np.max(snp_values))}\n")
        handle.write(f"snp_rate_median\t{np.nanmedian(rate_values):.6f}\n")
        handle.write(f"snp_rate_mean\t{np.nanmean(rate_values):.6f}\n")
        handle.write(f"snp_per_2000bp_median\t{np.nanmedian(per2000_values):.6f}\n")
        handle.write(f"snp_per_2000bp_mean\t{np.nanmean(per2000_values):.6f}\n")
        handle.write(f"snp_per_2000bp_min\t{np.nanmin(per2000_values):.6f}\n")
        handle.write(f"snp_per_2000bp_max\t{np.nanmax(per2000_values):.6f}\n")

    plt.style.use("default")

    max_per2000 = float(np.nanmax(per2000_values))
    histogram_bins = np.arange(0, max_per2000 + 1, 1)
    if histogram_bins[-1] <= max_per2000:
        histogram_bins = np.append(histogram_bins, histogram_bins[-1] + 1)

    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.hist(per2000_values, bins=histogram_bins, color="black", edgecolor="black", alpha=1.0)
    ax.set_xlabel("Pairwise SNP per 2000 bp")
    ax.set_ylabel("Number of Species Pairs")
    ax.set_title("Ribosomal42 SNP per 2000 bp Across All Species Pairs")
    fig.tight_layout()
    fig.savefig(OUTDIR / f"{PREFIX}_histogram.pdf")
    fig.savefig(OUTDIR / f"{PREFIX}_histogram.png", dpi=200)
    plt.close(fig)

    zoom_values = per2000_values[(per2000_values >= 0) & (per2000_values <= 20)]
    zoom_bins = np.arange(0, 20 + 1, 1)
    fig, ax = plt.subplots(figsize=(7.2, 5.5))
    ax.hist(zoom_values, bins=zoom_bins, color="black", edgecolor="black", alpha=1.0, rwidth=0.62)
    ax.set_xlim(0, 20)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.xaxis.set_major_formatter(StrMethodFormatter("{x:.0f}"))
    ax.set_xlabel("Pairwise SNP per 2000 bp")
    ax.set_ylabel("Number of Species Pairs")
    ax.set_title("Ribosomal42 SNP per 2000 bp Across All Species Pairs (0-20)")
    fig.tight_layout()
    fig.savefig(OUTDIR / f"{PREFIX}_histogram_0_20.pdf")
    fig.savefig(OUTDIR / f"{PREFIX}_histogram_0_20.png", dpi=200)
    plt.close(fig)

    save_integer_zoom_histogram(
        per2000_values,
        OUTDIR / f"{PREFIX}_histogram_0_10",
        ylabel="Number of Species Pairs",
        title="Ribosomal42 SNP per 2000 bp Across All Species Pairs",
        max_value=10,
    )

    print(f"[INFO] Wrote detailed rows: {detailed_path}")
    print(f"[INFO] Wrote summary: {summary_path}")
    print(f"[INFO] Wrote stats: {stats_path}")
    print(f"[INFO] Representative genome pairs used: {len(detailed_rows)}")
    print(f"[INFO] Unique species pairs: {len(grouped)}")


if __name__ == "__main__":
    main()
