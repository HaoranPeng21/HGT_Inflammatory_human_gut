#!/usr/bin/env python3
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import pandas as pd
from upsetplot import UpSet, from_indicators


BASE = Path("/scratch/p312334/project/10--HGT_isolates/1--General_infomation/Strain_transmission")
INFILE = BASE / "archive_20260409" / "family_eligible_between_individual_genome_pairs_HGT_ST_diff_species_regression_with_phylogeny.tsv"
OUTDIR = BASE / "figures_full_data"
OUTPNG = OUTDIR / "family_st100_hgt_upset_full_data.png"
OUTPDF = OUTDIR / "family_st100_hgt_upset_full_data.pdf"
OUTSVG = OUTDIR / "family_st100_hgt_upset_full_data.svg"


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(
        INFILE,
        sep="\t",
        usecols=["Family", "strain_transmission_ANI_eq_100", "HGT"],
    )
    df = df.rename(
        columns={
            "Family": "Same family",
            "strain_transmission_ANI_eq_100": "ST100",
            "HGT": "HGT",
        }
    )
    for col in df.columns:
        df[col] = df[col].astype(int).eq(1)

    total_n = len(df)
    combo_counts = (
        df.astype(int)
        .value_counts()
        .rename_axis(["Same family", "Strain transmission\n(ANI = 100)", "HGT"])
        .sort_values(ascending=False)
    )

    upset_data = from_indicators(df.columns.tolist(), data=df)

    plt.rcParams.update({
        "font.size": 12,
        "axes.titlesize": 15,
        "axes.labelsize": 12,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "figure.dpi": 300,
        "font.family": "Helvetica",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "pdf.use14corefonts": True,
    })
    fig = plt.figure(figsize=(11.5, 7.8))
    upset = UpSet(
        upset_data,
        subset_size="count",
        show_counts=True,
        sort_by="cardinality",
        sort_categories_by=None,
        facecolor="#9AA3AF",
        shading_color="#F5F6F8",
        other_dots_color="#D7DCE2",
        totals_plot_elements=3,
    )
    upset.plot(fig=fig)

    # Improve axis labels and formatting after the UpSet artists are created.
    for ax in fig.axes:
        if ax.get_ylabel() == "Intersection size":
            ax.set_ylabel("Genome-pair count")
            ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda x, _: f"{int(x):,}"))
            ax.grid(axis="y", color="#D9DEE7", linewidth=0.8)
            ax.set_axisbelow(True)
        elif ax.get_ylabel() == "Set size":
            ax.set_ylabel("Set size")
            ax.xaxis.set_major_formatter(mtick.FuncFormatter(lambda x, _: f"{int(x):,}"))
            ax.grid(axis="x", color="#D9DEE7", linewidth=0.8)
            ax.set_axisbelow(True)
        for spine in ax.spines.values():
            spine.set_color("#B8C2CC")

    subtitle = (
        f"different-species genome pairs, n = {total_n:,}; "
        f"largest all-zero group = {combo_counts.iloc[0]:,}"
    )
    fig.suptitle("Overlap of Same-Family, ST100, and HGT Genome Pairs", y=0.98, fontweight="bold")
    fig.text(0.5, 0.94, subtitle, ha="center", va="center", fontsize=11, color="#4C566A")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(OUTPNG, bbox_inches="tight")
    fig.savefig(OUTPDF, bbox_inches="tight")
    fig.savefig(OUTSVG, bbox_inches="tight")
    print(OUTPNG)
    print(OUTPDF)
    print(OUTSVG)


if __name__ == "__main__":
    main()
