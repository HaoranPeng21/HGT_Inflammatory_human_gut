#!/usr/bin/env python3

from pathlib import Path
import os

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import gridspec


INPUT = Path(os.getenv("INPUT_TABLE", "/scratch/p312334/project/10--HGT_isolates/data/all_within_individual_cross_species_genome_pairs_2000bp_hgt_flag_plasmid_func.csv"))
OUT_DIR = Path(os.getenv("OUT_DIR", "/scratch/p312334/project/10--HGT_isolates/1--General_infomation/function/result/function_sets"))
PRESENCE_OUT = OUT_DIR / "species_pair_ARG_CAZyme_Virulence_presence.csv"
COMBO_OUT = OUT_DIR / "species_pair_ARG_CAZyme_Virulence_upset_counts.csv"
PLOT_OUT = OUT_DIR / "species_pair_ARG_CAZyme_Virulence_upset.pdf"

COL_MAP = {
    "ARG": "Resfams_core_v1.2",
    "CAZyme": "CAZyme",
    "Virulence": "VFDB_setB_20260327",
}


def has_value(series: pd.Series) -> pd.Series:
    text = series.fillna("").astype(str).str.strip()
    return (~text.isin(["", "nan", "None"]))


def make_presence_table(df: pd.DataFrame) -> pd.DataFrame:
    work = df[["species_pair"] + list(COL_MAP.values())].copy()
    for out_col, in_col in COL_MAP.items():
        work[out_col] = has_value(work[in_col]).astype(int)

    presence = (
        work.groupby("species_pair")[list(COL_MAP.keys())]
        .max()
        .reset_index()
        .sort_values("species_pair")
    )
    return presence


def make_combo_table(presence: pd.DataFrame) -> pd.DataFrame:
    combo = (
        presence.groupby(list(COL_MAP.keys()))
        .size()
        .reset_index(name="species_pair_count")
        .sort_values(["species_pair_count", "ARG", "CAZyme", "Virulence"], ascending=[False, False, False, False])
        .reset_index(drop=True)
    )
    combo["combination"] = combo.apply(
        lambda r: "&".join([name for name in COL_MAP.keys() if r[name] == 1]) if any(r[name] == 1 for name in COL_MAP.keys()) else "None",
        axis=1,
    )
    return combo


def plot_upset(combo: pd.DataFrame, presence: pd.DataFrame) -> None:
    plot_df = combo[combo["species_pair_count"] > 0].copy()
    plot_df["x"] = range(len(plot_df))

    set_sizes = presence[list(COL_MAP.keys())].sum().reindex(list(COL_MAP.keys()))

    fig = plt.figure(figsize=(7.2, 5.2))
    gs = gridspec.GridSpec(
        2,
        2,
        width_ratios=[1.7, 4.8],
        height_ratios=[3.0, 1.6],
        wspace=0.08,
        hspace=0.05,
    )

    ax_blank = fig.add_subplot(gs[0, 0])
    ax_blank.axis("off")

    ax_bar = fig.add_subplot(gs[0, 1])
    ax_set = fig.add_subplot(gs[1, 0])
    ax_mat = fig.add_subplot(gs[1, 1], sharex=ax_bar)

    ax_bar.bar(plot_df["x"], plot_df["species_pair_count"], color="#4DBBD5FF", width=0.72)
    ax_bar.set_ylabel("Species-pair count")
    ax_bar.set_xticks([])
    ax_bar.spines["top"].set_visible(False)
    ax_bar.spines["right"].set_visible(False)

    y_pos = list(range(len(COL_MAP)))[::-1]
    names = list(COL_MAP.keys())
    ax_set.barh(y_pos, set_sizes.values, color="#B09C85FF", height=0.62)
    ax_set.set_yticks(y_pos)
    ax_set.set_yticklabels(names)
    ax_set.set_xlabel("Set size")
    ax_set.spines["top"].set_visible(False)
    ax_set.spines["right"].set_visible(False)

    ax_mat.set_ylim(-0.5, len(COL_MAP) - 0.5)
    ax_mat.set_yticks(y_pos)
    ax_mat.set_yticklabels([])
    ax_mat.set_xlabel("Intersection")
    ax_mat.spines["top"].set_visible(False)
    ax_mat.spines["right"].set_visible(False)
    ax_mat.spines["left"].set_visible(False)
    ax_mat.tick_params(axis="y", length=0)

    for _, row in plot_df.iterrows():
        x = row["x"]
        included_y = []
        for i, name in enumerate(names):
            y = len(names) - 1 - i
            present = row[name] == 1
            ax_mat.scatter(
                [x],
                [y],
                s=48,
                color="black" if present else "#d9d9d9",
                zorder=3,
            )
            if present:
                included_y.append(y)
        if len(included_y) >= 2:
            ax_mat.plot([x, x], [min(included_y), max(included_y)], color="black", linewidth=1.2, zorder=2)

    ax_mat.set_xticks(plot_df["x"])
    ax_mat.set_xticklabels(plot_df["combination"], rotation=45, ha="right")

    plt.tight_layout()
    fig.savefig(PLOT_OUT)
    plt.close(fig)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(INPUT, usecols=["species_pair"] + list(COL_MAP.values()), low_memory=False)

    presence = make_presence_table(df)
    presence.to_csv(PRESENCE_OUT, index=False)

    combo = make_combo_table(presence)
    combo.to_csv(COMBO_OUT, index=False)

    plot_upset(combo, presence)

    print(f"presence={PRESENCE_OUT}")
    print(f"combinations={COMBO_OUT}")
    print(f"plot={PLOT_OUT}")
    print(f"species_pairs={len(presence)}")
    for name in COL_MAP:
        print(f"{name}={int(presence[name].sum())}")


if __name__ == "__main__":
    main()
