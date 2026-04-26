#!/usr/bin/env python3

from pathlib import Path
import os

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import gridspec


PRESENCE_FILE = Path(os.getenv("PRESENCE_FILE", "/scratch/p312334/project/10--HGT_isolates/1--General_infomation/function/result/function_sets/species_pair_ARG_CAZyme_Virulence_presence.csv"))
COMBO_FILE = Path(os.getenv("COMBO_FILE", "/scratch/p312334/project/10--HGT_isolates/1--General_infomation/function/result/function_sets/species_pair_ARG_CAZyme_Virulence_upset_counts.csv"))
PLOT_FILE = Path(os.getenv("PLOT_FILE", "/scratch/p312334/project/10--HGT_isolates/1--General_infomation/function/result/function_sets/species_pair_ARG_CAZyme_Virulence_upset.pdf"))

SETS = ["ARG", "CAZyme", "Virulence"]
ROW_BG = {
    "ARG": "#FBE3DE",
    "CAZyme": "#DDF2F6",
    "Virulence": "#D9F1E8",
}
BAR_GRAY = "#D9D9D9"
ABSENT_GRAY = "#D3D3D3"
PRESENT_BLACK = "#000000"


def main():
    presence = pd.read_csv(PRESENCE_FILE)
    combos = pd.read_csv(COMBO_FILE, keep_default_na=False)
    combos = combos[(combos["combination"] != "None") & (combos["combination"] != "")].copy()
    combos = combos.sort_values(
        ["species_pair_count", "ARG", "CAZyme", "Virulence"],
        ascending=[False, False, False, False],
    ).reset_index(drop=True)

    set_sizes = {name: int(presence[name].sum()) for name in SETS}

    plt.rcParams.update({
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.family": "DejaVu Sans",
        "font.size": 10,
    })

    fig = plt.figure(figsize=(7.2, 5.4), dpi=300)
    gs = gridspec.GridSpec(
        2,
        2,
        width_ratios=[5.2, 1.6],
        height_ratios=[3.2, 1.5],
        wspace=0.38,
        hspace=0.16,
    )

    ax_top = fig.add_subplot(gs[0, 0])
    ax_empty = fig.add_subplot(gs[0, 1])
    ax_mat = fig.add_subplot(gs[1, 0], sharex=ax_top)
    ax_right = fig.add_subplot(gs[1, 1])

    ax_empty.axis("off")

    x = list(range(len(combos)))
    counts = combos["species_pair_count"].tolist()

    ax_top.bar(x, counts, color=BAR_GRAY, width=0.56)
    for xi, yi in zip(x, counts):
        ax_top.text(xi, yi + max(counts) * 0.015, str(yi), ha="center", va="bottom", fontsize=10)
    ax_top.set_ylabel("Species-pair count")
    ax_top.set_xticks([])
    ax_top.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    ax_top.spines["top"].set_visible(False)
    ax_top.spines["right"].set_visible(False)
    ax_top.grid(False)

    y_positions = {name: len(SETS) - 1 - i for i, name in enumerate(SETS)}

    for name in SETS:
        y = y_positions[name]
        ax_mat.axhspan(y - 0.5, y + 0.5, color=ROW_BG[name], zorder=0)

    for xi, row in combos.iterrows():
        included = []
        for name in SETS:
            y = y_positions[name]
            present = int(row[name]) == 1
            ax_mat.scatter(xi, y, s=18, color=PRESENT_BLACK if present else ABSENT_GRAY, zorder=3)
            if present:
                included.append(y)
        if len(included) >= 2:
            ax_mat.plot([xi, xi], [min(included), max(included)], color=PRESENT_BLACK, linewidth=0.8, zorder=2)

    ax_mat.set_xlim(-0.5, len(combos) - 0.5)
    ax_mat.set_ylim(-0.5, len(SETS) - 0.5)
    ax_mat.set_yticks([y_positions[name] for name in SETS])
    ax_mat.set_yticklabels(SETS)
    ax_mat.set_xticks(x)
    ax_mat.set_xticklabels(combos["combination"], rotation=35, ha="right")
    ax_mat.spines["top"].set_visible(False)
    ax_mat.spines["right"].set_visible(False)
    ax_mat.spines["left"].set_visible(False)
    ax_mat.tick_params(axis="y", length=0)
    ax_mat.tick_params(axis="x", length=0)
    ax_mat.grid(False)

    right_order = list(reversed(SETS))
    right_vals = [set_sizes[name] for name in right_order]
    right_y = list(range(len(right_order)))

    ax_right.barh(right_y, right_vals, color=BAR_GRAY, height=0.56, zorder=1)
    for yi, val in zip(right_y, right_vals):
        ax_right.text(val * 0.06, yi, str(val), va="center", ha="left", fontsize=10)
    ax_right.set_yticks(right_y)
    ax_right.set_yticklabels([])
    ax_right.set_xlabel("Set size")
    ax_right.spines["top"].set_visible(False)
    ax_right.spines["right"].set_visible(False)
    ax_right.spines["left"].set_visible(False)
    ax_right.tick_params(axis="y", length=0)
    ax_right.grid(False)

    fig.suptitle("Transfer of ARG, CAZyme, and Virulence across species pairs", y=0.98, fontsize=12)
    fig.subplots_adjust(left=0.08, right=0.98, top=0.90, bottom=0.17)
    fig.savefig(PLOT_FILE, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
