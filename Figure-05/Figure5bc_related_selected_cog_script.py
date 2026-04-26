#!/usr/bin/env python3

from __future__ import annotations

import textwrap
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd


BASE = Path(
    "/scratch/p312334/project/10--HGT_isolates/4--Function_IBD_association/"
    "COG_function_secondary_cluster"
)
FOCUS = BASE / "robust_category_focus"
OUT_PDF = BASE / "selected_C_H_COGs_ibd_fisher_fdr05_lollipop.pdf"
OUT_CSV = BASE / "selected_C_H_COGs_ibd_fisher_fdr05_lollipop.csv"

TARGETS = {
    "C": ["COG1853", "COG2764"],
    "H": ["COG0294", "COG0156", "COG4206", "COG1541"],
}

CAT_COLORS = {
    "C": "#D73027",
    "H": "#1A9850",
}


def wrap_label(value: str, width: int) -> str:
    return "\n".join(textwrap.wrap(str(value).replace("_", " "), width=width))


def load_selected() -> pd.DataFrame:
    frames = []
    for category, cogs in TARGETS.items():
        path = FOCUS / category / "cog_function_species_pair_ibd_fisher.csv"
        df = pd.read_csv(path)
        df = df[df["COG24_FUNCTION_ID"].isin(cogs)].copy()
        df = df[df["fisher_q_bh"] < 0.05].copy()
        df["target_category"] = category
        frames.append(df)

    out = pd.concat(frames, ignore_index=True)
    target_order = [cog for cogs in TARGETS.values() for cog in cogs]
    out["COG24_FUNCTION_ID"] = pd.Categorical(
        out["COG24_FUNCTION_ID"], categories=target_order, ordered=True
    )
    out = out.sort_values(["COG24_FUNCTION_ID", "log2_or_haldane"], ascending=[True, True])
    return out


def main() -> None:
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42
    mpl.rcParams["font.family"] = "DejaVu Sans"

    df = load_selected()
    df.to_csv(OUT_CSV, index=False)

    cogs = [cog for cog in df["COG24_FUNCTION_ID"].cat.categories if (df["COG24_FUNCTION_ID"] == cog).any()]
    fig, axes = plt.subplots(
        nrows=3,
        ncols=2,
        figsize=(12.5, 8.2),
        sharex=True,
        constrained_layout=False,
    )
    axes = axes.flatten()

    xmax = max(10.2, float(df["log2_or_haldane"].max()) + 0.7)
    xmin = min(-0.5, float(df["log2_or_haldane"].min()) - 0.4)

    for ax, cog in zip(axes, cogs):
        sub = df[df["COG24_FUNCTION_ID"] == cog].copy()
        sub = sub.sort_values("log2_or_haldane", ascending=True).reset_index(drop=True)
        y = range(len(sub))
        category = str(sub["target_category"].iloc[0])
        color = CAT_COLORS[category]

        ax.axvline(0, color="#888888", linestyle="--", linewidth=0.6, zorder=0)
        ax.hlines(y, 0, sub["log2_or_haldane"], color=color, linewidth=1.25, alpha=0.85)
        ax.scatter(
            sub["log2_or_haldane"],
            list(y),
            s=42,
            color=color,
            edgecolor="white",
            linewidth=0.5,
            zorder=3,
        )

        ax.set_yticks(list(y))
        ax.set_yticklabels([wrap_label(x, 34) for x in sub["species_pair_name"]], fontsize=9.2)
        ax.set_xlim(xmin, xmax)
        ax.grid(axis="x", color="#E6E6E6", linewidth=0.6)
        ax.tick_params(axis="x", labelsize=9.5)
        ax.tick_params(axis="y", length=0)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        function_name = str(sub["COG24_FUNCTION_NAME"].iloc[0])
        title = f"{category} | {cog}\n{wrap_label(function_name, 48)}"
        ax.set_title(title, fontsize=10.8, fontweight="bold", pad=8)

    for ax in axes[len(cogs):]:
        ax.axis("off")

    fig.suptitle(
        "Selected COG Functions: IBD-Enriched Species-Pair Fisher Tests",
        fontsize=16,
        fontweight="bold",
        y=0.985,
    )
    fig.supxlabel("log2 odds ratio (IBD vs non-IBD), Haldane corrected", fontsize=12.5, y=0.035)

    handles = [
        plt.Line2D([0], [0], color=CAT_COLORS["C"], marker="o", linestyle="-", label="C: Energy production and conversion"),
        plt.Line2D([0], [0], color=CAT_COLORS["H"], marker="o", linestyle="-", label="H: Coenzyme transport and metabolism"),
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.0),
        ncol=2,
        frameon=False,
        fontsize=10.5,
    )

    fig.subplots_adjust(left=0.24, right=0.985, top=0.88, bottom=0.11, wspace=0.42, hspace=0.78)
    fig.savefig(OUT_PDF, format="pdf")
    plt.close(fig)

    print(f"Wrote {OUT_PDF}")
    print(f"Wrote {OUT_CSV}")


if __name__ == "__main__":
    main()
