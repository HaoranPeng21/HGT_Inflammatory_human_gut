#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
import math
import re

import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Polygon
import pandas as pd


BASE = Path("/scratch/p312334/project/10--HGT_isolates/2--Phenotype_association/comparision_plot/Streptococcus parasanguinis D | Bifidobacterium longum")
FEATURES = BASE / "233_3-338_1_plasmid_exactdedup_feature_annotation.tsv"
REPS = BASE / "233_3-338_1_plasmid_exactdedup_representatives.csv"
OUT_PDF = BASE / "233_3-338_1_plasmid_exactdedup_circular_pretty.pdf"
OUT_PNG = BASE / "233_3-338_1_plasmid_exactdedup_circular_pretty.png"


COLOR_MAP = {
    "RepA": "#c0392b",
    "putative replication protein rep": "#d35400",
    "MobA": "#1f78b4",
    "Mobilization protein MobA": "#2a9d8f",
    "Conjugal transfer protein TraD": "#6a4c93",
    "DUF4328 domain-containing protein": "#f4a261",
    "Transmembrane protein": "#4d908e",
    "Memba protein": "#90be6d",
}

LABEL_MAP = {
    "RepA": "RepA",
    "putative replication protein rep": "rep",
    "MobA": "MobA",
    "Mobilization protein MobA": "MobA",
    "Conjugal transfer protein TraD": "TraD",
    "DUF4328 domain-containing protein": "DUF4328",
    "Transmembrane protein": "TM",
    "Memba protein": "Memba",
}


def parse_coords(text: str) -> tuple[int, int, int]:
    start, end, strand = text.split(":")
    return int(start), int(end), int(strand)


def polar_to_xy(radius: float, angle_deg: float) -> tuple[float, float]:
    ang = math.radians(angle_deg)
    return radius * math.cos(ang), radius * math.sin(ang)


def short_title(label: str) -> str:
    m = re.match(r"\d+_Module_([A-Z])_(\d+)bp_n(\d+)", label)
    if not m:
        return label
    mod, bp, n = m.groups()
    return f"Module {mod} ({bp} bp)"


def subtitle(label: str) -> str:
    m = re.match(r"\d+_Module_([A-Z])_(\d+)bp_n(\d+)", label)
    if not m:
        return ""
    return f"n={m.group(3)} exact plasmid-side copies"


def draw_gene(ax, start: int, end: int, strand: int, length: int, color: str, r_outer: float = 1.0, width: float = 0.24):
    theta1 = 90 - (start - 1) / length * 360
    theta2 = 90 - end / length * 360
    if theta2 > theta1:
        theta1, theta2 = theta2, theta1

    body_frac = 0.86
    delta = theta1 - theta2
    head = max(min(delta * (1 - body_frac), 10), 4)

    if strand == 1:
        body_theta1, body_theta2 = theta1, theta2 + head
        tip_theta = theta2
    else:
        body_theta1, body_theta2 = theta1 - head, theta2
        tip_theta = theta1

    ax.add_patch(
        Wedge((0, 0), r_outer, body_theta2, body_theta1, width=width, facecolor=color, edgecolor="white", lw=1.8)
    )

    r_mid = r_outer - width / 2
    r_in = r_outer - width
    r_out = r_outer
    if strand == 1:
        p1 = polar_to_xy(r_in, body_theta2)
        p2 = polar_to_xy(r_out, body_theta2)
    else:
        p1 = polar_to_xy(r_in, body_theta1)
        p2 = polar_to_xy(r_out, body_theta1)
    p3 = polar_to_xy(r_mid, tip_theta)
    ax.add_patch(Polygon([p1, p2, p3], closed=True, facecolor=color, edgecolor="white", lw=1.8))


def draw_module(ax, rep_label: str, length: int, n_members: int, sub: pd.DataFrame) -> None:
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_xlim(-1.55, 1.55)
    ax.set_ylim(-1.45, 1.45)

    ax.add_patch(Wedge((0, 0), 1.04, 0, 360, width=0.03, facecolor="#d9d9d9", edgecolor="none"))

    for _, row in sub.iterrows():
        start, end, strand = parse_coords(row["coordinates"])
        color = COLOR_MAP.get(row["name"], "#bbbbbb")
        draw_gene(ax, start, end, strand, length, color)

        mid = (start + end) / 2
        angle = 90 - mid / length * 360
        rx, ry = polar_to_xy(1.22, angle)
        label = LABEL_MAP.get(row["name"], row["name"])
        ha = "left" if rx >= 0 else "right"
        ax.text(rx, ry, label, ha=ha, va="center", fontsize=10.5, fontweight="bold", color=color)

    ax.text(0, 0.22, short_title(rep_label), ha="center", va="center", fontsize=16, fontweight="bold")
    ax.text(0, 0.02, subtitle(rep_label), ha="center", va="center", fontsize=11.5, color="#444444")
    ax.text(0, -0.18, "plasmid-side exact dedup representative", ha="center", va="center", fontsize=10, color="#666666")


def draw_legend(fig):
    items = [
        ("Replication", COLOR_MAP["RepA"]),
        ("Mobilization", COLOR_MAP["MobA"]),
        ("Conjugation", COLOR_MAP["Conjugal transfer protein TraD"]),
        ("Accessory / membrane", COLOR_MAP["DUF4328 domain-containing protein"]),
    ]
    x = 0.14
    y = 0.04
    for label, color in items:
        fig.text(x, y, "■", color=color, fontsize=14, va="center")
        fig.text(x + 0.018, y, label, fontsize=10.5, va="center", color="#333333")
        x += 0.19


def main() -> None:
    feat = pd.read_csv(FEATURES, sep="\t")
    reps = pd.read_csv(REPS)

    fig, axes = plt.subplots(2, 1, figsize=(9.5, 13))
    fig.patch.set_facecolor("white")

    for ax, rep in zip(axes, reps.itertuples(index=False)):
        sub = feat[feat["locus_id"] == rep.rep_label].copy()
        draw_module(ax, rep.rep_label, int(rep.length), int(rep.n_members), sub)

    draw_legend(fig)
    plt.tight_layout(rect=(0.02, 0.06, 0.98, 0.98), h_pad=2.6)
    fig.savefig(OUT_PDF, bbox_inches="tight", facecolor=fig.get_facecolor())
    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())


if __name__ == "__main__":
    main()
