#!/usr/bin/env python3
"""
Figure 1: assay x genome mismatch heatmap.

Reads step5_output/assay_coverage_per_genome.csv (worst_mm per assay per
genome, plus classification) and draws one cell per assay/genome pair,
colour = worst-component mismatch count. NOT EVALUABLE cells (binding site
too N-rich to call - worst_mm is not a real 0, it is "nothing scored") are
rendered as a distinct muted grey, never blended into the mismatch scale.

Palette: validated sequential blue ramp and status/categorical hexes from
the project's dataviz reference palette (single-hue sequential for
magnitude, fixed-order categorical for the outbreak-arm strip, muted grey
for the not-evaluable class) - not eyeballed.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

OUT = Path("figures")
OUT.mkdir(exist_ok=True)

# ---- palette (from dataviz reference palette.md) --------------------------
SEQ_BLUE = ["#cde2fb", "#9ec5f4", "#5598e7", "#2a78d6", "#184f95", "#0d366b"]
NOT_EVAL_GREY = "#c3c2b7"
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRID = "#e1e0d9"
ARM_COLORS = {"Bdbv-2007": "#eb6834", "Bdbv-2012": "#1baf7a", "Bdbv-2026": "#eda100"}  # cat slots 2/3/4, blue reserved for the sequential fill
SURFACE = "#fcfcfb"

ASSAY_ORDER = [
    "Greninger2024_PfiloL",
    "Ogawa2010_EBOV",
    "Panning2007_panfilo",
    "SanchezFiloAB",
    "Weidmann2004_ENZ",
    "Weidmann2004_ENS",
    "Weidmann2004_MN",
    "Panning2007_MARVprobe_control",
]
ASSAY_LABELS = {
    "Greninger2024_PfiloL": "Greninger 2024 (pan-filovirus L)",
    "Ogawa2010_EBOV": "Ogawa 2010 (EBOV NP)",
    "Panning2007_panfilo": "Panning 2007 (pan-filovirus L)",
    "SanchezFiloAB": "Sanchez 1999 (pan-filovirus L)",
    "Weidmann2004_ENZ": "Weidmann 2004 ENZ (EBOV-Zaire NP)",
    "Weidmann2004_ENS": "Weidmann 2004 ENS (EBOV-Sudan NP)",
    "Weidmann2004_MN": "Weidmann 2004 MN (Marburg NP, control)",
    "Panning2007_MARVprobe_control": "Panning MARV probe (control)",
}
ARM_ORDER = ["Bdbv-2007", "Bdbv-2012", "Bdbv-2026"]
ARM_LABELS = {"Bdbv-2007": "2007 (n=12)", "Bdbv-2012": "2012 (n=23)", "Bdbv-2026": "2026 (n=35)"}


def main():
    df = pd.read_csv("step5_output/assay_coverage_per_genome.csv")

    # column order: by arm, then by genome_id within arm
    df["outbreak"] = pd.Categorical(df["outbreak"], categories=ARM_ORDER, ordered=True)
    genome_order = (
        df[["genome_id", "outbreak"]]
        .drop_duplicates()
        .sort_values(["outbreak", "genome_id"])["genome_id"]
        .tolist()
    )

    pivot_mm = df.pivot(index="assay", columns="genome_id", values="worst_mm")
    pivot_cls = df.pivot(index="assay", columns="genome_id", values="classification")
    pivot_mm = pivot_mm.reindex(index=ASSAY_ORDER, columns=genome_order)
    pivot_cls = pivot_cls.reindex(index=ASSAY_ORDER, columns=genome_order)

    n_assays, n_genomes = pivot_mm.shape
    vmax = int(np.nanmax(pivot_mm.values[~(pivot_cls.values == "NOT EVALUABLE")]))

    cmap = LinearSegmentedColormap.from_list("seq_blue", SEQ_BLUE, N=256)

    fig_w = max(10, n_genomes * 0.11)
    fig_h = 3.6
    fig, (ax_strip, ax) = plt.subplots(
        2, 1, figsize=(fig_w, fig_h + 0.35),
        gridspec_kw={"height_ratios": [0.35, n_assays], "hspace": 0.04},
    )
    fig.patch.set_facecolor(SURFACE)

    # outbreak-arm strip
    arm_row = df[["genome_id", "outbreak"]].drop_duplicates().set_index("genome_id").loc[genome_order, "outbreak"]
    for i, gid in enumerate(genome_order):
        ax_strip.add_patch(plt.Rectangle((i, 0), 1, 1, facecolor=ARM_COLORS[arm_row[gid]], edgecolor="none"))
    ax_strip.set_xlim(0, n_genomes)
    ax_strip.set_ylim(0, 1)
    ax_strip.axis("off")
    for arm in ARM_ORDER:
        idxs = [i for i, gid in enumerate(genome_order) if arm_row[gid] == arm]
        mid = (min(idxs) + max(idxs) + 1) / 2
        ax_strip.text(mid, 1.35, ARM_LABELS[arm], ha="center", va="bottom",
                      fontsize=9, color=INK_SECONDARY)

    # main heatmap
    disp = pivot_mm.copy()
    mask_ne = (pivot_cls == "NOT EVALUABLE")
    for i in range(n_assays):
        for j in range(n_genomes):
            val = disp.iat[i, j]
            if mask_ne.iat[i, j] or pd.isna(val):
                color = NOT_EVAL_GREY
            else:
                color = cmap(min(val, vmax) / vmax)
            ax.add_patch(plt.Rectangle((j, n_assays - i - 1), 1, 1, facecolor=color,
                                        edgecolor=SURFACE, linewidth=0.6))

    # arm separators
    boundaries = []
    prev = arm_row[genome_order[0]]
    for i, gid in enumerate(genome_order):
        if arm_row[gid] != prev:
            boundaries.append(i)
            prev = arm_row[gid]
    for b in boundaries:
        ax.axvline(b, color=INK_PRIMARY, linewidth=1.1)
        ax_strip.axvline(b, color=SURFACE, linewidth=1.1)

    ax.set_xlim(0, n_genomes)
    ax.set_ylim(0, n_assays)
    ax.set_xticks([])
    ax.set_yticks([n_assays - i - 0.5 for i in range(n_assays)])
    ax.set_yticklabels([ASSAY_LABELS[a] for a in ASSAY_ORDER], fontsize=9, color=INK_PRIMARY)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_facecolor(SURFACE)
    ax.set_xlabel(f"{n_genomes} BDBV genomes, grouped by outbreak", fontsize=9, color=INK_MUTED)

    # colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=vmax))
    cbar = fig.colorbar(sm, ax=ax, fraction=0.025, pad=0.015, aspect=30)
    cbar.set_label("Worst-component mismatches", fontsize=8, color=INK_SECONDARY)
    cbar.ax.tick_params(labelsize=8, color=INK_MUTED, labelcolor=INK_SECONDARY)

    # legend: not evaluable
    legend_handles = [Patch(facecolor=NOT_EVAL_GREY, label="NOT EVALUABLE (binding site unresolved or outside the record)")]
    ax.legend(handles=legend_handles, loc="upper center", bbox_to_anchor=(0.5, -0.12),
              frameon=False, fontsize=8, labelcolor=INK_SECONDARY, handlelength=1.2)

    fig.suptitle("In silico primer/probe mismatch coverage across five RT-PCR assays and 70 BDBV genomes",
                 fontsize=11, color=INK_PRIMARY, y=1.0)

    fig.savefig(OUT / "fig1_mismatch_heatmap.png", dpi=300, bbox_inches="tight", facecolor=SURFACE)
    fig.savefig(OUT / "fig1_mismatch_heatmap.pdf", bbox_inches="tight", facecolor=SURFACE)
    print(f"wrote {OUT}/fig1_mismatch_heatmap.png (+ .pdf)  [{n_assays} assays x {n_genomes} genomes, vmax={vmax}]")


if __name__ == "__main__":
    main()
