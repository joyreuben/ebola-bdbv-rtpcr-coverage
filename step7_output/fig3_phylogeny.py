#!/usr/bin/env python3
"""
Figure 3: phylogeny of the 70 curated BDBV genomes (MAFFT alignment ->
FastTree GTR, both run once as part of this figure's pipeline - see
step7_output/mafft.log and fasttree.log), annotated with two tip strips:

  1. outbreak arm (2007 / 2012 / 2026) - categorical
  2. Panning/Sanchez FiloB worst-component mismatch count - the project's
     headline REDUCED-tier result.

WHY THIS ANNOTATION CHANGED (2026-08-28)
    This figure previously annotated tips with the Ogawa 2010 mismatch
    count, on the stated grounds that Ogawa was "the one assay with real
    per-genome variation (0-11 mm)" while everything else was flat. That
    variation was not real. It was a single record, MS811591 (PP_006X6B5,
    a 6,633 bp partial genome), whose Ogawa binding site is not contained
    in the record at all - Step 4 had scored the primers against unrelated
    sequence and returned a phantom 11 mm. With Step 4b's coverage check in
    place, that pair is NOT EVALUABLE and Ogawa is a uniform 2 mm across
    every evaluable genome in every arm.

    So no assay in this panel has per-genome texture. That is not a reason
    to drop the figure - it IS the finding, and it is stronger stated
    positively: the tree carries clear phylogenetic structure (three
    outbreak-pure clades, sampled ~19 years apart), and against that
    structure the assay mismatch count is invariant. A flat strip beside a
    structured tree is exactly the right picture for "species-level
    divergence, not lineage-specific drift". The annotation is therefore
    switched to the primer the paper actually reports on.

The question this figure answers visually: does the FiloB mismatch level
track clade structure or outbreak arm (which would suggest a lineage-
specific effect), or is it invariant across a tree that is itself clearly
structured (consistent with the "species-level gap, not a newly evolving
mutation" reading)?

Genome IDs are not shown as tip labels (70 leaves, unreadable) - only the
two colour strips plus a legend. Reads step1_output/curated_metadata.csv
for outbreak arm and step5_output/assay_coverage_per_genome.csv for the
Panning/Sanchez FiloB mismatch call.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle, Patch
import pandas as pd
from Bio import Phylo

OUT = Path("figures")
OUT.mkdir(exist_ok=True)

SEQ_BLUE = ["#cde2fb", "#9ec5f4", "#5598e7", "#2a78d6", "#184f95", "#0d366b"]
NOT_EVAL_GREY = "#c3c2b7"
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
SURFACE = "#fcfcfb"
ARM_COLORS = {"Bdbv-2007": "#eb6834", "Bdbv-2012": "#1baf7a", "Bdbv-2026": "#eda100"}
ARM_LABELS = {"Bdbv-2007": "2007", "Bdbv-2012": "2012", "Bdbv-2026": "2026"}


def main():
    tree = Phylo.read("step7_output/tree.nwk", "newick")
    tree.root_at_midpoint()
    tree.ladderize()

    meta = pd.read_csv("step1_output/curated_metadata.csv", keep_default_na=False)
    meta["genome_id"] = meta["accessionVersion"].str.split(".").str[0]
    arm_by_id = meta.set_index("genome_id")["outbreak"].to_dict()

    cov = pd.read_csv("step5_output/assay_coverage_per_genome.csv")
    assay = cov[cov.assay == "Panning2007_panfilo"].set_index("genome_id")
    mm_by_id = assay["worst_mm"].to_dict()
    cls_by_id = assay["classification"].to_dict()

    terminals = tree.get_terminals()
    n = len(terminals)
    # Scale to the classification thresholds, not to the observed range. The
    # observed range is a single value (the finding is that it is invariant),
    # so a data-driven vmax would map that constant onto the top of the ramp
    # and imply a spread that does not exist. Fixing the scale at the
    # LIKELY FAILURE boundary keeps the colour honest and comparable to Fig 1.
    vmax = 6

    cmap = LinearSegmentedColormap.from_list("seq_blue", SEQ_BLUE, N=256)

    fig_h = max(8, n * 0.14)
    fig, (ax_tree, ax_strip) = plt.subplots(
        1, 2, figsize=(11.5, fig_h),
        gridspec_kw={"width_ratios": [5, 0.4], "right": 0.82}, sharey=False,
    )
    fig.patch.set_facecolor(SURFACE)

    Phylo.draw(tree, do_show=False, axes=ax_tree, label_func=lambda x: "",
               show_confidence=False)
    for coll in ax_tree.collections + ax_tree.lines:
        coll.set_color(INK_SECONDARY)
        coll.set_linewidth(0.9)
    ax_tree.set_facecolor(SURFACE)
    ax_tree.set_title("Maximum-likelihood tree (FastTree, GTR, midpoint-rooted)\nfull-genome MAFFT alignment, n=70",
                       fontsize=10, color=INK_PRIMARY, loc="left")
    ax_tree.set_xlabel("substitutions / site", fontsize=8, color=INK_MUTED)
    ax_tree.set_ylabel("")
    for spine in ("top", "right", "left"):
        ax_tree.spines[spine].set_visible(False)
    ax_tree.tick_params(axis="y", length=0, labelleft=False)
    ax_tree.tick_params(axis="x", labelsize=7, colors=INK_MUTED)

    # tip y-positions come straight from Bio.Phylo.draw's own layout;
    # copy ax_tree's actual y-limits (rather than guessing the inversion)
    # so the strip lines up with the tree exactly as drawn.
    y_positions = {t: i + 1 for i, t in enumerate(terminals)}

    ax_strip.set_xlim(0, 2)
    ax_strip.set_ylim(ax_tree.get_ylim())
    for term, y in y_positions.items():
        gid = term.name.split(".")[0]
        arm = arm_by_id.get(gid)
        mm = mm_by_id.get(gid)
        cls = cls_by_id.get(gid)
        arm_color = ARM_COLORS.get(arm, "white")
        if cls == "NOT EVALUABLE" or mm is None:
            mm_color = NOT_EVAL_GREY
        else:
            mm_color = cmap(mm / vmax)
        ax_strip.add_patch(Rectangle((0, y - 0.5), 1, 1, facecolor=arm_color, edgecolor=SURFACE, linewidth=0.3))
        ax_strip.add_patch(Rectangle((1, y - 0.5), 1, 1, facecolor=mm_color, edgecolor=SURFACE, linewidth=0.3))

    ax_strip.set_xticks([0.5, 1.5])
    ax_strip.set_xticklabels(["arm", "FiloB\nmm"], fontsize=6.5, color=INK_SECONDARY, rotation=90)
    ax_strip.set_yticks([])
    ax_strip.tick_params(length=0, pad=2)
    for spine in ax_strip.spines.values():
        spine.set_visible(False)
    ax_strip.set_facecolor(SURFACE)

    # legend and colorbar are placed in figure-fraction coordinates (not
    # relative to ax_strip's own bbox) so they sit clearly clear of the
    # strip itself - anchoring to the axes here previously overlapped the
    # legend swatches with the real data cells.
    arm_handles = [Patch(facecolor=c, label=ARM_LABELS[a]) for a, c in ARM_COLORS.items()]
    ne_handle = [Patch(facecolor=NOT_EVAL_GREY, label="NOT EVALUABLE")]
    fig.legend(handles=arm_handles + ne_handle, loc="upper left",
               bbox_to_anchor=(0.855, 0.88), frameon=False, fontsize=8,
               labelcolor=INK_SECONDARY, handlelength=1.1, title="outbreak arm",
               title_fontsize=8.5)

    cbar_ax = fig.add_axes([0.90, 0.30, 0.02, 0.35])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=vmax))
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_label("Panning/Sanchez FiloB\nmismatches", fontsize=8, color=INK_SECONDARY)
    cbar.ax.tick_params(labelsize=8, labelcolor=INK_SECONDARY)

    fig.suptitle("BDBV genomes form three outbreak-pure clades, but assay mismatch is invariant across all of them",
                 fontsize=10.5, color=INK_PRIMARY, y=0.995)

    fig.savefig(OUT / "fig3_phylogeny.png", dpi=300, bbox_inches="tight", facecolor=SURFACE)
    fig.savefig(OUT / "fig3_phylogeny.pdf", bbox_inches="tight", facecolor=SURFACE)
    print(f"wrote {OUT}/fig3_phylogeny.png (+ .pdf)  [{n} tips]")


if __name__ == "__main__":
    main()
