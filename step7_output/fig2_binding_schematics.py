#!/usr/bin/env python3
"""
Figure 2: binding-site schematics for the primers that matter most to the
coverage story - the real headline gap (Panning/Sanchez shared FiloB), a
clean broad-range performer for contrast (Greninger PfiloL_F), and one of
the intentionally off-target species-typing controls (Weidmann ENZ_F) shown
for contrast so a reader can see "designed to fail here" looks different
from "unexpectedly weak".

For each primer, recomputes the same mismatch-tolerant alignment Step 4
uses (best_alignment from step4_oligo_alignment.py) against all 70 curated
genomes, then plots one column per primer position showing what fraction
of genomes mismatch there. This is a per-position view Step 4/5/6 never
produced (they only ever kept whole-primer mismatch counts), so it is
computed fresh here, not read from an existing CSV.
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
import pandas as pd
from Bio import SeqIO

sys.path.insert(0, ".")
from step4_oligo_alignment import best_alignment, base_matches, rc

OUT = Path("figures")
OUT.mkdir(exist_ok=True)

SEQ_BLUE = ["#cde2fb", "#9ec5f4", "#5598e7", "#2a78d6", "#184f95", "#0d366b"]
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
SURFACE = "#fcfcfb"
SKIP_GREY = "#e1e0d9"

PANEL_FILE = Path("step2_output/validated_panel.csv")
GENOMES_FILE = Path("step1_output/curated_sequences.fasta")

PRIMERS = [
    ("PANNING2007_FiloB", "Panning/Sanchez shared FiloB (reverse) — the real REDUCED-tier gap"),
    ("GRENINGER2024_PfiloL_F", "Greninger 2024 PfiloL_F (forward) — clean broad-range performer"),
    ("WEIDMANN2004_ENZ_F", "Weidmann 2004 ENZ_F (forward) — intentional EBOV-Zaire-only control"),
]


def per_position_stats(oligo: str, genomes: dict) -> list[dict]:
    """For each genome, realign and return per-primer-position outcome
    ('match'/'mismatch'/'skip'), oriented 5'->3' as the primer is written."""
    olen = len(oligo)
    rows = []
    for gid, gseq in genomes.items():
        b = best_alignment(oligo, gseq)
        strand = b["strand"]
        probe = oligo if strand == "as_written" else rc(oligo)
        window = gseq[b["start"]:b["end"]]
        three_prime_at_start = strand == "reverse_complement"
        outcomes = [None] * olen
        for i, (ob, gb) in enumerate(zip(probe, window)):
            res = base_matches(ob, gb)
            pos = (olen - i - 1) if three_prime_at_start else i  # 0-indexed, 5'->3'
            outcomes[pos] = res
        rows.append({"genome_id": gid, "outcomes": outcomes})
    return rows


def main():
    panel = pd.read_csv(PANEL_FILE, dtype=str)
    panel.columns = panel.columns.str.strip()
    genomes = {rec.id.split(".")[0]: str(rec.seq).upper()
               for rec in SeqIO.parse(str(GENOMES_FILE), "fasta")}

    cmap = LinearSegmentedColormap.from_list("seq_blue", SEQ_BLUE, N=256)

    fig, axes = plt.subplots(len(PRIMERS), 1, figsize=(9, 2.1 * len(PRIMERS) + 0.6),
                              gridspec_kw={"hspace": 1.5})
    fig.patch.set_facecolor(SURFACE)

    for ax, (assay_id, title) in zip(axes, PRIMERS):
        row = panel[panel.assay_id == assay_id].iloc[0]
        oligo = str(row.sequence_5to3).upper().strip()
        olen = len(oligo)
        rows = per_position_stats(oligo, genomes)
        n = len(rows)

        frac_mismatch = []
        frac_skip = []
        for pos in range(olen):
            mm = sum(1 for r in rows if r["outcomes"][pos] == "mismatch")
            sk = sum(1 for r in rows if r["outcomes"][pos] == "skip")
            frac_mismatch.append(mm / n)
            frac_skip.append(sk / n)

        for pos in range(olen):
            if frac_skip[pos] > 0.5:
                color = SKIP_GREY
            else:
                color = cmap(frac_mismatch[pos])
            ax.add_patch(Rectangle((pos, 0), 1, 1, facecolor=color, edgecolor=SURFACE, linewidth=1))
            ax.text(pos + 0.5, 0.5, oligo[pos], ha="center", va="center",
                     fontsize=9, color=INK_PRIMARY if frac_mismatch[pos] < 0.6 else "white",
                     family="monospace")

        ax.set_xlim(0, olen)
        ax.set_ylim(-0.65, 1.55)
        ax.set_yticks([])
        ax.set_xticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_facecolor(SURFACE)
        ax.text(0, 1.12, "5'", fontsize=8, color=INK_MUTED, ha="left", va="bottom")
        ax.text(olen, 1.12, "3'", fontsize=8, color=INK_MUTED, ha="right", va="bottom")
        # mark the 3' terminal window (last 5 nt) - the biologically critical zone
        ax.add_patch(Rectangle((olen - 5, -0.22), 5, 0.1, facecolor=INK_PRIMARY, edgecolor="none"))
        ax.text(olen - 2.5, -0.55, "3' window (last 5 nt)", fontsize=7, color=INK_MUTED, ha="center")
        ax.set_title(title, fontsize=9.5, color=INK_PRIMARY, loc="left", pad=20)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=1))
    cbar = fig.colorbar(sm, ax=axes, fraction=0.03, pad=0.02, aspect=25, location="right")
    cbar.set_label("Fraction of 70 genomes mismatched at this position", fontsize=8, color=INK_SECONDARY)
    cbar.ax.tick_params(labelsize=8, labelcolor=INK_SECONDARY)

    fig.suptitle("Per-position mismatch frequency across 70 BDBV genomes, three representative primers",
                 fontsize=11, color=INK_PRIMARY, x=0.42, y=1.0)

    fig.savefig(OUT / "fig2_binding_site_schematics.png", dpi=300, bbox_inches="tight", facecolor=SURFACE)
    fig.savefig(OUT / "fig2_binding_site_schematics.pdf", bbox_inches="tight", facecolor=SURFACE)
    print(f"wrote {OUT}/fig2_binding_site_schematics.png (+ .pdf)")


if __name__ == "__main__":
    main()
