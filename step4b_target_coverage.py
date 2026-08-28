#!/usr/bin/env python3
"""
Step 4b: Per-genome target-site COVERAGE detection.

WHY THIS EXISTS
    Step 4's aligner always returns a best-scoring window. It has no concept
    of "this oligo's target region is not in this sequence at all". Handed a
    partial genome that simply does not span an oligo's binding site, it
    slides the oligo across whatever sequence IS present and reports the
    least-bad window it finds. The resulting mismatch count is an artifact of
    an unrelated region, not a measurement of the real binding site.

    This is not hypothetical. MS811591.1 (record PP_006X6B5, 2007 arm) is a
    6,633 bp partial genome covering only NC_014373.1 positions 11,567-18,199.
    The Ogawa 2010 primers target positions 1,028-1,621 - roughly 10 kb
    outside that window, i.e. absent from the record entirely. Step 4 reported
    11 mismatches for them, and Step 5 classified the assay LIKELY FAILURE for
    that genome. Both numbers were measuring unrelated L-gene sequence. The
    manuscript then described this as a "genuine biological outlier". It was
    not; it was a coverage artifact.

    Note this is the exact confound already flagged in the manuscript's
    discussion of the 8 EXCLUDED records ("an elevated mismatch is also
    consistent with the true binding region simply not being present in the
    retrieved excerpt"). That caveat was written about excluded records and
    never checked against the retained 70. Running it finds one hit.

WHAT IT DOES
    1. Locates each oligo's target site on NC_014373.1 (the BDBV RefSeq
       entry) via Step 4's own best_alignment(), so coordinates agree with
       Step 4 / the MIQE table by construction.
    2. Determines, for every curated genome, which NC_014373.1 coordinate
       ranges that genome actually covers. Method: exact-match short anchors
       sampled along the genome into the reference, then merge the hit
       intervals. Exact matching is safe here because BDBV genomes are
       ~99% identical to the BDBV reference, so anchors land densely; small
       gaps (divergent spots where one anchor failed) are bridged.
    3. Emits, per (oligo, genome), whether that oligo's reference target
       site falls inside a covered range.

    Step 5 consumes this and forces NOT EVALUABLE when a required component's
    target site is absent, instead of trusting a phantom mismatch count.

WHY EXACT-MATCH ANCHORS RATHER THAN FULL PAIRWISE ALIGNMENT
    A 19 kb x 19 kb pairwise alignment per genome is slow and buys nothing
    here: we need coverage extent, not per-base homology. Anchors give the
    same answer in well under a second per genome, degrade gracefully on
    divergent records (more failed anchors, bridged by the gap tolerance),
    and are trivially auditable - every covered interval traces to a
    specific exact sequence match.

Requires: pandas, biopython
Inputs : step2_output/validated_panel.csv
         step1_output/curated_sequences.fasta
         references/BDBV_NC_014373.1.fasta
Outputs: step4_output/target_presence.csv
         step4_output/genome_reference_coverage.csv
         step4_output/step4b_report.txt
"""

import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO

import step4_oligo_alignment as s4

PANEL_FILE = Path("step2_output/validated_panel.csv")
GENOMES_FILE = Path("step1_output/curated_sequences.fasta")
REF_FILE = Path("references/BDBV_NC_014373.1.fasta")
OUT_DIR = Path("step4_output")

# --- anchor parameters ------------------------------------------------
# ANCHOR_LEN: long enough to be unique in a 19 kb genome (4^40 >> 19000),
#   short enough to survive local divergence.
# ANCHOR_STRIDE: sample density along the genome.
# BRIDGE_GAP: merge covered intervals separated by less than this. A gap
#   this small means a handful of consecutive anchors failed on local
#   divergence, not that the genome is missing a real chunk. Set well below
#   the ~10 kb scale of a genuinely absent region.
ANCHOR_LEN = 40
ANCHOR_STRIDE = 100
BRIDGE_GAP = 1000


def covered_intervals(genome: str, reference: str) -> list:
    """
    Return merged [start, end) reference intervals (0-based) that this
    genome demonstrably covers, via exact anchor matching.
    """
    hits = []
    for i in range(0, len(genome) - ANCHOR_LEN + 1, ANCHOR_STRIDE):
        anchor = genome[i:i + ANCHOR_LEN]
        if "N" in anchor:
            continue
        pos = reference.find(anchor)
        if pos == -1:
            continue
        # Reject repeat-derived anchors: a second occurrence means the
        # mapping is ambiguous and the interval is not evidence of coverage.
        if reference.find(anchor, pos + 1) != -1:
            continue
        hits.append((pos, pos + ANCHOR_LEN))

    if not hits:
        return []

    hits.sort()
    merged = [list(hits[0])]
    for start, end in hits[1:]:
        if start - merged[-1][1] <= BRIDGE_GAP:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [tuple(iv) for iv in merged]


def site_is_covered(start: int, end: int, intervals: list) -> bool:
    """
    True if reference range [start, end) lies within the genome's overall
    mapped span.

    DELIBERATELY SPAN-BASED, NOT INTERVAL-BASED. An earlier version required
    the site to sit inside a single merged interval, which also fired on
    INTERNAL gaps. Those gaps are not missing sequence - they are long N-runs
    (or local divergence) where anchors failed, in genomes that do span the
    coordinate. Verified: the internal gaps this flagged in PP_006X6QC /
    PP_006X6JQ are N-runs at 11829-13044, and Step 5 already routes those to
    NOT EVALUABLE via skip_fraction (0.91-0.96, far above the 0.15 cutoff).
    Treating them as "absent" here would be redundant at best; treating a
    merely-divergent internal region as absent would be actively wrong,
    discarding a real measurement.

    The span test isolates the one thing skip_fraction genuinely cannot see:
    a record that does not reach the coordinate at all, where Step 4 has
    silently scored some unrelated region instead.
    """
    if not intervals:
        return False
    span_start = min(iv_s for iv_s, _ in intervals)
    span_end = max(iv_e for _, iv_e in intervals)
    return span_start <= start and end <= span_end


def main() -> None:
    for f in (PANEL_FILE, GENOMES_FILE, REF_FILE):
        if not f.exists():
            sys.exit(f"Missing {f}")
    OUT_DIR.mkdir(exist_ok=True)

    report = []

    def log(line=""):
        print(line)
        report.append(line)

    log("=" * 70)
    log("STEP 4b: TARGET-SITE COVERAGE DETECTION")
    log("=" * 70)
    log()

    reference = str(next(SeqIO.parse(str(REF_FILE), "fasta")).seq).upper()
    panel = pd.read_csv(PANEL_FILE, dtype=str)
    panel.columns = panel.columns.str.strip()
    genomes = {rec.id.split(".")[0]: str(rec.seq).upper()
               for rec in SeqIO.parse(str(GENOMES_FILE), "fasta")}

    log(f"Reference : NC_014373.1 ({len(reference)} bp)")
    log(f"Oligos    : {len(panel)}")
    log(f"Genomes   : {len(genomes)}")
    log()

    # --- 1. oligo target sites on the reference -----------------------
    log("-" * 70)
    log("Oligo target sites on NC_014373.1")
    log("-" * 70)
    sites = {}
    for _, row in panel.iterrows():
        oligo = str(row["sequence_5to3"]).upper().strip()
        b = s4.best_alignment(oligo, reference)
        start = b["start"]
        sites[row["assay_id"]] = (start, start + len(oligo))
        log(f"  {row['assay_id']:<26} {start + 1:>6}-{start + len(oligo):<6} "
            f"({b['strand']}, mm={b['mismatches']})")
    log()

    # --- 2. per-genome reference coverage ------------------------------
    log("-" * 70)
    log("Per-genome reference coverage")
    log("-" * 70)
    coverage, cov_rows = {}, []
    for gid, gseq in genomes.items():
        intervals = covered_intervals(gseq, reference)
        coverage[gid] = intervals
        span = sum(e - s for s, e in intervals)
        cov_rows.append({
            "genome_id": gid,
            "genome_length": len(gseq),
            "reference_bp_covered": span,
            "fraction_of_reference": round(span / len(reference), 4),
            "n_intervals": len(intervals),
            "covered_ranges_1based": "; ".join(
                f"{s + 1}-{e}" for s, e in intervals) or "none",
        })

    cov_df = pd.DataFrame(cov_rows).sort_values("reference_bp_covered")
    cov_df.to_csv(OUT_DIR / "genome_reference_coverage.csv", index=False)

    partial = cov_df[cov_df["fraction_of_reference"] < 0.95]
    if partial.empty:
        log("  All genomes cover >=95% of the reference.")
    else:
        log(f"  {len(partial)} genome(s) cover <95% of the reference:")
        for _, r in partial.iterrows():
            log(f"    {r['genome_id']:<14} {r['genome_length']:>6} bp  "
                f"covers {r['fraction_of_reference']:.1%}  "
                f"[{r['covered_ranges_1based']}]")
    log()

    # --- 3. per (oligo, genome) target presence ------------------------
    rows = []
    for assay_id, (s, e) in sites.items():
        for gid in genomes:
            rows.append({
                "assay_id": assay_id,
                "genome_id": gid,
                "target_start_1based": s + 1,
                "target_end_1based": e,
                "target_present": site_is_covered(s, e, coverage[gid]),
            })
    pres = pd.DataFrame(rows)
    pres.to_csv(OUT_DIR / "target_presence.csv", index=False)

    absent = pres[~pres["target_present"]]
    log("-" * 70)
    log("ABSENT TARGET SITES  (Step 4 mismatch counts here are artifacts)")
    log("-" * 70)
    if absent.empty:
        log("  None. Every oligo target site is present in every genome.")
    else:
        log(f"  {len(absent)} of {len(pres)} oligo/genome pairs "
            f"({len(absent) / len(pres):.1%}) have an absent target site,")
        log(f"  across {absent['genome_id'].nunique()} genome(s):")
        for gid, grp in absent.groupby("genome_id"):
            log(f"    {gid}: {len(grp)} oligo(s) - "
                f"{', '.join(sorted(grp['assay_id']))}")
    log()
    log("=" * 70)
    log(f"WROTE  {OUT_DIR}/target_presence.csv")
    log(f"WROTE  {OUT_DIR}/genome_reference_coverage.csv")
    log("=" * 70)
    log()
    log("NEXT: Step 5 reads target_presence.csv and classifies an assay")
    log("NOT EVALUABLE for a genome when a required component's target site")
    log("is absent, rather than trusting a phantom mismatch count.")

    (OUT_DIR / "step4b_report.txt").write_text("\n".join(report))


if __name__ == "__main__":
    main()
