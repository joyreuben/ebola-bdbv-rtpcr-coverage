#!/usr/bin/env python3
"""
Step 6: primer3-py thermodynamics on the REDUCED calls.

WHY THIS EXISTS
    Step 5's classification is mismatch-count based (a proxy). Step 6 asks
    the more direct question for the borderline REDUCED calls: given the
    actual mismatches at the actual binding site, how much does the
    predicted duplex melting temperature (Tm) drop relative to a perfect
    match? A large Tm drop is what actually causes reduced amplification in
    the wet lab, so this is the confirmation step before writing any claim
    about assay performance for these borderline genomes.

SCOPE - REDUCED ONLY, NOT LIKELY FAILURE
    Tm modeling is run only on REDUCED-tier rows (3-5 mismatches, the
    borderline calls where thermodynamic confirmation actually adds
    information). LIKELY FAILURE rows (6+ mismatches, or a 3'-end
    mismatch) are reported by mismatch count alone, NOT modeled with
    primer3, for a methodological reason discovered while building this
    script (see "WHY LIKELY FAILURE IS NOT TM-MODELED" below) - this was
    not a scope simplification, it was a correctness finding.
    (AT RISK is empty post phantom-10 fix; LIKELY DETECTED and
    NOT EVALUABLE need no thermodynamic confirmation either way.)

METHOD (REDUCED rows)
    For each chosen primer/probe component of that row (forward, reverse,
    probe - whichever the assay actually uses), reconstruct the real primer-
    template duplex from Step 4's alignment coordinates (best_start, strand,
    oligo_length) and Step 1's genome sequences, then compare:
      - ref_tm    = calc_heterodimer(primer, rc(primer)).tm   perfect-match Tm
      - duplex_tm = calc_heterodimer(primer, template).tm     actual, mismatch-aware Tm
      - tm_shift  = ref_tm - duplex_tm

    Orientation of the template strand (see build_template() docstring)
    matters for calc_heterodimer, since it is not a simple identity check -
    it is a genuine antiparallel nearest-neighbor duplex calculation.

WHY LIKELY FAILURE IS NOT TM-MODELED
    primer3's calc_heterodimer() does not score the fixed, full-length
    duplex at Step 4's known binding site - it SEARCHES for whichever local
    partial alignment is thermodynamically most favorable, allowing large
    unpaired dangling ends. With few mismatches (REDUCED tier) the full
    intended site remains that optimum, so the numbers are sound. With many
    mismatches (LIKELY FAILURE tier, verified up to 8 mm here), primer3
    instead finds a short, mostly-unrelated stable sub-stretch and reports
    ITS Tm - producing physically nonsensical results (e.g. a duplex Tm of
    -25C, an "89C shift", for one WEIDMANN2004_ENZ_P/genome pair checked by
    hand). primer3 has no mode to force the known, fixed-frame alignment,
    so for these rows a Tm number would misrepresent, not confirm, the
    mismatch-based call. Mismatch count/3'-position is already decisive at
    this tier regardless (that is exactly why Step 5 calls it LIKELY
    FAILURE), so nothing is lost by not modeling it.

CAVEATS ON THE REDUCED-TIER NUMBERS (see step6_report.txt too - keep this
in the writeup)
    - Several primers/probes use real IUPAC degenerate bases (Nakayama,
      Greninger, the Panning/Weidmann probes), not just inosine-N. primer3's
      nearest-neighbor model only accepts ACGT, so every degenerate position
      is first "realized" against the actual genome base using Step 4's own
      match/mismatch/skip rule (see realize_duplex()): match/skip positions
      become a non-penalizing perfect pair, mismatch positions become a
      genuine concrete mismatch. This keeps the mismatch count driving each
      Tm number identical to Step 4/5's own count - flagged per-row via
      primer_has_degenerate_base for transparency, not because it's a gap.
    - A handful of genome windows still contain a real (non-inosine) N from
      sequencing ambiguity even though they passed the <=15% skip-fraction
      gate at Step 5. That position is realized as non-penalizing (same
      "skip, not scoring" convention as Step 4) rather than invented as a
      mismatch. Flagged per-row via window_has_N.
    - Tm shift is a predictor of risk, not proof of assay failure - same
      caveat as Step 5. State this plainly when reporting results.

RUN:
    python step6_thermodynamics.py

Requires: pandas, biopython, primer3-py
    pip install pandas biopython primer3-py --break-system-packages
"""

import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO
import primer3

import step4_oligo_alignment as s4  # reuse the exact match/mismatch/skip rules

PANEL_FILE = Path("assay_panel.csv")
GENOMES_FILE = Path("step1_output/curated_sequences.fasta")
STEP4_FILE = Path("step4_output/oligo_vs_genome_full.csv")
STEP5_FILE = Path("step5_output/assay_coverage_per_genome.csv")
OUT_DIR = Path("step6_output")

COMPONENT_ROLES = ["forward", "reverse", "probe"]

# Heuristic Tm-shift bands, consistent in spirit with the Step 5 thresholds
# (declared here, not derived from data). Not a claim of assay failure by
# itself - see caveats above.
CONFIRM_STRONG = 5.0   # >= this many degC drop: thermodynamics confirms risk
CONFIRM_MARGINAL = 2.0  # >= this: marginal support; below: thermo does not add support


def realize_duplex(oligo: str, genome_window: str, strand: str) -> dict:
    """
    Turn a (possibly IUPAC-degenerate) primer/probe and its genome window
    into a concrete all-ACGT primer + template pair that primer3 can score,
    using the SAME per-position match/mismatch/skip calls as Step 4
    (`s4.base_matches`) so the mismatch count here always agrees with
    Step 4/5's own numbers. Concretely:

      - "match" or "skip" positions (incl. inosine-derived N on either
        side) are realized as a perfect Watson-Crick pair - non-penalizing,
        consistent with Step 2/4's declared "inosine/N is non-scoring" rule.
      - "mismatch" positions are realized as a genuine concrete mismatch:
        the genome's actual base opposite a fixed, different member of the
        primer symbol's IUPAC-allowed set.

    Returns the concrete primer (5'->3') and the concrete template strand
    (5'->3') it actually anneals to - see build_template()'s docstring
    below for the strand-orientation logic, applied here to the realized
    (all-ACGT) genome window rather than the raw one.
    """
    probe = oligo if strand == "as_written" else s4.rc(oligo)
    concrete_probe, concrete_genome = [], []
    window_has_N = "N" in genome_window.upper()

    for ob, gb in zip(probe, genome_window):
        ob_u, gb_u = ob.upper(), gb.upper()
        result = s4.base_matches(ob_u, gb_u)
        if result == "mismatch":
            concrete_genome.append(gb_u)  # guaranteed ACGT (base_matches requires it)
            allowed = s4.IUPAC.get(ob_u, {ob_u})
            concrete_probe.append(sorted(allowed)[0])  # any member differs from gb_u by definition
        else:  # match or skip -> realize as a perfect pair, non-penalizing
            if gb_u in "ACGT":
                realized = gb_u
            elif ob_u in "ACGT":
                realized = ob_u
            else:
                realized = "A"  # both sides unresolved (rare); arbitrary, non-penalizing on both sides
            concrete_probe.append(realized)
            concrete_genome.append(realized)

    concrete_probe_str = "".join(concrete_probe)
    concrete_genome_str = "".join(concrete_genome)
    concrete_primer = concrete_probe_str if strand == "as_written" else s4.rc(concrete_probe_str)
    return {"primer": concrete_primer, "genome": concrete_genome_str,
            "window_has_N": window_has_N}


def build_template(concrete_genome_window: str, strand: str) -> str:
    """
    Reconstruct the actual template strand (5'->3') the primer physically
    anneals to, from the (now fully-realized, all-ACGT) genome window.

    Step 4's aligner compares the oligo (or its reverse complement, per
    `strand`) against the genome's stored (sense/top) strand base-for-base.
    That is an identity check, not a duplex check - it tells us where the
    binding site is, not what strand a primer actually hybridizes to.

    - strand == "as_written": the primer sequence matched the top strand
      directly, meaning the primer binds the strand complementary to it,
      i.e. template = reverse_complement(genome_window).
    - strand == "reverse_complement": rc(primer) matched the top strand,
      meaning the primer binds the top strand directly as stored,
      i.e. template = genome_window (unchanged).

    In both cases, when there are zero mismatches this reduces to
    template == rc(primer), the exact complement - the correct physical
    duplex. Verified against primer3 directly before use in this script.
    """
    if strand == "as_written":
        return s4.rc(concrete_genome_window)
    if strand == "reverse_complement":
        return concrete_genome_window
    raise ValueError(f"unknown strand value: {strand}")


def tm_shift_for(primer_seq: str, genome_window: str, strand: str) -> dict:
    realized = realize_duplex(primer_seq, genome_window, strand)
    concrete_primer = realized["primer"]
    template = build_template(realized["genome"], strand)
    ref = primer3.bindings.calc_heterodimer(concrete_primer, s4.rc(concrete_primer))
    duplex = primer3.bindings.calc_heterodimer(concrete_primer, template)
    return {"ref_tm": round(ref.tm, 2), "duplex_tm": round(duplex.tm, 2),
             "tm_shift": round(ref.tm - duplex.tm, 2),
             "window_has_N": realized["window_has_N"]}


def confirm_flag(tm_shift):
    if tm_shift >= CONFIRM_STRONG:
        return "THERMO CONFIRMS"
    if tm_shift >= CONFIRM_MARGINAL:
        return "THERMO MARGINAL"
    return "THERMO MINIMAL"


def main() -> None:
    OUT_DIR.mkdir(exist_ok=True)
    report = []

    def log(line=""):
        print(line)
        report.append(line)

    log("=" * 70)
    log("STEP 6: PRIMER3-PY THERMODYNAMICS ON REDUCED / LIKELY FAILURE CALLS")
    log("=" * 70)
    log()

    for f in (PANEL_FILE, GENOMES_FILE, STEP4_FILE, STEP5_FILE):
        if not f.exists():
            sys.exit(f"Required input not found: {f}")

    panel = pd.read_csv(PANEL_FILE, dtype=str)
    panel.columns = panel.columns.str.strip()
    primer_seq_by_id = dict(zip(panel["assay_id"], panel["sequence_5to3"].str.upper()))

    genomes = {rec.id.split(".")[0]: str(rec.seq).upper()
               for rec in SeqIO.parse(str(GENOMES_FILE), "fasta")}

    step4 = pd.read_csv(STEP4_FILE, dtype={"genome_id": str})
    step4_idx = step4.set_index(["assay_id", "genome_id"])

    step5 = pd.read_csv(STEP5_FILE, dtype={"genome_id": str})
    reduced = step5[step5["classification"] == "REDUCED"].copy()
    failure = step5[step5["classification"] == "LIKELY FAILURE"].copy()
    log(f"REDUCED rows (Tm-modeled): {len(reduced)} of {len(step5)}")
    log(f"LIKELY FAILURE rows (mismatch-count only, see module docstring")
    log(f" 'WHY LIKELY FAILURE IS NOT TM-MODELED'): {len(failure)} of {len(step5)}")
    log("(AT RISK is empty post phantom-10 fix; LIKELY DETECTED / NOT EVALUABLE")
    log(" need no confirmation either way.)")
    log()

    rows = []
    for _, row in reduced.iterrows():
        genome_seq = genomes.get(row["genome_id"])
        if genome_seq is None:
            log(f"  ! genome {row['genome_id']} not found in {GENOMES_FILE}, skipped")
            continue

        for role in COMPONENT_ROLES:
            oligo_id = row.get(f"chosen_{role}")
            if pd.isna(oligo_id):
                continue
            primer_seq = primer_seq_by_id.get(oligo_id)
            if primer_seq is None:
                log(f"  ! oligo {oligo_id} not found in {PANEL_FILE}, skipped")
                continue
            key = (oligo_id, row["genome_id"])
            if key not in step4_idx.index:
                log(f"  ! no Step 4 alignment for {key}, skipped")
                continue
            a = step4_idx.loc[key]
            window = genome_seq[a["best_start"]: a["best_start"] + a["oligo_length"]]

            t = tm_shift_for(primer_seq, window, a["strand"])
            rows.append({
                "assay": row["assay"],
                "genome_id": row["genome_id"],
                "outbreak": row["outbreak"],
                "classification": row["classification"],
                "is_worst_component": role == row["worst_component"],
                "component_role": role,
                "oligo_id": oligo_id,
                "primer_has_degenerate_base": any(b not in "ACGT" for b in primer_seq),
                "mismatches": int(a["mismatches"]),
                "three_prime_mismatches": int(a["three_prime_mismatches"]),
                "window_has_N": t["window_has_N"],
                "ref_tm": t["ref_tm"],
                "duplex_tm": t["duplex_tm"],
                "tm_shift": t["tm_shift"],
                "confirm_flag": confirm_flag(t["tm_shift"]),
            })

    out = pd.DataFrame(rows)
    out.to_csv(OUT_DIR / "thermodynamics_per_component.csv", index=False)
    log(f"Wrote {len(out)} component-level rows to thermodynamics_per_component.csv")
    log("(REDUCED tier only - see module docstring for why LIKELY FAILURE is excluded)")
    log()

    # --- assay-level summary, driven by the worst_component row only ------
    worst_only = out[out["is_worst_component"]].copy()
    summary = (worst_only.groupby(["assay", "outbreak", "classification"])
               .agg(n=("tm_shift", "size"),
                    mean_tm_shift=("tm_shift", "mean"),
                    min_tm_shift=("tm_shift", "min"),
                    max_tm_shift=("tm_shift", "max"))
               .round(2).reset_index())
    summary.to_csv(OUT_DIR / "thermodynamics_summary.csv", index=False)

    log("-" * 70)
    log("SUMMARY  Tm shift on the classification-driving (worst) component")
    log("REDUCED tier only")
    log("-" * 70)
    log(summary.to_string(index=False))
    log()

    confirm_counts = (worst_only.groupby(["assay", "outbreak"])["confirm_flag"]
                       .value_counts().unstack(fill_value=0))
    log("-" * 70)
    log("CONFIRM-FLAG COUNTS  (>=5C THERMO CONFIRMS, 2-5C MARGINAL, <2C MINIMAL)")
    log("-" * 70)
    log(confirm_counts.to_string())
    log()

    # --- LIKELY FAILURE: mismatch-count summary only, no Tm claim ---------
    fail_worst = failure[["assay", "outbreak", "worst_mm", "worst_3p_mm"]].copy()
    fail_summary = (fail_worst.groupby(["assay", "outbreak"])
                     .agg(n=("worst_mm", "size"),
                          mean_worst_mm=("worst_mm", "mean"),
                          min_worst_mm=("worst_mm", "min"),
                          max_worst_mm=("worst_mm", "max"))
                     .round(2).reset_index())
    fail_summary.to_csv(OUT_DIR / "likely_failure_mismatch_summary.csv", index=False)
    log("-" * 70)
    log("LIKELY FAILURE  mismatch-count summary (NOT Tm-modeled, see docstring)")
    log("-" * 70)
    log(fail_summary.to_string(index=False))
    log()

    n_degen = worst_only["primer_has_degenerate_base"].sum()
    n_window_n = worst_only["window_has_N"].sum()
    log("-" * 70)
    log("CAVEATS")
    log("-" * 70)
    log("  LIKELY FAILURE rows are NOT Tm-modeled: primer3's calc_heterodimer")
    log("  searches for the best local partial alignment rather than scoring")
    log("  Step 4's known fixed-frame binding site, and at high mismatch")
    log("  counts it finds a short unrelated sub-stretch instead, producing")
    log("  physically nonsensical Tm values (verified: -25C duplex Tm, an")
    log("  '89C shift', for one WEIDMANN2004_ENZ_P case). Mismatch count/")
    log("  3'-position is already decisive at this tier - that is why Step 5")
    log("  calls it LIKELY FAILURE in the first place.")
    log(f"  {n_degen} REDUCED worst-component rows use a primer/probe with an")
    log("  IUPAC degenerate base (incl. inosine-derived N). Those positions")
    log("  were realized against the actual genome base using Step 4's own")
    log("  match/mismatch/skip rule before Tm modeling, so mismatch counts")
    log("  here always agree with Step 4/5 - see realize_duplex() in this script.")
    log(f"  {n_window_n} REDUCED worst-component rows have a real N in the genome")
    log("  window (sequencing ambiguity). Realized as a non-penalizing perfect")
    log("  pair (Step 4's 'skip, not scoring' rule for genome N), not a mismatch.")
    log("  Tm shift predicts risk; it does not prove assay failure. Wet-lab")
    log("  confirmation remains the final word.")

    (OUT_DIR / "step6_report.txt").write_text("\n".join(report))
    log()
    log("=" * 70)
    log(f"WROTE  {OUT_DIR}/thermodynamics_per_component.csv")
    log(f"WROTE  {OUT_DIR}/thermodynamics_summary.csv")
    log(f"WROTE  {OUT_DIR}/likely_failure_mismatch_summary.csv")
    log(f"WROTE  {OUT_DIR}/step6_report.txt")
    log("=" * 70)


if __name__ == "__main__":
    main()
