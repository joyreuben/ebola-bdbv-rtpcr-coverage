#!/usr/bin/env python3
"""
Step 5: Combine per-oligo results into per-ASSAY coverage per genome.

WHY THIS EXISTS
    Step 4 scored each oligo separately. But a real PCR assay succeeds or
    fails as a unit, and several of these assays use PRIMER MIXES:
      - Panning FiloA2.2 / 2.3 / 2.4 are three forward primers mixed in one
        tube. The assay's forward performance = the BEST of the three.
      - Panning FiloB / FiloB-Ra are two reverse primers mixed. Reverse
        performance = the BEST (minimum mismatch) of the two.
    Scoring these individually OVERSTATES failure. This script applies the
    correct "best-of-mix" logic, per genome, to get true assay coverage.

WHAT IT DOES
    1. Groups oligos into functional roles per assay (forward / reverse /
       probe), honouring primer mixes defined in ASSAY_STRUCTURE below.
    2. For each assay and genome, computes:
         - forward_mm  = min mismatch across the forward primer mix
         - reverse_mm  = min mismatch across the reverse primer mix
         - probe_mm    = min mismatch across probe(s), if any
         - worst_component_mm = the limiting component (the assay is only
           as good as its weakest required part)
         - 3'-end mismatches for the chosen best forward and reverse
    3. Classifies each assay-genome pair by a transparent, cited rule set
       into: LIKELY DETECTED / REDUCED / AT RISK / LIKELY FAILURE, or
       NOT EVALUABLE when a required component's binding site is too
       unresolved (N/gap) in that genome to call at all.
    4. Summarises coverage per assay across the three outbreak arms.

    Classification thresholds are declared up front and are EDITABLE. They
    are decision aids, not measurements - state them in your methods and
    report the raw mismatch numbers alongside every call.

Requires: pandas
Inputs : step4_output/oligo_vs_genome_full.csv
         step1_output/curated_metadata.csv   (for outbreak arm per genome)
Outputs: step5_output/assay_coverage_per_genome.csv
         step5_output/assay_coverage_summary.csv
         step5_output/step5_report.txt
"""

import sys
from pathlib import Path
import pandas as pd

FULL = Path("step4_output/oligo_vs_genome_full.csv")
META = Path("step1_output/curated_metadata.csv")
OUT_DIR = Path("step5_output")

# --------------------------------------------------------------------
# ASSAY STRUCTURE
# Define which oligos belong to which assay and functional role.
# Oligos listed together in a role are a MIX (best-of scoring).
# 'required' lists the roles that must work for the assay to detect.
#   Conventional PCR: forward + reverse required, probe optional.
#   TaqMan: forward + reverse + probe all required.
# --------------------------------------------------------------------
ASSAY_STRUCTURE = {
    "Panning2007_panfilo": {
        "type": "TaqMan (EBOV/SUDV target)",
        "forward": ["PANNING2007_FiloA2.2", "PANNING2007_FiloA2.3",
                    "PANNING2007_FiloA2.4"],
        "reverse": ["PANNING2007_FiloB", "PANNING2007_FiloB-Ra"],
        # EBOV-targeting probes are the relevant ones for BDBV (an ebolavirus).
        # The Marburg probe is a specificity control - see
        # Panning2007_MARVprobe_control below, same primer mix.
        "probe":   ["PANNING2007_FAMEBOSu", "PANNING2007_FAMEBOg"],
        "required": ["forward", "reverse", "probe"],
    },
    "Panning2007_MARVprobe_control": {
        "type": "TaqMan (MARV probe control - should NOT detect BDBV)",
        # Same panfilovirus primer mix as Panning2007_panfilo; only the
        # probe differs. This checks that the Marburg-specific probe
        # channel does not fire on an ebolavirus sample.
        "forward": ["PANNING2007_FiloA2.2", "PANNING2007_FiloA2.3",
                    "PANNING2007_FiloA2.4"],
        "reverse": ["PANNING2007_FiloB", "PANNING2007_FiloB-Ra"],
        "probe":   ["PANNING2007_FAMMBG"],
        "required": ["forward", "reverse", "probe"],
    },
    "Ogawa2010_EBOV": {
        "type": "Conventional (EBOV typing)",
        "forward": ["OGAWA2010_EBOV_F"],
        "reverse": ["OGAWA2010_EBOV_R"],
        "probe":   [],
        "required": ["forward", "reverse"],
    },
    "SanchezFiloAB": {
        "type": "Conventional (pan-filovirus L)",
        "forward": ["SANCHEZ1999_FILO_A"],
        "reverse": ["PANNING2007_FiloB"],  # FILO-B == FiloB, reused
        "probe":   [],
        "required": ["forward", "reverse"],
    },
    "Weidmann2004_ENZ": {
        "type": "TaqMan (EBOV-Zaire typing)",
        "forward": ["WEIDMANN2004_ENZ_F"],
        "reverse": ["WEIDMANN2004_ENZ_R"],
        "probe":   ["WEIDMANN2004_ENZ_P"],
        "required": ["forward", "reverse", "probe"],
    },
    "Weidmann2004_ENS": {
        "type": "TaqMan (EBOV-Sudan typing)",
        "forward": ["WEIDMANN2004_ENS_F"],
        "reverse": ["WEIDMANN2004_ENS_R"],
        "probe":   ["WEIDMANN2004_ENS_P"],
        "required": ["forward", "reverse", "probe"],
    },
    "Weidmann2004_MN": {
        "type": "TaqMan (MARV control - should NOT detect BDBV)",
        "forward": ["WEIDMANN2004_MN_F"],
        "reverse": ["WEIDMANN2004_MN_R"],
        "probe":   ["WEIDMANN2004_MN_P"],
        "required": ["forward", "reverse", "probe"],
    },
    "Greninger2024_PfiloL": {
        "type": "Conventional/qPCR (pan-filovirus L)",
        "forward": ["GRENINGER2024_PfiloL_F"],
        "reverse": ["GRENINGER2024_PfiloL_R"],
        "probe":   [],
        "required": ["forward", "reverse"],
    },
}

# --------------------------------------------------------------------
# CLASSIFICATION RULES  (editable; declare these in your methods)
# Based on the widely used principle that total mismatch load AND
# 3'-terminal mismatches both degrade PCR, with 3'-end mismatches
# far more damaging (Stadhouders et al. 2010; Lefever et al. 2013).
# These are decision aids, not measurements.
# --------------------------------------------------------------------

def classify(worst_mm: int, worst_3p: int, worst_skip_fraction: float) -> str:
    if worst_skip_fraction > 0.15:
        return "NOT EVALUABLE"       # binding site too unresolved (N/gap) to call
    if worst_3p >= 2 or worst_mm >= 6:
        return "LIKELY FAILURE"          # multiple 3'-end mismatches, or heavy total load
    if worst_3p == 1 and worst_mm >= 3:
        return "AT RISK"                 # 3'-end mismatch plus load
    if worst_mm >= 3:
        return "REDUCED"                 # moderate load, internal
    return "LIKELY DETECTED"             # none, or minor internal mismatches


def best_component(df, oligo_ids):
    """
    Given the full per-oligo/per-genome frame and a list of oligo ids that
    form a mix, return, per genome, the best (min mismatch) member, its
    3'-end mismatches, and its skipped-position fraction (N/gap bases in
    the binding site - a low mismatch count there is not evidence of a
    match, it's an absence of data). Returns a DataFrame indexed by
    genome_id.
    """
    sub = df[df["assay_id"].isin(oligo_ids)].copy()
    if sub.empty:
        return None
    sub["skip_frac"] = sub["skipped_positions"] / sub["oligo_length"]
    # Rank within each genome by (mismatches, three_prime_mismatches)
    sub = sub.sort_values(["genome_id", "mismatches", "three_prime_mismatches"])
    best = sub.groupby("genome_id").first()
    return best[["mismatches", "three_prime_mismatches", "skip_frac", "assay_id"]].rename(
        columns={"mismatches": "mm",
                 "three_prime_mismatches": "mm_3p",
                 "assay_id": "best_member"})


def main():
    if not FULL.exists():
        sys.exit(f"Missing {FULL}. Run Step 4 first.")
    OUT_DIR.mkdir(exist_ok=True)
    report = []

    def log(s=""):
        print(s)
        report.append(s)

    df = pd.read_csv(FULL)

    # Attach outbreak arm from Step 1 metadata, if available
    arm = {}
    if META.exists():
        meta = pd.read_csv(META)
        acc_col = next((c for c in meta.columns
                        if "accession" in c.lower()), meta.columns[0])
        ob_col = next((c for c in meta.columns
                       if "outbreak" in c.lower()), None)
        if ob_col:
            for _, r in meta.iterrows():
                key = str(r[acc_col]).split(".")[0]
                arm[key] = r[ob_col]

    log("=" * 70)
    log("STEP 5: MIXED-PRIMER-AWARE PER-ASSAY COVERAGE")
    log("=" * 70)
    log()
    log("Primer mixes scored best-of-members per genome. Assay limited by")
    log("its weakest REQUIRED component. Thresholds are decision aids -")
    log("raw mismatch numbers reported alongside every call.")
    log()

    rows = []
    for assay_name, struct in ASSAY_STRUCTURE.items():
        comp_best = {}
        for role in ("forward", "reverse", "probe"):
            ids = struct.get(role, [])
            if ids:
                comp_best[role] = best_component(df, ids)

        genomes = df["genome_id"].unique()
        for g in genomes:
            required = struct["required"]
            comp_mm, comp_3p, comp_skip, chosen = {}, {}, {}, {}
            usable = True
            for role in required:
                cb = comp_best.get(role)
                if cb is None or g not in cb.index:
                    usable = False
                    break
                comp_mm[role] = int(cb.loc[g, "mm"])
                comp_3p[role] = int(cb.loc[g, "mm_3p"])
                comp_skip[role] = float(cb.loc[g, "skip_frac"])
                chosen[role] = cb.loc[g, "best_member"]
            if not usable:
                continue

            worst_role = max(comp_mm, key=lambda r: (comp_mm[r], comp_3p[r]))
            worst_mm = comp_mm[worst_role]
            worst_3p = comp_3p[worst_role]
            # NOT EVALUABLE if ANY required component's best member sits on
            # a mostly-unresolved (N/gap) part of the genome - a low
            # mismatch count there reflects missing data, not a real match.
            worst_skip = max(comp_skip.values())
            call = classify(worst_mm, worst_3p, worst_skip)

            rows.append({
                "assay": assay_name,
                "assay_type": struct["type"],
                "genome_id": g,
                "outbreak": arm.get(g, "unknown"),
                "forward_mm": comp_mm.get("forward"),
                "reverse_mm": comp_mm.get("reverse"),
                "probe_mm": comp_mm.get("probe"),
                "worst_component": worst_role,
                "worst_mm": worst_mm,
                "worst_3p_mm": worst_3p,
                "worst_skip_fraction": round(worst_skip, 3),
                "classification": call,
                "chosen_forward": chosen.get("forward"),
                "chosen_reverse": chosen.get("reverse"),
                "chosen_probe": chosen.get("probe"),
            })

    cov = pd.DataFrame(rows)
    cov.to_csv(OUT_DIR / "assay_coverage_per_genome.csv", index=False)

    # ---- summary per assay x outbreak --------------------------------
    log("-" * 70)
    log("COVERAGE SUMMARY  (% of genomes in each class, per assay x arm)")
    log("-" * 70)
    log()
    summ_rows = []
    for (assay, ob), grp in cov.groupby(["assay", "outbreak"]):
        n = len(grp)
        counts = grp["classification"].value_counts()
        # NOT EVALUABLE rows carry worst_mm=0 as a placeholder (skipped
        # positions, not a real measurement) - averaging them in with real
        # mismatch counts silently drags the mean toward zero. Exclude them.
        evaluable = grp[grp["classification"] != "NOT EVALUABLE"]
        mean_mm = round(evaluable["worst_mm"].mean(), 2) if len(evaluable) else float("nan")
        row = {
            "assay": assay,
            "outbreak": ob,
            "n_genomes": n,
            "mean_worst_mm": mean_mm,
            "detected": counts.get("LIKELY DETECTED", 0),
            "reduced": counts.get("REDUCED", 0),
            "at_risk": counts.get("AT RISK", 0),
            "failure": counts.get("LIKELY FAILURE", 0),
            "not_evaluable": counts.get("NOT EVALUABLE", 0),
        }
        summ_rows.append(row)
    summ = pd.DataFrame(summ_rows).sort_values(["assay", "outbreak"])
    summ.to_csv(OUT_DIR / "assay_coverage_summary.csv", index=False)
    log(summ.to_string(index=False))
    log()

    # ---- focused readout: the 2026 arm -------------------------------
    log("-" * 70)
    log("FOCUS  2026 outbreak arm - the key question")
    log("-" * 70)
    log()
    arm2026 = cov[cov["outbreak"].astype(str).str.contains("2026", na=False)]
    if len(arm2026):
        for assay, grp in arm2026.groupby("assay"):
            n = len(grp)
            calls = grp["classification"].value_counts().to_dict()
            evaluable = grp[grp["classification"] != "NOT EVALUABLE"]
            mean_mm = round(evaluable["worst_mm"].mean(), 2) if len(evaluable) else float("nan")
            log(f"  {assay} (n={n}, mean worst mm={mean_mm})")
            for cls in ["LIKELY DETECTED", "REDUCED", "AT RISK", "LIKELY FAILURE", "NOT EVALUABLE"]:
                if cls in calls:
                    log(f"      {cls}: {calls[cls]}")
            log()
    else:
        log("  No 2026-labelled genomes found - check outbreak metadata.")

    (OUT_DIR / "step5_report.txt").write_text("\n".join(report))
    log("=" * 70)
    log(f"WROTE  {OUT_DIR}/assay_coverage_per_genome.csv")
    log(f"WROTE  {OUT_DIR}/assay_coverage_summary.csv")
    log(f"WROTE  {OUT_DIR}/step5_report.txt")
    log("=" * 70)
    log()
    log("READ THIS: MARV control assay showing high failure vs BDBV is")
    log("CORRECT (good specificity). A pan-filovirus assay showing AT RISK")
    log("or FAILURE in the 2026 arm is your actual finding. Confirm those")
    log("with Step 6 thermodynamics before making claims.")


if __name__ == "__main__":
    main()