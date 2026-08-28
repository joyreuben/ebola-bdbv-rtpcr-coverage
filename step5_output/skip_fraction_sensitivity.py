#!/usr/bin/env python3
"""
Sensitivity check on the 0.15 worst_skip_fraction cutoff used to call
NOT EVALUABLE in step5_combine_assays.py's classify().

WHY THIS EXISTS
    The 0.15 threshold is not cited to an external source - it's a
    judgment call ("more than ~1 in 7 of the binding-site bases are N/gap
    is too unresolved to call"). Reviewers will ask whether the tier
    counts reported in Table 1 are an artifact of exactly where that line
    is drawn. This reruns the SAME classify() logic already used in Step 5
    (mismatch/3'-mismatch rules untouched) at 0.10 / 0.15 / 0.20 instead,
    using assay_coverage_per_genome.csv's already-computed worst_mm,
    worst_3p_mm, worst_skip_fraction per assay/genome row - no need to
    redo the Step 4/5 alignment or best-of-mix logic, since those numbers
    do not depend on the NOT EVALUABLE threshold at all.

Input : step5_output/assay_coverage_per_genome.csv
Output: step5_output/skip_fraction_sensitivity.csv, printed summary
"""
import pandas as pd

df = pd.read_csv("step5_output/assay_coverage_per_genome.csv")


def classify(worst_mm, worst_3p, worst_skip, cutoff):
    if worst_skip > cutoff:
        return "NOT EVALUABLE"
    if worst_3p >= 2 or worst_mm >= 6:
        return "LIKELY FAILURE"
    if worst_3p == 1 and worst_mm >= 3:
        return "AT RISK"
    if worst_mm >= 3:
        return "REDUCED"
    return "LIKELY DETECTED"


cutoffs = [0.10, 0.15, 0.20]
for c in cutoffs:
    df[f"classification_{c}"] = df.apply(
        lambda r: classify(r["worst_mm"], r["worst_3p_mm"], r["worst_skip_fraction"], c),
        axis=1)

df.to_csv("step5_output/skip_fraction_sensitivity.csv", index=False)

print("Tier counts per assay, at each skip-fraction cutoff")
print("=" * 70)
for assay, g in df.groupby("assay"):
    print(f"\n{assay}  (n={len(g)})")
    tab = pd.DataFrame({c: g[f"classification_{c}"].value_counts() for c in cutoffs}).fillna(0).astype(int)
    tab = tab.reindex(["LIKELY DETECTED", "REDUCED", "AT RISK", "LIKELY FAILURE", "NOT EVALUABLE"]).fillna(0).astype(int)
    print(tab.to_string())

print("\n" + "=" * 70)
print("Rows whose call CHANGES between 0.10 and 0.20:")
changed = df[df["classification_0.1"] != df["classification_0.2"]]
print(f"{len(changed)} of {len(df)} total assay-genome rows")
if len(changed):
    print(changed[["assay", "genome_id", "worst_skip_fraction",
                    "classification_0.1", "classification_0.15", "classification_0.2"]].to_string(index=False))
