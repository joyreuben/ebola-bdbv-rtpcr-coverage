#!/usr/bin/env python3
"""
Supplementary analysis: refetch the 8 records excluded at Step 1 and run
them through the same Step 4/5 pipeline logic, so a supplementary
"all 78 records" version of Table 1 can be shown alongside the main
70-genome result. This does NOT change any primary output (step1/4/5
directories are untouched) - it reuses their functions against a
separate, smaller genome set and writes to step1_output/supp_*.

WHY THIS EXISTS
    A reviewer should be able to see exactly what the patent-record
    exclusion (CLAUDE.md, Step 1) changed, not just take the exclusion's
    correctness on faith. The 8 excluded sequences were never saved to
    disk after exclusion (only the metadata audit-trail rows survive), so
    this refetches them from NCBI by accession.

Requires: pandas, biopython, network access to eutils.ncbi.nlm.nih.gov
Reuses  : best_alignment() from step4_oligo_alignment.py
          ASSAY_STRUCTURE, classify(), best_component() from step5_combine_assays.py
Outputs : step1_output/supp_excluded_sequences.fasta
          step1_output/supp_all78_coverage.csv
          step1_output/supp_all78_report.txt
"""
import sys
import time
from pathlib import Path

import pandas as pd
from Bio import SeqIO, Entrez

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from step4_oligo_alignment import best_alignment  # noqa: E402
from step5_combine_assays import ASSAY_STRUCTURE, classify, best_component  # noqa: E402

Entrez.email = "samigold8650@gmail.com"

EXCLUDED = {
    # accession: (pathoplexus_id, exclusion_reason_short)
    "LQ278392.1": ("PP_006X66F.1", "133bp fragment, INSDC automated ingest"),
    "HC918458.1": ("PP_006X5S9.1", "GenBank Patent-division (HC)"),
    "HC918461.1": ("PP_006X5T7.1", "GenBank Patent-division (HC)"),
    "HV753272.1": ("PP_006X5U5.1", "GenBank Patent-division (HV)"),
    "HV753275.1": ("PP_006X5V3.1", "GenBank Patent-division (HV)"),
    "HZ164395.1": ("PP_006X5W1.1", "GenBank Patent-division (HZ)"),
    "HZ164398.1": ("PP_006X5XZ.1", "GenBank Patent-division (HZ)"),
    "LQ278396.1": ("PP_006X67D.1", "173bp fragment, GenBank Patent-division (LQ)"),
}

OUT_DIR = Path("step1_output")
FASTA_OUT = OUT_DIR / "supp_excluded_sequences.fasta"
PANEL_FILE = Path("step2_output/validated_panel.csv")


def fetch_excluded():
    if FASTA_OUT.exists():
        print(f"Using cached {FASTA_OUT}")
        return {rec.id.split(".")[0]: str(rec.seq).upper()
                for rec in SeqIO.parse(str(FASTA_OUT), "fasta")}
    seqs = {}
    with open(FASTA_OUT, "w") as fh:
        for acc in EXCLUDED:
            print(f"  fetching {acc} ...")
            handle = Entrez.efetch(db="nuccore", id=acc, rettype="fasta", retmode="text")
            text = handle.read()
            handle.close()
            fh.write(text)
            if not text.startswith(">"):
                sys.exit(f"Unexpected response for {acc}: {text[:200]}")
            time.sleep(0.4)  # NCBI rate-limit courtesy
    return {rec.id.split(".")[0]: str(rec.seq).upper()
            for rec in SeqIO.parse(str(FASTA_OUT), "fasta")}


def main():
    print("Fetching 8 excluded records from NCBI ...")
    excluded_seqs = fetch_excluded()
    print(f"  got {len(excluded_seqs)} sequences")
    for acc, seq in excluded_seqs.items():
        print(f"    {acc}: {len(seq)} bp")

    panel = pd.read_csv(PANEL_FILE, dtype=str)
    panel.columns = panel.columns.str.strip()

    print("\nAligning all 22 oligos against the 8 excluded records ...")
    rows = []
    for _, row in panel.iterrows():
        oligo = str(row["sequence_5to3"]).upper().strip()
        for gid, gseq in excluded_seqs.items():
            b = best_alignment(oligo, gseq)
            rows.append({
                "assay_id": row["assay_id"],
                "genome_id": gid,
                "best_start": b["start"],
                "strand": b["strand"],
                "mismatches": b["mismatches"],
                "three_prime_mismatches": b["three_prime_mismatches"],
                "skipped_positions": b["skipped_positions"],
                "oligo_length": len(oligo),
            })
    excl_full = pd.DataFrame(rows)

    print("Classifying via the same Step 5 ASSAY_STRUCTURE / classify() logic ...")
    cov_rows = []
    for assay_name, struct in ASSAY_STRUCTURE.items():
        comp_best = {}
        for role in ("forward", "reverse", "probe"):
            ids = struct.get(role, [])
            if ids:
                comp_best[role] = best_component(excl_full, ids)
        for g in excluded_seqs:
            required = struct["required"]
            comp_mm, comp_3p, comp_skip = {}, {}, {}
            usable = True
            for role in required:
                cb = comp_best.get(role)
                if cb is None or g not in cb.index:
                    usable = False
                    break
                comp_mm[role] = int(cb.loc[g, "mm"])
                comp_3p[role] = int(cb.loc[g, "mm_3p"])
                comp_skip[role] = float(cb.loc[g, "skip_frac"])
            if not usable:
                continue
            worst_role = max(comp_mm, key=lambda r: (comp_mm[r], comp_3p[r]))
            worst_mm = comp_mm[worst_role]
            worst_3p = comp_3p[worst_role]
            worst_skip = max(comp_skip.values())
            call = classify(worst_mm, worst_3p, worst_skip)
            cov_rows.append({
                "assay": assay_name, "genome_id": g,
                "worst_mm": worst_mm, "worst_3p_mm": worst_3p,
                "worst_skip_fraction": round(worst_skip, 3),
                "classification": call,
            })
    excl_cov = pd.DataFrame(cov_rows)

    # ---- merge with the main 70-genome coverage for a full 78-record view
    main_cov = pd.read_csv("step5_output/assay_coverage_per_genome.csv")
    main_cov = main_cov.rename(columns={"worst_mm": "worst_mm", "worst_3p_mm": "worst_3p_mm"})
    main_slim = main_cov[["assay", "genome_id", "worst_mm", "worst_3p_mm",
                           "worst_skip_fraction", "classification"]].copy()
    main_slim["excluded"] = False
    excl_cov["excluded"] = True
    id_to_reason = {v[0].split(".")[0]: (k, v[1]) for k, v in EXCLUDED.items()}
    excl_cov["accession"] = excl_cov["genome_id"].map(lambda g: id_to_reason.get(g, (g, ""))[0])
    excl_cov["exclusion_reason"] = excl_cov["genome_id"].map(lambda g: id_to_reason.get(g, (g, ""))[1])

    all78 = pd.concat([main_slim, excl_cov], ignore_index=True, sort=False)
    all78.to_csv(OUT_DIR / "supp_all78_coverage.csv", index=False)

    # ---- report: per assay, tier counts with vs without the 8 excluded records
    report = []

    def log(s=""):
        print(s)
        report.append(s)

    log("=" * 78)
    log("SUPPLEMENTARY: all 78 downloaded records (70 curated + 8 excluded),")
    log("per-assay classification, side by side")
    log("=" * 78)
    for assay in ASSAY_STRUCTURE:
        log(f"\n{assay}")
        sub70 = main_slim[main_slim["assay"] == assay]
        sub8 = excl_cov[excl_cov["assay"] == assay]
        log("  70 curated (main dataset):")
        log("    " + sub70["classification"].value_counts().to_string().replace("\n", "\n    "))
        log("  8 excluded (this supplementary check):")
        if sub8.empty:
            log("    (no evaluable rows)")
        else:
            for _, r in sub8.sort_values("genome_id").iterrows():
                log(f"    {r['accession']:<14} mm={r['worst_mm']:>2}  3p_mm={r['worst_3p_mm']}  "
                    f"skip={r['worst_skip_fraction']:.3f}  -> {r['classification']}")

    (OUT_DIR / "supp_all78_report.txt").write_text("\n".join(report))
    log(f"\nWROTE {OUT_DIR}/supp_excluded_sequences.fasta")
    log(f"WROTE {OUT_DIR}/supp_all78_coverage.csv")
    log(f"WROTE {OUT_DIR}/supp_all78_report.txt")


if __name__ == "__main__":
    main()
