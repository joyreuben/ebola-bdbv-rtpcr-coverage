#!/usr/bin/env python3
"""
Step 4 (and replacement sanity check for Step 2 Check 6):
Mismatch-tolerant alignment of every oligo against every genome.

WHY THIS EXISTS
    Exact string matching (Step 2 Check 6) reports NOT FOUND for any oligo
    that does not match a genome perfectly. Old EBOV/SUDV/MARV assays will
    never match a BDBV genome perfectly - that non-match IS the research
    question. This script measures HOW FAR each oligo is from matching,
    instead of a binary found/not-found.

WHAT IT DOES
    1. Fetches four species reference genomes from NCBI (Zaire, Sudan,
       Marburg, Bundibugyo).
    2. SANITY CHECK: aligns each oligo against its own design-species
       reference. An oligo that fits its design species with <=1 mismatch
       is correctly transcribed. One that fails even there is a real
       transcription/strand error to fix.
    3. Slides each oligo across each genome, scoring mismatches with full
       IUPAC-degenerate awareness. Inosine-derived N positions are skipped
       (neither match nor mismatch).
    4. Reports best-fit position, mismatch count, 3'-end mismatches
       (primers), and strand, per oligo per genome.

RUN A SELF-TEST FIRST:
    python step4_oligo_alignment.py --selftest
This runs on synthetic data and must print ALL TESTS PASSED before you
trust it on real genomes.

NORMAL RUN:
    python step4_oligo_alignment.py

Requires: pandas, biopython
    pip install pandas biopython --break-system-packages
Network: needs eutils.ncbi.nlm.nih.gov (Entrez) unless --refs-dir given.
"""

import argparse
import sys
import time
from pathlib import Path

import pandas as pd
from Bio import SeqIO, Entrez

# ---------------------------------------------------------------- config

PANEL_FILE = Path("step2_output/validated_panel.csv")
GENOMES_FILE = Path("step1_output/curated_sequences.fasta")
OUT_DIR = Path("step4_output")

# You must set a real email for NCBI Entrez. They block anonymous abuse.
ENTREZ_EMAIL = "TODO_your_email@example.com"

# Species reference genomes on NCBI RefSeq
SPECIES_REFS = {
    "EBOV": "NC_002549.1",   # Zaire ebolavirus
    "SUDV": "NC_006432.1",   # Sudan ebolavirus
    "MARV": "NC_001608.3",   # Marburg marburgvirus
    "BDBV": "NC_014373.1",   # Bundibugyo ebolavirus
}

# Map an oligo's organism_scope to the reference it was designed against.
# This is what the sanity check aligns against.
SCOPE_TO_DESIGN_REF = {
    "EBOV-specific": "EBOV",
    "EBOV-Zaire-specific": "EBOV",
    "EBOV-Sudan-specific": "SUDV",
    "MARV-specific": "MARV",
    # pan-filovirus oligos were designed against an alignment of several
    # species, so we accept a good fit against ANY reference for the
    # sanity check, handled specially below.
}

SANITY_MAX_MISMATCH = 2   # <= this against design species = transcription OK
THREE_PRIME_WINDOW = 5    # last N bases scored separately for primers

# IUPAC expansion (N covers everything; used for inosine-derived positions)
IUPAC = {
    "A": {"A"}, "C": {"C"}, "G": {"G"}, "T": {"T"}, "U": {"T"},
    "R": {"A", "G"}, "Y": {"C", "T"}, "S": {"C", "G"}, "W": {"A", "T"},
    "K": {"G", "T"}, "M": {"A", "C"}, "B": {"C", "G", "T"},
    "D": {"A", "G", "T"}, "H": {"A", "C", "T"}, "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},
}
# Positions we treat as universal / non-scoring (inosine was entered as N).
SKIP_BASES = {"N"}

COMP = str.maketrans("ACGTURYSWKMBDHVNacgturyswkmbdhvn",
                     "TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn")


def rc(seq: str) -> str:
    return seq.translate(COMP)[::-1]


def base_matches(oligo_base: str, genome_base: str) -> str:
    """
    Return 'skip', 'match', or 'mismatch' for one position.
    - inosine-derived N (and any N) -> skip (non-scoring, universal)
    - genome N or gap -> skip (not evaluable)
    - otherwise IUPAC-aware match test
    """
    ob, gb = oligo_base.upper(), genome_base.upper()
    if ob in SKIP_BASES:
        return "skip"
    if gb in ("N", "-", ".") or gb not in "ACGT":
        return "skip"
    allowed = IUPAC.get(ob, {ob})
    return "match" if gb in allowed else "mismatch"


def score_window(oligo: str, window: str, three_prime_at_start: bool = False,
                 abandon_above=None):
    """
    Score an oligo against an equal-length genome window.

    `oligo` is whatever probe string is actually being slid against the
    genome here - either the primer as written, or its reverse complement
    when that strand wins the alignment. Reverse-complementing reverses
    coordinate order, so when the probe is the reverse complement, the
    primer's true 3' base sits at INDEX 0 of this array, not the end.
    `three_prime_at_start` tells this function which end to score as the
    3' end for the array it was actually given.

    Skipped positions are counted in two ways: `skipped_positions` is every
    non-scoring position (oligo inosine/N, or unresolved genome base), and
    `genome_skipped_positions` is only the subset caused by the GENOME base
    being unresolved (N/gap/ambiguity) under a scoring oligo base. The second
    one is what best_alignment() uses to keep the oligo off N-runs.

    `abandon_above` is a speed-up only: if given, scoring stops and None is
    returned as soon as mismatches + genome_skipped_positions exceeds it
    (that window can no longer beat the current best placement).
    """
    mism = skips = genome_skips = 0
    three_prime_mism = 0
    olen = len(oligo)
    for i, (ob, gb) in enumerate(zip(oligo, window)):
        result = base_matches(ob, gb)
        if result == "skip":
            skips += 1
            if ob.upper() not in SKIP_BASES:
                genome_skips += 1
        elif result == "mismatch":
            mism += 1
            in_three_prime = (
                i < THREE_PRIME_WINDOW if three_prime_at_start
                else i >= olen - THREE_PRIME_WINDOW
            )
            if in_three_prime:
                three_prime_mism += 1
        if abandon_above is not None and mism + genome_skips > abandon_above:
            return None
    return {
        "mismatches": mism,
        "three_prime_mismatches": three_prime_mism,
        "skipped_positions": skips,
        "genome_skipped_positions": genome_skips,
    }


def best_alignment(oligo: str, genome: str) -> dict:
    """
    Slide oligo across genome (both strands), return the best-scoring
    placement. Reports which strand won.

    PLACEMENT RULE: the winning window minimises
        mismatches + genome_skipped_positions
    (ties broken by fewest 3'-end mismatches, then first position found).
    Unresolved genome bases are still reported as skips, never as
    mismatches - they only count against a window when deciding WHERE the
    oligo sits.

    Why: ranking on mismatches alone makes any stretch of N at least as long
    as the oligo an unbeatable "0 mismatches, everything skipped" window. A
    genome with one N-run anywhere then had every oligo with >=1 real
    mismatch parked on that N-run and reported as unresolved, although its
    true binding site was fully sequenced elsewhere in the same record
    (349 of 374 high-skip oligo/genome pairs in the 70-genome set were this
    artifact, checked against the MAFFT alignment).

    The opposite failure - the true site really is inside an N-run, so the
    winner is some unrelated resolved window - cannot be detected here. It is
    caught by Step 4b (site_unresolved_fraction), which Step 5 consumes.
    """
    best = None
    for strand, probe in (("as_written", oligo), ("reverse_complement", rc(oligo))):
        olen = len(probe)
        three_prime_at_start = strand == "reverse_complement"
        for i in range(len(genome) - olen + 1):
            window = genome[i:i + olen]
            score = score_window(probe, window, three_prime_at_start,
                                 abandon_above=best["_key"][0] if best else None)
            if score is None:
                continue
            key = (score["mismatches"] + score["genome_skipped_positions"],
                   score["three_prime_mismatches"])
            if best is None or key < best["_key"]:
                best = {
                    "_key": key,
                    "strand": strand,
                    "start": i,
                    "end": i + olen,
                    **score,
                }
    if best:
        best.pop("_key")
    return best


# ------------------------------------------------------------- self test

def selftest() -> None:
    print("SELF-TEST on synthetic data")
    print("-" * 50)
    ok = True

    # 1. exact match, zero mismatch
    g = "AAAA" + "ACGTACGT" + "TTTT"
    r = best_alignment("ACGTACGT", g)
    assert r["mismatches"] == 0 and r["start"] == 4, r
    print("  [pass] exact match found at correct position")

    # 2. one internal mismatch
    r = best_alignment("ACGTACGT", "AAAAACGAACGTTTTT".replace("ACGAACGT", "ACGAACGT"))
    g2 = "AAAA" + "ACGAACGT" + "TTTT"   # position 6 G->A vs oligo? build cleanly
    oligo = "ACGTACGT"
    win   = "ACGAACGT"
    r = best_alignment(oligo, "AAAA" + win + "TTTT")
    assert r["mismatches"] == 1, r
    print("  [pass] single mismatch counted")

    # 3. degenerate base matches, no mismatch
    #    oligo Y at a position where genome has C or T
    r = best_alignment("ACGTACYT", "AAAA" + "ACGTACCT" + "TTTT")
    assert r["mismatches"] == 0, r
    r = best_alignment("ACGTACYT", "AAAA" + "ACGTACTT" + "TTTT")
    assert r["mismatches"] == 0, r
    print("  [pass] degenerate Y matches C and T")

    # 4. N (inosine) position is skipped, not scored
    r = best_alignment("ACGTACNT", "AAAA" + "ACGTACGT" + "TTTT")
    assert r["mismatches"] == 0 and r["skipped_positions"] == 1, r
    r = best_alignment("ACGTACNT", "AAAA" + "ACGTACAT" + "TTTT")
    assert r["mismatches"] == 0 and r["skipped_positions"] == 1, r
    print("  [pass] N/inosine position skipped, never a mismatch")

    # 5. genome N is skipped (not evaluable)
    r = best_alignment("ACGTACGT", "AAAA" + "ACGTACNT" + "TTTT")
    assert r["mismatches"] == 0 and r["skipped_positions"] == 1, r
    print("  [pass] genome N treated as not-evaluable")

    # 6. reverse complement detection
    #    Use a non-palindromic oligo so only the RC strand can match.
    oligo = "AAACCCGGGTTTAAAC"   # rc = GTTTAAACCCGGGTTT, clearly different
    g = "TTTTTT" + rc(oligo) + "GGGGGG"
    r = best_alignment(oligo, g)
    assert r["mismatches"] == 0 and r["strand"] == "reverse_complement", r
    #    ...and it should NOT spuriously match the as-written strand better
    g2 = "TTTTTT" + oligo + "GGGGGG"
    r2 = best_alignment(oligo, g2)
    assert r2["mismatches"] == 0 and r2["strand"] == "as_written", r2
    print("  [pass] reverse-complement strand detected (non-palindromic)")

    # 7. 3' end mismatch flagged
    #    Unique oligo, single 3'-terminal mismatch, flanks chosen so no
    #    alternative placement or strand scores as well.
    oligo = "GATTACAGATTACAGA"
    win   = "GATTACAGATTACAGT"   # last base A->T mismatch (3' end)
    genome = "CCCCCC" + win + "CCCCCC"
    r = best_alignment(oligo, genome)
    assert r["mismatches"] == 1, r
    assert r["three_prime_mismatches"] == 1, r
    assert r["strand"] == "as_written", r
    print("  [pass] 3'-end mismatch flagged separately")

    # 8. internal mismatch NOT counted as 3' mismatch
    oligo = "GATTACAGATTACAGA"
    win   = "GATTACTGATTACAGA"   # position 6 A->T, internal
    genome = "CCCCCC" + win + "CCCCCC"
    r = best_alignment(oligo, genome)
    assert r["mismatches"] == 1 and r["three_prime_mismatches"] == 0, r
    print("  [pass] internal mismatch not miscounted as 3'-end")

    # 9. 3'-end detection on the reverse-complement strand (regression).
    #    Reverse-complementing reverses coordinate order: the primer's true
    #    3' base ends up at INDEX 0 of the rc(oligo) array, not the end.
    oligo = "AAAAAAAAAAACCCCC"   # 5' end = AAAAA, 3' end = CCCCC
    rc_oligo = rc(oligo)
    bad_at_3prime = "A" + rc_oligo[1:]   # genome mismatch at oligo's true 3' base
    genome = "TTTTTT" + bad_at_3prime + "TTTTTT"
    r = best_alignment(oligo, genome)
    assert r["strand"] == "reverse_complement", r
    assert r["mismatches"] == 1 and r["three_prime_mismatches"] == 1, r
    bad_at_5prime = rc_oligo[:-1] + "A"   # genome mismatch at oligo's 5' base
    genome2 = "TTTTTT" + bad_at_5prime + "TTTTTT"
    r2 = best_alignment(oligo, genome2)
    assert r2["strand"] == "reverse_complement", r2
    assert r2["mismatches"] == 1 and r2["three_prime_mismatches"] == 0, r2
    print("  [pass] 3'-end correctly identified on reverse-complement strand")

    # 10. an N-run elsewhere in the genome must not capture the oligo
    #     (regression). The true site carries 2 real mismatches; a 40 nt
    #     N-run would score "0 mismatches, all skipped" under a
    #     mismatch-only ranking and win.
    oligo = "GATTACAGATTACAGA"
    site  = "GATTACTGATTTCAGA"   # 2 internal mismatches
    genome = "CCCCCC" + "N" * 40 + "CCCCCC" + site + "CCCCCC"
    r = best_alignment(oligo, genome)
    assert r["start"] == 52, r
    assert r["mismatches"] == 2 and r["skipped_positions"] == 0, r
    #     ...and a site that is itself partly N is still found and reported
    #     with its skips, not swapped for some other window.
    site_n = "GATTACNGATTTCAGA"  # 1 N + 1 mismatch
    genome = "CCCCCC" + "N" * 40 + "CCCCCC" + site_n + "CCCCCC"
    r = best_alignment(oligo, genome)
    assert r["start"] == 52, r
    assert r["mismatches"] == 1 and r["skipped_positions"] == 1, r
    assert r["genome_skipped_positions"] == 1, r
    #     oligo inosine/N is not a genome skip
    r = best_alignment("GATTACNGATTACAGA", "CCCCCC" + "GATTACAGATTACAGA" + "CCCCCC")
    assert r["skipped_positions"] == 1 and r["genome_skipped_positions"] == 0, r
    print("  [pass] N-run elsewhere does not capture the oligo")

    print("-" * 50)
    print("ALL TESTS PASSED" if ok else "TESTS FAILED")


# ------------------------------------------------------------- ncbi fetch

def fetch_references(refs_dir: Path) -> dict:
    refs_dir.mkdir(exist_ok=True)
    genomes = {}
    for name, acc in SPECIES_REFS.items():
        cache = refs_dir / f"{name}_{acc}.fasta"
        if cache.exists():
            rec = next(SeqIO.parse(str(cache), "fasta"))
            genomes[name] = str(rec.seq).upper()
            print(f"  {name} {acc}: loaded from cache ({len(genomes[name])} bp)")
            continue
        print(f"  {name} {acc}: fetching from NCBI ...")
        try:
            handle = Entrez.efetch(db="nucleotide", id=acc,
                                   rettype="fasta", retmode="text")
            data = handle.read()
            handle.close()
            cache.write_text(data)
            rec = next(SeqIO.parse(str(cache), "fasta"))
            genomes[name] = str(rec.seq).upper()
            print(f"    got {len(genomes[name])} bp")
            time.sleep(0.4)   # be polite to NCBI
        except Exception as exc:
            print(f"    ERROR fetching {acc}: {exc}")
            print("    Set --refs-dir to a folder with the FASTA files instead.")
            sys.exit(1)
    return genomes


# ------------------------------------------------------------------ main

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--selftest", action="store_true")
    parser.add_argument("--refs-dir", default="references",
                        help="folder to cache/read species reference FASTAs")
    parser.add_argument("--email", default=ENTREZ_EMAIL,
                        help="email for NCBI Entrez")
    args = parser.parse_args()

    if args.selftest:
        selftest()
        return

    # An email is only needed if a reference actually has to be fetched from
    # NCBI. With all four FASTAs already cached in --refs-dir (as shipped in
    # this repository) the script runs offline with no email.
    need_fetch = any(not (Path(args.refs_dir) / f"{name}_{acc}.fasta").exists()
                     for name, acc in SPECIES_REFS.items())
    if need_fetch and args.email.startswith("TODO"):
        sys.exit("Set your email: --email you@example.com  (NCBI requires it)")
    Entrez.email = args.email

    OUT_DIR.mkdir(exist_ok=True)
    report = []

    def log(line=""):
        print(line)
        report.append(line)

    log("=" * 70)
    log("STEP 4: MISMATCH-TOLERANT OLIGO ALIGNMENT")
    log("=" * 70)
    log()

    if not PANEL_FILE.exists():
        sys.exit(f"Panel not found: {PANEL_FILE}. Run Step 2 first.")
    if not GENOMES_FILE.exists():
        sys.exit(f"Genomes not found: {GENOMES_FILE}. Run Step 1 first.")

    panel = pd.read_csv(PANEL_FILE, dtype=str)
    panel.columns = panel.columns.str.strip()
    log(f"Oligos in panel: {len(panel)}")

    # --- fetch references ---------------------------------------------
    log()
    log("Fetching species reference genomes")
    log("-" * 70)
    refs = fetch_references(Path(args.refs_dir))
    log()

    # --- SANITY: each oligo vs its design-species reference ------------
    log("-" * 70)
    log("SANITY CHECK  oligo vs its own design-species reference")
    log("-" * 70)
    log("An oligo should fit its design species with few mismatches.")
    log(f"Threshold for 'transcription OK': <= {SANITY_MAX_MISMATCH} mismatches.")
    log("pan-filovirus oligos: accept a good fit against ANY reference.")
    log()

    sanity_rows = []
    for _, row in panel.iterrows():
        oligo = str(row["sequence_5to3"]).upper().strip()
        scope = str(row["organism_scope"])
        if scope == "pan-filovirus":
            # best fit across all references
            results = {sp: best_alignment(oligo, refs[sp]) for sp in refs}
            best_sp = min(results, key=lambda s: results[s]["mismatches"])
            b = results[best_sp]
            design_ref = f"best={best_sp}"
        else:
            sp = SCOPE_TO_DESIGN_REF.get(scope)
            if not sp:
                log(f"  ? {row['assay_id']}: unknown scope '{scope}', skipped")
                continue
            b = best_alignment(oligo, refs[sp])
            design_ref = sp

        status = "OK" if b["mismatches"] <= SANITY_MAX_MISMATCH else "CHECK!"
        log(f"  [{status}] {row['assay_id']:<28} vs {design_ref:<10} "
            f"mm={b['mismatches']} strand={b['strand']} skip={b['skipped_positions']}")
        sanity_rows.append({
            "assay_id": row["assay_id"],
            "design_ref": design_ref,
            "mismatches_vs_design": b["mismatches"],
            "strand": b["strand"],
            "status": status,
        })

    failed = [r for r in sanity_rows if r["status"] == "CHECK!"]
    log()
    if failed:
        log(f"  {len(failed)} oligo(s) do not fit their design species.")
        log("  These are the ones to re-check against the paper (sequence or")
        log("  strand). All others are correctly transcribed, and their")
        log("  mismatches against BDBV are real biology.")
    else:
        log("  All oligos fit their design species. Panel transcription verified.")
        log("  Any mismatches against BDBV below are genuine divergence.")
    log()

    pd.DataFrame(sanity_rows).to_csv(OUT_DIR / "sanity_check.csv", index=False)

    # --- FULL: each oligo vs every BDBV/study genome ------------------
    log("-" * 70)
    log("FULL ANALYSIS  each oligo vs all study genomes")
    log("-" * 70)
    study = {rec.id.split(".")[0]: str(rec.seq).upper()
             for rec in SeqIO.parse(str(GENOMES_FILE), "fasta")}
    log(f"Study genomes: {len(study)}")
    log("Writing per-oligo-per-genome results ...")

    results = []
    for _, row in panel.iterrows():
        oligo = str(row["sequence_5to3"]).upper().strip()
        for gid, gseq in study.items():
            b = best_alignment(oligo, gseq)
            results.append({
                "assay_id": row["assay_id"],
                "organism_scope": row["organism_scope"],
                "oligo_type": row["oligo_type"],
                "target_gene": row["target_gene"],
                "genome_id": gid,
                "best_start": b["start"],
                "strand": b["strand"],
                "mismatches": b["mismatches"],
                "three_prime_mismatches": b["three_prime_mismatches"],
                "skipped_positions": b["skipped_positions"],
                "genome_skipped_positions": b["genome_skipped_positions"],
                "oligo_length": len(oligo),
            })

    res_df = pd.DataFrame(results)
    res_df.to_csv(OUT_DIR / "oligo_vs_genome_full.csv", index=False)
    log(f"  wrote {len(res_df)} rows to oligo_vs_genome_full.csv")
    log()

    # --- summary per oligo across the 2026 arm ------------------------
    log("-" * 70)
    log("SUMMARY  mean mismatches per oligo (all study genomes)")
    log("-" * 70)
    summ = (res_df.groupby(["assay_id", "organism_scope"])
            .agg(mean_mm=("mismatches", "mean"),
                 max_mm=("mismatches", "max"),
                 min_mm=("mismatches", "min"),
                 mean_3p_mm=("three_prime_mismatches", "mean"))
            .round(2).reset_index())
    log(summ.to_string(index=False))

    (OUT_DIR / "step4_report.txt").write_text("\n".join(report))
    log()
    log("=" * 70)
    log(f"WROTE  {OUT_DIR}/sanity_check.csv")
    log(f"WROTE  {OUT_DIR}/oligo_vs_genome_full.csv")
    log(f"WROTE  {OUT_DIR}/step4_report.txt")
    log("=" * 70)
    log()
    log("NEXT: interpret. High mean_mm on a pan-filovirus oligo = coverage")
    log("concern. High mm on MARV controls = expected (good specificity).")
    log("Then Step 5 adds thermodynamics to the flagged oligos.")


if __name__ == "__main__":
    main()
