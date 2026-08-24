#!/usr/bin/env python3
"""
Step 2: Validate the assay panel CSV.

Run this script TWICE:
    Pass 1 — immediately, on the template, to see what is missing.
    Pass 2 — after you have filled in all the sequences, to confirm
              everything is ready for alignment.

Input:
    assay_panel.csv   (copy of template with sequences filled in)

Output:
    step2_output/validated_panel.csv   clean panel ready for Step 3
    step2_output/step2_report.txt

Requires: pandas, biopython
    pip install pandas biopython --break-system-packages
"""

import sys
import re
from pathlib import Path
from itertools import product

import pandas as pd
from Bio import SeqIO

# ---------------------------------------------------------------- config

PANEL_FILE = Path("assay_panel.csv")
REFERENCE_FASTA = Path("step1_output/curated_sequences.fasta")
OUT_DIR = Path("step2_output")

# Valid IUPAC nucleotide codes (upper case)
IUPAC_CODES = set("ACGTURYSWKMBDHVN")

# IUPAC degenerate base expansion
IUPAC_EXPAND = {
    "A": ["A"], "C": ["C"], "G": ["G"], "T": ["T"], "U": ["T"],
    "R": ["A", "G"], "Y": ["C", "T"], "S": ["C", "G"],
    "W": ["A", "T"], "K": ["G", "T"], "M": ["A", "C"],
    "B": ["C", "G", "T"], "D": ["A", "G", "T"],
    "H": ["A", "C", "T"], "V": ["A", "C", "G"],
    "N": ["A", "C", "G", "T"],
}

EXPECTED_TYPES = {"forward", "reverse", "probe"}

# Sanity thresholds (not hard cutoffs — flag only)
MIN_PRIMER_LEN = 16
MAX_PRIMER_LEN = 35
MIN_PROBE_LEN  = 18
MAX_PROBE_LEN  = 40

# ---------------------------------------------------------------- helpers

def expand_degenerate(seq: str) -> list:
    """Return all concrete sequences an IUPAC-coded oligo can match."""
    options = [IUPAC_EXPAND.get(b.upper(), [b.upper()]) for b in seq]
    return ["".join(combo) for combo in product(*options)]


def reverse_complement(seq: str) -> str:
    comp = str.maketrans("ACGTURYSWKMBDHVNacgturyswkmbdhvn",
                         "TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn")
    return seq.translate(comp)[::-1]


def iupac_match(oligo: str, target: str) -> bool:
    """True if any expansion of oligo matches target (case-insensitive)."""
    target = target.upper()
    oligo  = oligo.upper()
    if len(oligo) != len(target):
        return False
    for exp in expand_degenerate(oligo):
        if exp == target:
            return True
    return False


def locate_in_genome(oligo: str, genome: str, is_reverse: bool = False) -> list:
    """
    Return list of (start, end) positions where oligo matches genome.
    Checks both the oligo and its reverse complement.
    Returns empty list if not found.
    """
    genome = genome.upper()
    search_seq = reverse_complement(oligo.upper()) if is_reverse else oligo.upper()
    expansions = expand_degenerate(search_seq)
    positions  = []
    olen = len(oligo)
    for i in range(len(genome) - olen + 1):
        window = genome[i : i + olen]
        for exp in expansions:
            if exp == window:
                positions.append((i, i + olen))
                break
    return positions


def validate_sequence(seq: str, label: str) -> list:
    """Return list of error strings, empty if clean."""
    errors = []
    seq = seq.strip().upper()
    bad = set(seq) - IUPAC_CODES
    if bad:
        errors.append(f"{label}: illegal characters {bad}")
    if not seq:
        errors.append(f"{label}: empty sequence")
    return errors


# ------------------------------------------------------------------ main

def main() -> None:
    OUT_DIR.mkdir(exist_ok=True)
    report = []

    def log(line: str = "") -> None:
        print(line)
        report.append(line)

    log("=" * 70)
    log("STEP 2: ASSAY PANEL VALIDATION")
    log("=" * 70)
    log()

    if not PANEL_FILE.exists():
        sys.exit(f"Panel file not found: {PANEL_FILE.resolve()}\n"
                 "Copy assay_panel_template.csv to assay_panel.csv and fill in sequences.")

    df = pd.read_csv(PANEL_FILE, dtype=str)
    df.columns = df.columns.str.strip()
    log(f"Rows loaded: {len(df)}")
    log()

    required_cols = ["assay_id", "publication_doi", "year_designed",
                     "target_gene", "oligo_type", "sequence_5to3"]
    missing_cols = [c for c in required_cols if c not in df.columns]
    if missing_cols:
        sys.exit(f"Missing columns: {missing_cols}")

    # ---- CHECK 1: TODO placeholders ----------------------------------
    log("-" * 70)
    log("CHECK 1  sequences filled in")
    log("-" * 70)
    todo_rows = df[df["sequence_5to3"].astype(str).str.startswith("TODO")]
    if len(todo_rows):
        log(f"  {len(todo_rows)} rows still have TODO placeholders:")
        for _, row in todo_rows.iterrows():
            log(f"    {row['assay_id']}  doi:{row['publication_doi']}")
        log()
        log("  Fill these in from the papers before proceeding.")
        log("  Run the script again when done.")
        log()
        log("  CANNOT CONTINUE — panel incomplete.")
        (OUT_DIR / "step2_report.txt").write_text("\n".join(report))
        sys.exit(0)
    else:
        log("  All sequences filled in.")
    log()

    # ---- CHECK 2: IUPAC validity -------------------------------------
    log("-" * 70)
    log("CHECK 2  IUPAC character validity")
    log("-" * 70)
    errors_found = False
    for _, row in df.iterrows():
        errs = validate_sequence(str(row["sequence_5to3"]), row["assay_id"])
        for e in errs:
            log(f"  ERROR: {e}")
            errors_found = True
    if not errors_found:
        log("  All sequences contain only valid IUPAC characters.")
    log()

    # ---- CHECK 3: oligo type and length ------------------------------
    log("-" * 70)
    log("CHECK 3  oligo types and lengths")
    log("-" * 70)
    bad_type = df[~df["oligo_type"].isin(EXPECTED_TYPES)]
    if len(bad_type):
        for _, row in bad_type.iterrows():
            log(f"  ERROR: unknown oligo_type '{row['oligo_type']}' "
                f"in {row['assay_id']}. Use: forward / reverse / probe")
    else:
        log("  All oligo types valid.")
    log()

    for _, row in df.iterrows():
        seq = str(row["sequence_5to3"]).upper().strip()
        otype = str(row["oligo_type"])
        length = len(seq)
        if otype in ("forward", "reverse"):
            if length < MIN_PRIMER_LEN or length > MAX_PRIMER_LEN:
                log(f"  FLAG: {row['assay_id']} primer length {length} "
                    f"(expected {MIN_PRIMER_LEN}-{MAX_PRIMER_LEN}). Verify from paper.")
        elif otype == "probe":
            if length < MIN_PROBE_LEN or length > MAX_PROBE_LEN:
                log(f"  FLAG: {row['assay_id']} probe length {length} "
                    f"(expected {MIN_PROBE_LEN}-{MAX_PROBE_LEN}). Verify from paper.")

    log()

    # ---- CHECK 4: each assay has forward + reverse -------------------
    log("-" * 70)
    log("CHECK 4  primer pairs complete")
    log("-" * 70)
    assays = df["assay_id"].str.replace(r"_(F|R|P)$", "", regex=True)
    df["assay_root"] = assays
    for root, group in df.groupby("assay_root"):
        types = set(group["oligo_type"])
        if "forward" not in types:
            log(f"  FLAG: {root} has no forward primer")
        if "reverse" not in types:
            log(f"  FLAG: {root} has no reverse primer")
    log("  (Probes are optional for conventional PCR)")
    log()

    # ---- CHECK 5: degenerate base count ------------------------------
    log("-" * 70)
    log("CHECK 5  degenerate base summary")
    log("-" * 70)
    log("  This tells you how many concrete sequences each oligo can")
    log("  represent. Very high counts (>1000) suggest over-degenerate")
    log("  design — flag but do not automatically exclude.")
    log()
    for _, row in df.iterrows():
        seq = str(row["sequence_5to3"]).upper().strip()
        n_degen = sum(1 for b in seq if b not in "ACGT")
        try:
            n_expand = len(expand_degenerate(seq))
        except MemoryError:
            n_expand = "TOO_MANY"
        if n_degen > 0 or (isinstance(n_expand, int) and n_expand > 100):
            log(f"  {row['assay_id']}: {n_degen} degenerate positions, "
                f"{n_expand} expansions")
    log()

    # ---- CHECK 6: sanity — locate each oligo in its design genome ----
    log("-" * 70)
    log("CHECK 6  sanity — locate each oligo in its design genome")
    log("-" * 70)
    log("  CRITICAL: if an oligo cannot be found in the genome it was")
    log("  designed against, your sequence or your search is wrong.")
    log("  Fix this before running Step 3.")
    log()

    if not REFERENCE_FASTA.exists():
        log("  ! curated_sequences.fasta not found. Skipping this check.")
        log("    Run step1 first, or set REFERENCE_FASTA to the correct path.")
    else:
        # Load all 77 genomes for the sanity search
        genomes = {rec.id.split(".")[0]: str(rec.seq).upper()
                   for rec in SeqIO.parse(str(REFERENCE_FASTA), "fasta")}
        log(f"  Genomes loaded for sanity check: {len(genomes)}")
        log()

        # Use the BDBV reference NC_014373.1 if present, otherwise any 2007 genome
        ref_id = next(
            (k for k in genomes if "NC_014373" in k),
            next((k for k in genomes), None)
        )
        if not ref_id:
            log("  ! No genome available for sanity check.")
        else:
            log(f"  Sanity genome: {ref_id}")
            log()
            ref_seq = genomes[ref_id]
            all_ok = True
            for _, row in df.iterrows():
                seq   = str(row["sequence_5to3"]).upper().strip()
                otype = str(row["oligo_type"])
                is_rev = otype == "reverse"
                hits = locate_in_genome(seq, ref_seq, is_reverse=is_rev)
                if not hits:
                    # Also try the other strand
                    hits_alt = locate_in_genome(seq, ref_seq, is_reverse=not is_rev)
                    if hits_alt:
                        log(f"  NOTE: {row['assay_id']} found on opposite strand "
                            f"at {hits_alt} — check oligo_type assignment")
                    else:
                        log(f"  ! NOT FOUND: {row['assay_id']} — sequence or "
                            f"oligo_type may be wrong. Check paper.")
                    all_ok = False
                else:
                    log(f"  OK  {row['assay_id']}: found at position(s) {hits}")
            if all_ok:
                log()
                log("  All oligos located. Sanity check passed.")
    log()

    # ---- Summary table -----------------------------------------------
    summary = (df.groupby(["assay_root", "target_gene", "organism_scope",
                            "year_designed", "publication_doi"])["oligo_type"]
               .apply(lambda x: "/".join(sorted(set(x))))
               .reset_index()
               .rename(columns={"oligo_type": "components"}))
    log("=" * 70)
    log("ASSAY PANEL SUMMARY")
    log("=" * 70)
    log(summary.to_string(index=False))
    log()
    log(f"Total assays: {len(summary)}")
    log(f"Total oligos: {len(df)}")

    df.to_csv(OUT_DIR / "validated_panel.csv", index=False)
    (OUT_DIR / "step2_report.txt").write_text("\n".join(report))

    log()
    log("=" * 70)
    log(f"WROTE  {OUT_DIR/'validated_panel.csv'}")
    log("=" * 70)
    log()
    log("NEXT: Add the BDBV reference genome (NC_014373.1) to your sequences")
    log("and proceed to Step 3: MAFFT multiple sequence alignment.")


if __name__ == "__main__":
    main()
