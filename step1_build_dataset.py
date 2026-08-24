#!/usr/bin/env python3
"""
Step 1: Build and verify the open BDBV genome dataset.

Inputs (download from Pathoplexus with Data use terms = OPEN):
    - metadata TSV  ("Download all entries" -> metadata)
    - sequences FASTA ("Download all entries" -> sequences)

Outputs:
    - curated_metadata.csv   one row per genome, cleaned and checked
    - curated_sequences.fasta  sequences that passed the checks
    - step1_report.txt       everything you need for the Methods section

This script does NOT drop anything on length or completeness. Per-assay
coverage is decided later, after alignment. Here we only verify, describe
and flag.

Requires: pandas, biopython
    pip install pandas biopython --break-system-packages
"""

import hashlib
import re
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO

# ---------------------------------------------------------------- config

METADATA_FILE = Path("metadata.tsv")
FASTA_FILE = Path("sequences.fasta")
OUT_DIR = Path("step1_output")

# Recorded for the Methods section. Update every time you re-download.
DOWNLOAD_DATE = "2026-08-22"
QUERY_DESCRIPTION = "Pathoplexus ebola-bdbv, Data use terms = OPEN, all outbreaks"
SEQSET_ID = "TODO_create_a_seqset_and_paste_its_id_here"

# BDBV reference genome length (NC_014373.1), used only as a sanity anchor.
REFERENCE_LENGTH = 18940

# ------------------------------------------------------- column handling

# Pathoplexus column names can change. We look for these ideas, not exact
# strings, and report anything we cannot find rather than failing silently.
COLUMN_PATTERNS = {
    "accession": r"accession.*version|^accession$",
    "outbreak": r"outbreak",
    "collection_date": r"collection.*date",
    "country": r"^country$|geo.*country",
    "admin1": r"admin.*level.*1|division|province",
    "data_use_terms": r"data.*use.*terms",
    "submitting_group": r"submitting.*group",
    "insdc_accession": r"insdc",
    "length": r"^length$|sequence.*length",
    "authors": r"^authors$",
}


def normalise(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", str(name).strip().lower()).strip()


def map_columns(df: pd.DataFrame) -> dict:
    found, missing = {}, []
    normalised = {normalise(c): c for c in df.columns}
    for key, pattern in COLUMN_PATTERNS.items():
        match = next(
            (orig for norm, orig in normalised.items() if re.search(pattern, norm)),
            None,
        )
        if match:
            found[key] = match
        else:
            missing.append(key)
    if missing:
        print(f"  ! could not find columns for: {', '.join(missing)}")
        print(f"    available columns: {list(df.columns)}")
    return found


# ------------------------------------------------------------- sequences


def ambiguous_fraction(seq: str) -> float:
    """Fraction of positions that are not an unambiguous A/C/G/T."""
    seq = seq.upper()
    if not seq:
        return 1.0
    unambiguous = sum(seq.count(b) for b in "ACGT")
    return 1.0 - (unambiguous / len(seq))


def sequence_hash(seq: str) -> str:
    return hashlib.md5(seq.upper().encode()).hexdigest()


def load_sequences(path: Path) -> pd.DataFrame:
    rows = []
    for record in SeqIO.parse(str(path), "fasta"):
        seq = str(record.seq)
        rows.append(
            {
                "fasta_id": record.id,
                "seq_length": len(seq),
                "n_count": seq.upper().count("N"),
                "ambiguous_fraction": round(ambiguous_fraction(seq), 5),
                "seq_md5": sequence_hash(seq),
                "sequence": seq,
            }
        )
    if not rows:
        sys.exit(f"No sequences parsed from {path}")
    return pd.DataFrame(rows)


def join_key(value: str) -> str:
    """Strip version suffix so PP_007ESMY.1 matches PP_007ESMY."""
    return str(value).split()[0].split(".")[0].strip()


# ------------------------------------------------------------------ main


def main() -> None:
    for path in (METADATA_FILE, FASTA_FILE):
        if not path.exists():
            sys.exit(f"Missing input file: {path.resolve()}")

    OUT_DIR.mkdir(exist_ok=True)
    report = []

    def log(line: str = "") -> None:
        print(line)
        report.append(line)

    log("=" * 70)
    log("STEP 1: BDBV OPEN DATASET")
    log("=" * 70)
    log(f"Downloaded : {DOWNLOAD_DATE}")
    log(f"Query      : {QUERY_DESCRIPTION}")
    log(f"SeqSet     : {SEQSET_ID}")
    log()

    meta = pd.read_csv(METADATA_FILE, sep="\t", dtype=str)
    log(f"Metadata rows read: {len(meta)}")
    cols = map_columns(meta)
    log()

    seqs = load_sequences(FASTA_FILE)
    log(f"Sequences read: {len(seqs)}")
    log()

    # --- join metadata to sequences -----------------------------------
    acc_col = cols.get("accession")
    if not acc_col:
        sys.exit("Cannot proceed without an accession column.")

    meta["_key"] = meta[acc_col].map(join_key)
    seqs["_key"] = seqs["fasta_id"].map(join_key)

    merged = meta.merge(seqs, on="_key", how="outer", indicator=True)
    only_meta = (merged["_merge"] == "left_only").sum()
    only_seq = (merged["_merge"] == "right_only").sum()
    if only_meta or only_seq:
        log(f"  ! {only_meta} metadata rows without a sequence")
        log(f"  ! {only_seq} sequences without metadata")
        log("    Investigate before continuing. Do not ignore this.")
    df = merged[merged["_merge"] == "both"].drop(columns=["_merge"]).copy()
    log(f"Matched records: {len(df)}")
    log()

    # --- CHECK 1: every record really is open -------------------------
    log("-" * 70)
    log("CHECK 1  data use terms")
    log("-" * 70)
    dut_col = cols.get("data_use_terms")
    if dut_col:
        counts = df[dut_col].fillna("MISSING").value_counts()
        for value, n in counts.items():
            log(f"  {value}: {n}")
        non_open = df[~df[dut_col].astype(str).str.upper().str.contains("OPEN")]
        if len(non_open):
            log()
            log(f"  ! {len(non_open)} record(s) are NOT open. Publishing with")
            log("    these needs submitter permission. Accessions:")
            for acc in non_open[acc_col].head(20):
                log(f"      {acc}")
    else:
        log("  ! no data use terms column found. Add it in Customize columns")
        log("    and re-download. Do not assume open.")
    log()

    # --- CHECK 2: composition by outbreak and country -----------------
    log("-" * 70)
    log("CHECK 2  dataset composition")
    log("-" * 70)
    ob_col, ctry_col = cols.get("outbreak"), cols.get("country")
    if ob_col and ctry_col:
        table = pd.crosstab(
            df[ob_col].fillna("UNLABELLED"), df[ctry_col].fillna("UNKNOWN")
        )
        log(table.to_string())
        log()
        unlabelled = df[df[ob_col].isna()]
        if len(unlabelled):
            log(f"  ! {len(unlabelled)} record(s) have no outbreak label:")
            for acc in unlabelled[acc_col]:
                log(f"      {acc}")
            log("    These may still be usable. Check them individually.")
    elif ob_col:
        log(df[ob_col].fillna("UNLABELLED").value_counts().to_string())
    log()

    # --- CHECK 3: collection date spread ------------------------------
    log("-" * 70)
    log("CHECK 3  collection date spread")
    log("-" * 70)
    date_col = cols.get("collection_date")
    if date_col:
        df["_date"] = pd.to_datetime(df[date_col], errors="coerce")
        missing_dates = df["_date"].isna().sum()
        log(f"  records with unparseable or missing date: {missing_dates}")
        if ob_col:
            for outbreak, group in df.groupby(df[ob_col].fillna("UNLABELLED")):
                valid = group["_date"].dropna()
                if len(valid):
                    span = (valid.max() - valid.min()).days
                    log(
                        f"  {outbreak}: n={len(group)}, "
                        f"{valid.min().date()} to {valid.max().date()} "
                        f"({span} days)"
                    )
        log()
        log("  A narrow span is not a problem, but it is a limitation you")
        log("  must state rather than discover at review.")
    log()

    # --- CHECK 4: length and ambiguity --------------------------------
    log("-" * 70)
    log("CHECK 4  sequence quality")
    log("-" * 70)
    log(f"  reference length (NC_014373.1): {REFERENCE_LENGTH}")
    log(f"  observed lengths: min {df['seq_length'].min()}, "
        f"median {int(df['seq_length'].median())}, "
        f"max {df['seq_length'].max()}")
    log()
    log("  ambiguity distribution (non-ACGT fraction):")
    for label, bounds in [
        ("clean (<0.1%)", (0, 0.001)),
        ("low (0.1-1%)", (0.001, 0.01)),
        ("moderate (1-5%)", (0.01, 0.05)),
        ("high (>5%)", (0.05, 1.01)),
    ]:
        n = df["ambiguous_fraction"].between(bounds[0], bounds[1], "left").sum()
        log(f"    {label}: {n}")
    log()
    log("  Nothing is dropped on ambiguity here. What matters is whether an")
    log("  ambiguous base sits inside an oligo binding site, which is a")
    log("  Step 5 question. Genomes with an N in a binding site are scored")
    log("  NOT EVALUABLE for that assay, never as a match.")
    log()

    # --- CHECK 5: duplicates ------------------------------------------
    log("-" * 70)
    log("CHECK 5  duplicates")
    log("-" * 70)
    dup_seq = df[df.duplicated("seq_md5", keep=False)].sort_values("seq_md5")
    if len(dup_seq):
        log(f"  {len(dup_seq)} record(s) share an identical sequence with another.")
        log("  Identical genomes are expected within an outbreak and are NOT")
        log("  errors. Only merge them if they are the same sample submitted")
        log("  twice. Check the accessions below by hand:")
        for _, row in dup_seq.iterrows():
            log(f"      {row[acc_col]}  md5={row['seq_md5'][:8]}")
    else:
        log("  no identical sequences found")
    log()

    insdc_col = cols.get("insdc_accession")
    if insdc_col:
        with_insdc = df[insdc_col].notna().sum()
        log(f"  records carrying an INSDC accession: {with_insdc}")
        log("  Use this column to deduplicate if you also pull from GenBank.")
        log("  Never merge NCBI and Pathoplexus downloads without it.")
    log()

    # --- outputs -------------------------------------------------------
    keep = [c for c in cols.values() if c in df.columns]
    keep += ["seq_length", "n_count", "ambiguous_fraction", "seq_md5"]
    curated = df[keep].copy()
    curated.to_csv(OUT_DIR / "curated_metadata.csv", index=False)

    with open(OUT_DIR / "curated_sequences.fasta", "w") as handle:
        for _, row in df.iterrows():
            handle.write(f">{row[acc_col]}\n{row['sequence']}\n")

    log("=" * 70)
    log(f"WROTE  {OUT_DIR/'curated_metadata.csv'}  ({len(curated)} rows)")
    log(f"WROTE  {OUT_DIR/'curated_sequences.fasta'}")
    log("=" * 70)
    log()
    log("NEXT: Step 2, build the assay panel. Then align, and confirm your")
    log("2026 vs 2007 divergence lands near the published 219 nucleotide")
    log("differences. If it does not, the alignment is wrong, not the data.")

    (OUT_DIR / "step1_report.txt").write_text("\n".join(report))


if __name__ == "__main__":
    main()
