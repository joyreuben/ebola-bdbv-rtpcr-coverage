#!/usr/bin/env python3
"""
MIQE-required oligonucleotide table.

WHY THIS EXISTS
    MIQE guidelines require every primer/probe sequence, its target gene,
    and its source to be reported. assay_panel.csv has the sequences and
    sources; it does not have a genomic coordinate. This script adds one:
    best_start (and strand) on NC_014373.1, the BDBV RefSeq entry, which is
    itself one of the 70 curated genomes (record PP_006X6ZU, see Methods'
    note on this record's duplicate status).

    An earlier version of this script tried to report a *modal* best_start
    across all 70 genomes instead of a single reference coordinate, on the
    assumption that 70 BDBV genomes would share one coordinate frame. They
    don't: curated_metadata.csv's seq_length ranges from 133 bp to 18,942 bp
    (most genomes are partial/differently-trimmed assemblies, only a
    minority are complete full-length genomes), so a raw string-index
    best_start is not comparable across genomes and the modal position was
    shared by only 15-17/70 genomes per oligo — noise from heterogeneous
    assembly completeness, not a signal worth reporting. A single named
    reference genome is the standard, unambiguous way to report a primer's
    genomic position and is used instead.

Inputs : assay_panel.csv
         step4_output/oligo_vs_genome_full.csv
Output : step7_output/miqe_oligo_table.csv
"""
from pathlib import Path
import pandas as pd

REF_GENOME_ID = "PP_006X6ZU"  # NC_014373.1, the BDBV RefSeq entry

panel = pd.read_csv("assay_panel.csv", dtype=str)
panel.columns = panel.columns.str.strip()
full = pd.read_csv("step4_output/oligo_vs_genome_full.csv")
ref = full[full["genome_id"] == REF_GENOME_ID].set_index("assay_id")

rows = []
for _, r in panel.iterrows():
    oid = r["assay_id"]
    if oid not in ref.index:
        pos_str, strand_str, mm_str = "n/a", "n/a", "n/a"
    else:
        rr = ref.loc[oid]
        pos_str = f"{int(rr['best_start']) + 1}-{int(rr['best_start']) + int(rr['oligo_length'])}"
        strand_str = rr["strand"]
        mm_str = str(int(rr["mismatches"]))

    rows.append({
        "oligo_id": oid,
        "assay": r["notes"].split(".")[0].split(" probe")[0].split(" diagnostic")[0],
        "gene": r["target_gene"],
        "type": r["oligo_type"],
        "target_species_scope": r["organism_scope"],
        "sequence_5to3": r["sequence_5to3"],
        "length_nt": len(r["sequence_5to3"]),
        "source_doi": r["publication_doi"],
        "year_designed": r["year_designed"],
        "position_on_NC_014373.1_1based": pos_str,
        "strand": strand_str,
        "mismatches_vs_NC_014373.1": mm_str,
    })

out = pd.DataFrame(rows)
Path("step7_output").mkdir(exist_ok=True)
out.to_csv("step7_output/miqe_oligo_table.csv", index=False)
print(out.to_string(index=False))
print(f"\nwrote step7_output/miqe_oligo_table.csv ({len(out)} rows)")
