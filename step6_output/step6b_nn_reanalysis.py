#!/usr/bin/env python3
"""
Step 6b: independent, pinned-alignment Tm reanalysis of the REDUCED tier.

WHY THIS EXISTS
    Step 6's caveat (see step6_thermodynamics.py docstring, and the
    manuscript's Methods "Thermodynamic modeling" section) is that
    primer3-py's calc_heterodimer() performs its own internal search for
    the best local alignment - it is not constrained to score Step 4's
    fixed, known binding site. At low mismatch counts the search is
    expected to recover the intended site, but this was only actually
    hand-verified to break down at HIGH mismatch counts (LIKELY FAILURE
    tier); it was never independently confirmed to hold at the REDUCED
    tier the way it should have been.

    Bio.SeqUtils.MeltingTemp.Tm_NN, given an explicit c_seq, does not
    search - it scores nearest-neighbor thermodynamics for EXACTLY the
    two aligned sequences handed to it, using Allawi & SantaLucia (1997)
    internal-mismatch parameters (the same class of pinned-alignment,
    internal-mismatch-aware nearest-neighbor method the manuscript's
    Methods section flagged as the fix for this gap, e.g. MELTING5 /
    DECIPHER::CalculateEfficiencyPCR). It requires no new dependency -
    Biopython is already a project requirement.

METHOD
    Reuses step6_thermodynamics.realize_duplex() / build_template() to
    get the exact same concrete (all-ACGT, mismatch-realized) primer and
    template Step 6 already computed for each REDUCED row's classification-
    driving component. c_seq for Tm_NN must be complement(primer), i.e.
    reverse(template) - see inline comment; verified empirically against
    Tm_NN's own no-c_seq default before use below.

    ref_tm_nn    = Tm_NN(primer, c_seq=complement(primer))   perfect duplex
    duplex_tm_nn = Tm_NN(primer, c_seq=complement_of_template) actual duplex
    tm_shift_nn  = ref_tm_nn - duplex_tm_nn

    Absolute Tm values from Tm_NN are NOT expected to match primer3's
    calc_heterodimer (different default salt/oligo-concentration
    parameters, different NN tables) - only the SHIFT is comparable in
    kind, and only the shift is what Step 6 uses to classify
    THERMO CONFIRMS / MARGINAL / MINIMAL.

Requires: pandas, biopython (already a project dependency; no new install)
Inputs  : assay_panel.csv, step1_output/curated_sequences.fasta,
          step4_output/oligo_vs_genome_full.csv,
          step5_output/assay_coverage_per_genome.csv
Outputs : step6_output/nn_reanalysis_per_component.csv
          step6_output/step6b_report.txt
"""
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import MeltingTemp as mt

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import step4_oligo_alignment as s4          # noqa: E402
import step6_thermodynamics as s6           # noqa: E402

PANEL_FILE = Path("assay_panel.csv")
GENOMES_FILE = Path("step1_output/curated_sequences.fasta")
STEP4_FILE = Path("step4_output/oligo_vs_genome_full.csv")
STEP5_FILE = Path("step5_output/assay_coverage_per_genome.csv")
PRIMER3_FILE = Path("step6_output/thermodynamics_per_component.csv")
OUT_DIR = Path("step6_output")


def nn_tm_shift(concrete_primer: str, template: str) -> dict:
    """
    c_seq convention for Bio.SeqUtils.MeltingTemp.Tm_NN: c_seq must be
    given such that c_seq[i] is the base pairing with seq[i], read in the
    SAME left-to-right order as seq. Empirically verified: for a perfect
    duplex, Tm_NN(seq, c_seq=complement(seq)) == Tm_NN(seq) (its own
    default), while Tm_NN(seq, c_seq=reverse_complement(seq)) raises a
    "no thermodynamic data for neighbors" error (wrong orientation).
    template (from build_template()) is the actual template strand
    5'->3', i.e. reverse_complement(primer) in the perfect-match case, so
    c_seq = reverse(template) = complement(primer) exactly - the correct
    Tm_NN input - and this holds position-for-position in the mismatch
    case too, since reversal is a pure string operation independent of
    base identity.
    """
    c_seq_ref = s4.rc(concrete_primer)[::-1]
    c_seq_actual = template[::-1]
    ref = mt.Tm_NN(concrete_primer, c_seq=c_seq_ref, strict=False)
    dup = mt.Tm_NN(concrete_primer, c_seq=c_seq_actual, strict=False)
    return {"ref_tm_nn": round(ref, 2), "duplex_tm_nn": round(dup, 2),
            "tm_shift_nn": round(ref - dup, 2)}


def main():
    OUT_DIR.mkdir(exist_ok=True)
    report = []

    def log(s=""):
        print(s)
        report.append(s)

    panel = pd.read_csv(PANEL_FILE, dtype=str)
    panel.columns = panel.columns.str.strip()
    primer_seq_by_id = dict(zip(panel["assay_id"], panel["sequence_5to3"].str.upper()))

    genomes = {rec.id.split(".")[0]: str(rec.seq).upper()
               for rec in SeqIO.parse(str(GENOMES_FILE), "fasta")}
    step4_idx = pd.read_csv(STEP4_FILE, dtype={"genome_id": str}).set_index(["assay_id", "genome_id"])
    step5 = pd.read_csv(STEP5_FILE, dtype={"genome_id": str})
    reduced = step5[step5["classification"] == "REDUCED"].copy()

    log("=" * 70)
    log("STEP 6b: PINNED-ALIGNMENT NEAREST-NEIGHBOR Tm REANALYSIS (Bio.SeqUtils)")
    log("REDUCED tier only, classification-driving (worst) component only")
    log("=" * 70)
    log(f"REDUCED rows: {len(reduced)}")
    log()

    rows = []
    for _, row in reduced.iterrows():
        genome_seq = genomes.get(row["genome_id"])
        role = row["worst_component"]
        oligo_id = row.get(f"chosen_{role}")
        if genome_seq is None or pd.isna(oligo_id):
            continue
        primer_seq = primer_seq_by_id.get(oligo_id)
        key = (oligo_id, row["genome_id"])
        if key not in step4_idx.index:
            continue
        a = step4_idx.loc[key]
        window = genome_seq[a["best_start"]: a["best_start"] + a["oligo_length"]]

        realized = s6.realize_duplex(primer_seq, window, a["strand"])
        concrete_primer = realized["primer"]
        template = s6.build_template(realized["genome"], a["strand"])
        nn = nn_tm_shift(concrete_primer, template)

        rows.append({
            "assay": row["assay"], "genome_id": row["genome_id"],
            "outbreak": row["outbreak"], "oligo_id": oligo_id,
            "mismatches": int(a["mismatches"]),
            **nn,
        })

    out = pd.DataFrame(rows)
    out.to_csv(OUT_DIR / "nn_reanalysis_per_component.csv", index=False)
    log(f"Wrote {len(out)} rows to nn_reanalysis_per_component.csv")
    log()

    log("-" * 70)
    log("Per-oligo summary (Bio.SeqUtils Tm_NN, pinned alignment)")
    log("-" * 70)
    summ = (out.groupby("oligo_id")
            .agg(n=("tm_shift_nn", "size"),
                 mean_shift=("tm_shift_nn", "mean"),
                 min_shift=("tm_shift_nn", "min"),
                 max_shift=("tm_shift_nn", "max"),
                 std_shift=("tm_shift_nn", "std"))
            .round(3))
    log(summ.to_string())
    log()

    if PRIMER3_FILE.exists():
        p3 = pd.read_csv(PRIMER3_FILE)
        p3_worst = p3[p3["is_worst_component"]][["assay", "genome_id", "tm_shift"]]
        merged = out.merge(p3_worst, on=["assay", "genome_id"], how="inner",
                            suffixes=("_nn", "_primer3"))
        merged = merged.rename(columns={"tm_shift": "tm_shift_primer3"})
        merged.to_csv(OUT_DIR / "nn_vs_primer3_comparison.csv", index=False)
        log("-" * 70)
        log("primer3 calc_heterodimer() shift vs Bio.SeqUtils Tm_NN pinned shift")
        log("-" * 70)
        log(f"n compared: {len(merged)}")
        log(f"primer3 shift: mean={merged['tm_shift_primer3'].mean():.2f}  "
            f"min={merged['tm_shift_primer3'].min():.2f}  max={merged['tm_shift_primer3'].max():.2f}")
        log(f"Tm_NN shift  : mean={merged['tm_shift_nn'].mean():.2f}  "
            f"min={merged['tm_shift_nn'].min():.2f}  max={merged['tm_shift_nn'].max():.2f}")
        both_confirm = (merged['tm_shift_primer3'] >= 5.0) & (merged['tm_shift_nn'] >= 5.0)
        log(f"Rows where BOTH methods exceed the 5°C 'THERMO CONFIRMS' band: "
            f"{both_confirm.sum()} / {len(merged)}")

    (OUT_DIR / "step6b_report.txt").write_text("\n".join(report))
    log()
    log(f"WROTE {OUT_DIR}/nn_reanalysis_per_component.csv")
    log(f"WROTE {OUT_DIR}/nn_vs_primer3_comparison.csv")
    log(f"WROTE {OUT_DIR}/step6b_report.txt")


if __name__ == "__main__":
    main()
