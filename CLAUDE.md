# CLAUDE.md — Ebola BDBV RT-PCR Assay Coverage Project

## Engineering principle
Whenever possible, use well-tested packages to solve the problem. Do not re-invent the wheel.

---

## Who I am
Joy Okebugwu. MSc student, Molecular Biology and Genomics, Nigeria.
Python on WSL2/Windows 11. Conda environment: `amr-multiomics`.
Python path: `/home/joy/miniconda3/envs/amr-multiomics/bin/python`
Working directory: `~/ebola project/`

---

## What this project is

An in silico inclusivity and exclusivity analysis asking whether five published
filovirus RT-PCR assays still match the Bundibugyo virus (BDBV) genomes from
the 2026 DRC/Uganda outbreak. No laboratory work. Purely computational.
Target output: journal paper (Eurosurveillance, Journal of Clinical Virology,
IJID Regions, Viruses, or Microbiology Spectrum).

Methodological framework: AOAC Appendix Q, FDA molecular diagnostic template,
MIQE guidelines. Study type has an established name: in silico inclusivity
and exclusivity analysis.

---

## Project status: what is done and what is not

### Step 1 — Dataset (COMPLETE)
Source: Pathoplexus ebola-bdbv, Data use terms = OPEN, downloaded 2026-08-22.

**Files:**
- `step1_output/curated_sequences.fasta` — 77 genomes
- `step1_output/curated_metadata.csv` — 78 rows (one flagged excluded)

**Dataset breakdown:**
- Bdbv-2026: 35 genomes (DRC, Uganda, Germany)
- Bdbv-2012: 23 genomes (DRC)
- Bdbv-2007: 19 genomes (Uganda)

**Exclusions and fixes documented:**
- PP_006X66F.1 excluded: 133 bp fragment, automated INSDC ingest, not a patient genome. Row kept in metadata with `excluded` column flagged.
- 10 Bdbv-2007 records had no country metadata. Assigned Uganda (2007 outbreak was geographically confined to Bundibugyo District, Uganda). Documented under `country_source` column as "inferred from 2007 outbreak record".
- 20 identical-sequence pairs retained as independent patient samples (expected in outbreak clusters).
- 6 high-ambiguity genomes retained; N positions in binding sites are scored NOT EVALUABLE at Step 5, never as a match.
- 70 of 78 records carry INSDC accessions; use for deduplication if pulling from GenBank.

**Species reference genomes** (fetched from NCBI, cached in `references/`):
- EBOV: NC_002549.1 (18959 bp)
- SUDV: NC_006432.1 (18875 bp)
- MARV: NC_001608.3 (19111 bp)
- BDBV: NC_014373.1 (18940 bp, the reference)

---

### Step 2 — Assay panel (COMPLETE)
**File:** `assay_panel.csv` — 22 oligos, validated by `step2_validate_panel.py`

**Five assays:**

| Assay | Gene | Type | Oligos | DOI |
|---|---|---|---|---|
| Panning 2007 | L | pan-filovirus TaqMan | 3F + 2R + 3 probes (primer mix) | 10.1086/520539 |
| Nakayama 2010 | NP | EBOV conventional | 1F + 1R | 10.1128/JCM.02641-09 |
| Sanchez 1999 | L | pan-filovirus conventional | 1F + 1R (shares FiloB with Panning) | via Nakayama |
| Weidmann 2004 | NP | species-typing TaqMan | 3 assays: ENZ (Zaire), ENS (Sudan), MN (Marburg) | 10.1016/j.jcv.2003.09.004 |
| Greninger 2024 | L | pan-filovirus broad | 1F + 1R | 10.1128/spectrum.01010-24 |

**Critical design notes:**
- Panning uses a PRIMER MIX: 3 forward primers (FiloA2.2, 2.3, 2.4) and 2 reverse primers (FiloB, FiloB-Ra) are mixed in the tube. Assay performance = best-of-forwards AND best-of-reverses. Never report individual primers as the assay result.
- Inosine (I) in Panning probes (FAMEBOSu: 2 positions, FAMEBOg: 1 position) and Greninger primers (PfiloL_F: 2 positions, PfiloL_R: 1 position) was replaced with N for computation. Inosine positions are non-scoring — they count as neither match nor mismatch.
- Marburg oligos (PANNING2007_FAMMBG and all WEIDMANN2004_MN) are SPECIFICITY CONTROLS. They should NOT match BDBV. High mismatches against BDBV for these is correct and expected.
- FiloB reverse primer is shared between Panning and Sanchez. Entered once under Panning only.
- FiloB-Ra flagged CHECK! in sanity check (3 mm vs Marburg best). Confirmed correct by hand against paper — it is a deliberate mix member not designed to match any single species alone.

---

### Step 4 — Mismatch-tolerant alignment (COMPLETE)
**Script:** `step4_oligo_alignment.py`

Replaced exact-match Check 6 with a proper mismatch-tolerant sliding window aligner.
Full IUPAC degenerate base awareness. N (inosine) positions skipped. Both strands tested.
Three-prime window: last 5 bases of each primer scored separately.

**Sanity check result:** 21/22 oligos fit their design species at ≤2 mismatches.
Panel transcription confirmed correct. All mismatches against BDBV are real biology.

**Output:** `step4_output/oligo_vs_genome_full.csv` — 1,694 rows (22 oligos × 77 genomes)
Columns: assay_id, genome_id, mismatches, three_prime_mismatches, skipped_positions,
oligo_length, strand, best_start.

**Key per-oligo mean mismatches (all 77 genomes):**
- Greninger PfiloL_F: 0.45, PfiloL_R: 0.53
- Nakayama EBOV_F: 2.47, EBOV_R: 2.03
- Panning FiloA2.2: 1.69, FiloA2.3: 3.01, FiloA2.4: 3.05, FiloB: 4.05, FiloB-Ra: 6.13
- Weidmann ENZ_F: 3.92, ENZ_R: 4.27, ENZ_P: 6.48
- Weidmann ENS_F: 6.48, ENS_R: 7.31, ENS_P: 6.90
- Weidmann MN_F: 5.26, MN_R: 6.36, MN_P: 8.73
- Panning FAMMBG (Marburg probe): 7.23

---

### Step 5 — Mixed-primer combining (RUN ONCE, FIX PENDING)
**Script:** `step5_combine_assays.py`

Scores each assay by its weakest required component, using best-of-mix logic for Panning.
Outputs: `step5_output/assay_coverage_per_genome.csv`, `step5_output/assay_coverage_summary.csv`

**Classification thresholds (declared in methods):**
```
worst_3p >= 2                    → LIKELY FAILURE
worst_3p == 1 AND worst_mm >= 3  → AT RISK
worst_mm >= 6                    → LIKELY FAILURE
worst_mm >= 3                    → REDUCED
else                             → LIKELY DETECTED
```

**2026 arm results (PENDING FIX — see red flag below):**
- Greninger: 35/35 LIKELY DETECTED ✓
- Nakayama: 35/35 LIKELY DETECTED ✓
- Panning: 10 DETECTED + 25 AT RISK
- Sanchez: 10 DETECTED + 25 AT RISK
- Weidmann ENS: 10 DETECTED + 22 AT RISK + 3 FAILURE
- Weidmann ENZ: 10 DETECTED + 23 AT RISK + 2 FAILURE
- Weidmann MN (Marburg control): 10 DETECTED + 1 AT RISK + 24 FAILURE

---

## RED FLAG — OPEN, MUST FIX BEFORE WRITING RESULTS

**The phantom-10 problem.** Every assay, including the Marburg negative control
that should never detect BDBV, shows exactly 10 LIKELY DETECTED in the 2026 arm.
A Marburg primer detecting 10 Bundibugyo genomes is wrong.

**Cause:** Some genomes have long N-runs in the binding region. N positions are
skipped rather than scored as mismatches. A heavily N-rich binding site can produce
mismatches=0 simply because most positions were unevaluable, and the current code
misclassifies this as LIKELY DETECTED.

**The fix:** Add a NOT EVALUABLE class. In `step5_combine_assays.py`, update
`best_component()` to track `skip_fraction = skipped_positions / oligo_length`,
then update `classify()` to check skip_fraction first:

```python
def classify(worst_mm, worst_3p, worst_skip_fraction):
    if worst_skip_fraction > 0.15:   # >15% of binding site unreadable
        return "NOT EVALUABLE"
    if worst_3p >= 2:
        return "LIKELY FAILURE"
    if worst_3p == 1 and worst_mm >= 3:
        return "AT RISK"
    if worst_mm >= 6:
        return "LIKELY FAILURE"
    if worst_mm >= 3:
        return "REDUCED"
    return "LIKELY DETECTED"
```

Pass `worst_skip_fraction` through from `best_component()` into the classification
call in the main loop. Rerun Step 5. Marburg control should drop to near-zero
LIKELY DETECTED. The phantom 10 should convert to NOT EVALUABLE.

**Diagnostic script:** `diagnose_phantom10.py` — run this first to confirm the
suspect genomes have higher skip_fraction than the rest before applying the fix.

---

## What comes next (in order)

1. Run `diagnose_phantom10.py` — confirm phantom-10 is skip-fraction driven
2. Apply the fix to `step5_combine_assays.py` as described above
3. Rerun Step 5 — confirm Marburg control near-zero DETECTED
4. Step 6: thermodynamics using `primer3-py` on all AT RISK oligos — predict Tm shift per mismatch position to confirm AT RISK calls before writing claims
5. Figures: heatmap (assays × genomes, coloured by mismatch), binding site schematics for worst assays, phylogeny annotated with mismatch lineages
6. Write up results

---

## Key interpretation notes (do not lose these)

- The temporal signal: 2012 and 2026 genomes show similar divergence from the assay designs. This is a species-level gap (BDBV vs EBOV/SUDV), not a new 2026 variant escape. Frame it this way — reviewers will trust you more for it.
- Weidmann ENZ (Zaire, mean 5.63) sits closer to BDBV than ENS (Sudan, mean 6.34). Fits known filovirus phylogeny. Worth a sentence.
- All Weidmann high mismatches against BDBV are expected — species-typing primers for a third species. Report as specificity confirmation, not failure.
- In silico mismatch predicts risk, it does not prove assay failure. State this plainly every time you make a claim. Wet-lab confirmation is the final word.
- Co-author needed: someone with diagnostics-lab experience, ideally through GHIDI. Also send courtesy email to INRB Kinshasa (submitting group on the Pathoplexus records) even though the data is open — they may want to be involved.

---

## Scripts in project folder

| Script | Purpose | Status |
|---|---|---|
| `step1_build_dataset.py` | Load, verify, curate genomes | Complete |
| `step2_validate_panel.py` | Validate assay panel CSV | Complete |
| `step4_oligo_alignment.py` | Mismatch-tolerant oligo vs genome | Complete |
| `step5_combine_assays.py` | Mixed-primer per-assay coverage | Needs fix |
| `diagnose_phantom10.py` | Confirm cause of phantom-10 | Ready to run |

---

## Install dependencies if needed
```
pip install pandas biopython primer3-py --break-system-packages
```
