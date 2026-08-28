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

### Step 1 — Dataset (COMPLETE, patent-record exclusion applied 2026-08-25)
Source: Pathoplexus ebola-bdbv, Data use terms = OPEN, downloaded 2026-08-22.

**Files:**
- `step1_output/curated_sequences.fasta` — 70 genomes
- `step1_output/curated_metadata.csv` — 78 rows (8 flagged excluded)

**Dataset breakdown (evaluable genomes):**
- Bdbv-2026: 35 genomes (DRC, Uganda, Germany)
- Bdbv-2012: 23 genomes (DRC)
- Bdbv-2007: 12 genomes (Uganda)

**Exclusions and fixes documented:**
- PP_006X66F.1 excluded: 133 bp fragment (LQ278392.1), automated INSDC ingest, not a patient genome.
- **2026-08-25 correction — 7 more patent-division GenBank records excluded from the Bdbv-2007 arm**, found while verifying a downstream Greninger result:
  - PP_006X5S9.1 / HC918458.1, PP_006X5T7.1 / HC918461.1, PP_006X5U5.1 / HV753272.1, PP_006X5V3.1 / HV753275.1, PP_006X5W1.1 / HZ164395.1, PP_006X5XZ.1 / HZ164398.1 — six 2031 bp fragments, single submitter "Steigerwald, R.", no country/date.
  - PP_006X67D.1 / LQ278396.1 — 173 bp fragment, same submitter ("Vanatta, R.") as the already-excluded PP_006X66F.1, missed the first time.
  - **Why excluded:** confirmed against NCBI's official accession-prefix table that `HC` = EBI Patent division, `HV`/`HZ` = DDBJ Patent (JPO), `LQ` = EBI Patent — these are GenBank patent-office sequence deposits, not clinical/outbreak isolates. Same exclusion basis as PP_006X66F.1, just missed on the first pass. All 8 excluded rows are kept in metadata with `excluded` column flagged (audit trail), sequences removed from the FASTA.
  - **Why it mattered:** every one of these 8 records was inflating or fabricating apparent BDBV divergence for at least one assay (see Step 4/5/6 below — the entire "Greninger fails in the 2007 arm" result was 3 copies of one patent sequence, not real biology). **Any future re-download of this dataset must re-screen for GenBank Patent-division prefixes (2-letter codes — check https://www.ncbi.nlm.nih.gov/genbank/acc_prefix/), not just short/fragment length.**
- 10 Bdbv-2007 records had no country metadata. Assigned Uganda (2007 outbreak was geographically confined to Bundibugyo District, Uganda).
- Identical-sequence duplicate groups: of the 20 rows originally described as "retained as independent patient samples," 14 are genuinely legitimate (full-genome accessions from Towner et al. 2007 and Hulseberg et al. 2012, biologically plausible clonality within a confined outbreak) — the other 6 were the patent-triplicate records above and are now excluded, not merely deduplicated.
- 6 high-ambiguity genomes retained; N positions in binding sites are scored NOT EVALUABLE at Step 5, never as a match.
- Most records carry INSDC accessions; use for deduplication if pulling from GenBank.
- **2026-08-26 — full prefix audit completed (previously only re-checked the original 8).** A reviewer pointed out that the "must re-screen for GenBank Patent-division prefixes" instruction above had never actually been run against the full curated 70, only re-confirmed on the 8 already excluded — a narrative claim, not a completed check. Running it found 2 more: **JA489018.1** (PP_006X5YX.1, 2007 arm) and **MS811591.1** (PP_006X6B5.1, 2007 arm) both carry EBI Patent-division prefixes. Neither was excluded — the original 8 were excluded on prefix *plus* a patent-agent-style single submitter *plus* content divergent from every other record, and these two fail that combination: JA489018.1 is submitted by Towner JS et al. (the real CDC group, not a patent agent) and is byte-identical (MD5-matched) to NC_014373.1/FJ217161.1 — genuine BDBV data, already covered by the duplicate-triplet note above. MS811591.1 is submitted by Manuguerra/Escadafal/Kwasiborski/Vanhomwegen, a named Institut Pasteur hemorrhagic-fever-diagnostics group also linked to patent WO2016189490A1 and a companion Sci Rep paper — a substantial (6,633 bp) partial genome with unique sequence content, not a short single-submitter fragment. **This determination is provisional, not confirmed with the submitting groups.** It also means the manuscript's earlier description of MS811591.1 as "not a patent record" was wrong and has been corrected to "does not match the exclusion pattern used for the other 8." **Lesson: "check the prefix list" is not itself a completed action — it has to actually be run against the full retained set, not just re-verified on records already known to be excluded, or the instruction just restates itself without adding information.**
- **2026-08-26 — "54 independently sequenced genomes" was an overclaim.** Of the 54 genomes classified REDUCED for Panning/Sanchez, 6 identical-hash duplicate groups (including the JA489018.1/FJ217161.1/NC_014373.1 triplet above) account for 14 of those rows — only **46 unique underlying sequences**, not 54. The manuscript's thermodynamic-confirmation section previously cited "54 independently sequenced genomes" as replication evidence without checking this; corrected throughout to 46.
  - **Two different unique-sequence counts exist — do not mix them up.** **46** = unique sequences within the **54-genome Panning/Sanchez REDUCED subset** (Results, thermodynamic section). **62** = unique sequences across **all 70 curated records** (6 duplicate-hash groups covering 14 rows; Abstract/Methods). A draft abstract briefly described the full 70 as "representing 46 unique sequences," which is wrong — 46 is the REDUCED-subset figure. Whenever quoting a unique-sequence count, state which set it is over.

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
| Panning 2007 | L | pan-filovirus TaqMan | 3F + 2R + 3 probes (primer mix) | 10.1086/520600 |
| Ogawa 2010 | NP | EBOV conventional | 1F + 1R | 10.1016/j.jviromet.2010.11.010 |
| Sanchez 1999 | L | pan-filovirus conventional | 1F + 1R (shares FiloB with Panning) | via Ogawa |
| Weidmann 2004 | NP | species-typing TaqMan | 3 assays: ENZ (Zaire), ENS (Sudan), MN (Marburg) | 10.1016/j.jcv.2003.09.004 |
| Greninger 2024 | L | pan-filovirus broad | 1F + 1R | 10.1128/spectrum.01010-24 |

**Citation correction (2026-08-25):** the NP-gene assay (primers "FiloNP-Fe"/"FiloNP-Re")
was previously mislabeled "Nakayama 2010" with a non-resolving DOI
(10.1128/JCM.02641-09). Verified correct source: **Ogawa H, Miyamoto H,
Ebihara H, Ito K, Morikawa S, Feldmann H, Takada A. "Detection of all known
filovirus species by reverse transcription-polymerase chain reaction using a
primer set specific for the viral nucleoprotein gene." J Virol Methods.
2011;171(1):310-313 (epub 2010 Nov 19). DOI: 10.1016/j.jviromet.2010.11.010.** Renamed
`NAKAYAMA2010_EBOV_F/R` → `OGAWA2010_EBOV_F/R` in `assay_panel.csv`,
`step2_output/validated_panel.csv`, and the `ASSAY_STRUCTURE` dict in
`step5_combine_assays.py`. No real "Nakayama 2010 JCM filovirus RT-PCR"
paper could be found to exist — do not reintroduce that name.

**Citation correction (2026-08-26) — Panning 2007 DOI was wrong.** All 8
Panning oligo rows in `assay_panel.csv` and `step2_output/validated_panel.csv`
carried DOI `10.1086/520539`. That DOI does not resolve to the Panning
filovirus RT-PCR paper — it resolves to a one-page front-matter dedication
("Dedication: Jim Orzechowski (1944–2003) and Michael Kiley (1942–2004)",
*J Infect Dis* 196(Suppl 2):S127–S128) in the same 2007 supplement issue.
Verified correct source: **Panning M, Laue T, Ölschläger S, Eickmann M,
Becker S, Raith S, et al. "Diagnostic reverse-transcription polymerase chain
reaction kit for filoviruses based on the strain collections of all European
biosafety level 4 laboratories." J Infect Dis. 2007;196(Suppl 2):S199–S204.
DOI: 10.1086/520600.** Fixed in `assay_panel.csv`, `step2_output/validated_panel.csv`,
and this file. Not caught earlier because the wrong DOI still resolved to a
real, plausibly-adjacent-looking article (same journal, same year, same
supplement) rather than a dead link — always spot-check that a DOI resolves
to the *named* paper, not just that it resolves to *something*.

**Panning 2007 / BDBV reconciliation (2026-08-26, researched on request):**
Panning et al. 2007 was published 2007-11-15 — the same week Uganda's MoH was
first notified of the outbreak that would turn out to be BDBV (2007-11-28)
— and its validation panel (Zaire, Sudan, Côte d'Ivoire, Reston ebolavirus;
Marburg virus, several strains each) does not include Bundibugyo virus,
which did not exist as a described taxon yet. **There is no published
empirical Panning-vs-BDBV result to reconcile with the predicted REDUCED
call** — the assay predates the virus. Independently, the actual 2007
outbreak diagnostic history (Towner et al. 2008, PLoS Pathog 4(11):e1000212,
doi:10.1371/journal.ppat.1000212) is corroborating, not conflicting: species-
specific real-time RT-PCR assays for Zaire/Sudan ebolavirus and Marburg virus
were initially **negative** on real BDBV specimens; detection succeeded only
via "a traditionally less sensitive... but more broadly reactive filovirus L
gene-specific RT-PCR assay," followed by sequencing that showed the virus was
a new, divergent species. That is exactly this project's core pattern
(species-typing NP assays fail, broad L-gene pan-filovirus assays retain
reduced-but-real detection capability) playing out in the field in 2007, from
an independent source. Worth a Discussion paragraph citing Towner 2008
directly — strengthens the paper, doesn't need reconciling against it.

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

**Output:** `step4_output/oligo_vs_genome_full.csv` — 1,540 rows (22 oligos × 70 genomes,
rerun 2026-08-25 after excluding the 8 patent-division records)
Columns: assay_id, genome_id, mismatches, three_prime_mismatches, skipped_positions,
oligo_length, strand, best_start.

**Key per-oligo mean mismatches (all 70 genomes, post patent-exclusion):**
- Greninger PfiloL_F: **0.00**, PfiloL_R: **0.00** — was 0.45/0.53 when the 8 patent
  records were still in the dataset. The entire non-zero signal on this assay came
  from those records; every real BDBV genome scores a perfect match.
- Ogawa EBOV_F: 1.66, EBOV_R: 0.91
- Panning FiloA2.2: 0.77, FiloA2.3: 2.31, FiloA2.4: 2.31, FiloB: 3.09, FiloB-Ra: 5.40
- Weidmann ENZ_F: 3.64, ENZ_R: 3.87, ENZ_P: 6.13
- Weidmann ENS_F: 6.10, ENS_R: 6.91, ENS_P: 6.56
- Weidmann MN_F: 4.96, MN_R: 6.11, MN_P: 8.37
- Panning FAMMBG (Marburg probe): 6.89

---

### Step 4b — Target-site coverage detection (ADDED 2026-08-28)
**Script:** `step4b_target_coverage.py`

**The bug this fixes.** Step 4's aligner always returns a best-scoring window.
It has no representation of "this oligo's binding site is not in this record
at all." Handed a partial genome that does not span a binding site, it slides
the oligo across whatever sequence IS present and reports the least-bad
window. The resulting mismatch count measures unrelated sequence.

**It bit exactly once, and the manuscript had already written it up as
biology.** MS811591.1 (PP_006X6B5, 2007 arm) is a 6,633 bp partial genome
covering only NC_014373.1 **11,567–18,199**. The Ogawa primers target
**1,028–1,621** — ~10 kb outside the record. Step 4 returned 11 mismatches;
Step 5 called LIKELY FAILURE; the draft called it "a single real outlier...
a genuine biological/sequencing outlier." It was an artifact. 9 of 1,540
oligo/genome pairs were affected, all in this one record (Ogawa F/R, Panning
FAMMBG, Weidmann ENS F/R/P, ENZ F/R, MN F).

**How it works.** Locate each oligo's target on NC_014373.1 with Step 4's own
`best_alignment()` (so coordinates agree with the MIQE table by construction);
establish each genome's reference coverage by exact-matching 40 nt anchors
sampled every 100 nt (BDBV is ~99% identical to the BDBV reference so anchors
land densely; anchors hitting the reference twice are discarded as ambiguous);
flag an oligo/genome pair when the target falls outside the genome's mapped span.

**The span test is deliberately NOT interval-based.** An interval-based version
also fired on *internal* anchor gaps, which are not missing sequence — they are
long N-runs. Verified: the internal gaps in PP_006X6QC / PP_006X6JQ are N-runs
at 11,829–13,044, which Step 5's `skip_fraction` already scores at 0.91–0.96,
far above its 0.15 cutoff. Calling those "absent" would be redundant; calling a
merely *divergent* internal region absent would discard a real measurement.
Span-based isolates the one thing skip_fraction cannot see: a record that does
not reach the coordinate at all.

**Flagged pairs are dropped BEFORE the best-of-mix pick in Step 5**, not
filtered after — a phantom count can be spuriously low and would otherwise be
selected as a primer mix's "best" member over a genuine higher-scoring one.
`best_component()` takes `absent_pairs`; a genome missing from a role's frame
means every member of that role is absent → NOT EVALUABLE.

Outputs: `step4_output/target_presence.csv`, `genome_reference_coverage.csv`,
`step4b_report.txt`. New Step 5 column: `target_absent`.

**Lesson (same shape as the prefix-audit lesson above):** the manuscript's
Supplementary Table 2 already stated this exact confound — "an elevated
mismatch is also consistent with the true binding region simply not being
present in the retrieved excerpt" — about the *excluded* records, and never
checked it against the *retained* 70. Writing a caveat is not running a check.
**Any new partial/fragmentary record must be coverage-checked, not just
length-checked.**

---

### Step 5 — Mixed-primer combining (COMPLETE, phantom-10 fix applied 2026-08-25, coverage-aware since 2026-08-28)
**Script:** `step5_combine_assays.py`

Scores each assay by its weakest required component, using best-of-mix logic for Panning.
Outputs: `step5_output/assay_coverage_per_genome.csv`, `step5_output/assay_coverage_summary.csv`

**Classification thresholds (declared in methods):**
```
worst_skip_fraction > 0.15       → NOT EVALUABLE   (binding site too N-rich to call)
worst_3p >= 2                    → LIKELY FAILURE
worst_3p == 1 AND worst_mm >= 3  → AT RISK
worst_mm >= 6                    → LIKELY FAILURE
worst_mm >= 3                    → REDUCED
else                             → LIKELY DETECTED
```

**2026 arm results (n=35, unaffected by the patent-record exclusion — all excluded
records were in the 2007 arm):**
- Greninger PfiloL: 25 LIKELY DETECTED, 10 NOT EVALUABLE
- Ogawa EBOV: 25 LIKELY DETECTED, 10 NOT EVALUABLE
- Panning panfilo: 25 REDUCED, 10 NOT EVALUABLE
- Sanchez FiloAB: 25 REDUCED, 10 NOT EVALUABLE
- Weidmann ENZ: 23 LIKELY FAILURE, 12 NOT EVALUABLE
- Weidmann ENS: 22 LIKELY FAILURE, 13 NOT EVALUABLE
- Weidmann MN (Marburg control): 0 LIKELY DETECTED, 22 LIKELY FAILURE, 13 NOT EVALUABLE ✓ phantom-10 gone

**2007 arm results (n=12, post patent-exclusion — this arm changed completely):**
- Greninger PfiloL: **12 LIKELY DETECTED (100%)** — was 12 LIKELY DETECTED + 3 REDUCED
  + 3 LIKELY FAILURE + 1 more LIKELY FAILURE (LQ278396.1) before exclusion. Every
  non-DETECTED result in this arm came from the 8 removed patent records.
- Ogawa EBOV: 11 LIKELY DETECTED, 1 LIKELY FAILURE (MS811591.1, 11 mm — a real
  6633 bp partial genome, legitimate GenBank direct submission, genuine outlier
  worth a footnote, not another artifact)
- Panning panfilo: 12 REDUCED (100%), mean 4 mismatches, uniform
- Sanchez FiloAB: 12 REDUCED (100%), mean 4 mismatches, uniform
- Weidmann ENZ: 12 LIKELY FAILURE (100%)
- Weidmann ENS: 12 LIKELY FAILURE (100%)
- Weidmann MN (Marburg control): 12 LIKELY FAILURE (100%) ✓ correct specificity

Full per-arm breakdown (2007/2012/2026 x all 8 assay components) in `step5_output/step5_report.txt`.

**Note:** NOT EVALUABLE now accounts for a meaningful share of every assay/arm
(29% of the 2026 arm for several assays). This is genome N-content, not assay
performance — worth a sentence in methods/limitations so reviewers don't read
it as a coverage gap in the assay itself.

**Correctness fix 2026-08-25 (caught while drafting the manuscript):**
`assay_coverage_summary.csv`'s `mean_worst_mm` column previously averaged
`worst_mm` over *all* genomes in an arm, including NOT EVALUABLE ones — whose
`worst_mm` is a placeholder 0 from skipped positions, not a real measurement.
This silently understated every mean and made it look like arms differed
(e.g. Panning/Sanchez appeared to run 2.86–4.00 across arms) when in fact
every evaluable genome in every arm has the exact same mismatch count. Fixed
in `step5_combine_assays.py` to exclude NOT EVALUABLE rows from the mean;
rerun and confirmed. **Corrected means (evaluable genomes only, identical
across 2007/2012/2026 unless noted):** Greninger 0.00; Ogawa 2.00 (2.75 in
2007 only, due to the single MS811591 outlier); Panning panfilo 4.00; Sanchez
FiloAB 4.00; Weidmann ENZ 8.00; Weidmann ENS 9.00 (9.08 in 2007); Weidmann MN
11.00; Panning MARV control 9.00 — all with zero variance except the one
named exception. This makes the "historical, not new" finding *stronger*,
not weaker: the mismatch level isn't just similar across arms, it is
identical.

---

### Step 6 — Thermodynamics (COMPLETE, scope revised 2026-08-25)
**Script:** `step6_thermodynamics.py`

Originally scoped (per the old plan below) to run `primer3-py` Tm-shift modeling
on all AT RISK oligos — but AT RISK is empty post phantom-10 fix (Step 5's
thresholds route everything into REDUCED or LIKELY FAILURE instead). Rescoped
to model **REDUCED-tier rows only** (111 of 616 assay-genome rows). This was
a correctness finding, not a scope shortcut — see below.

**Why LIKELY FAILURE rows are NOT Tm-modeled:** `primer3.calc_heterodimer()`
does not score the fixed, known binding site from Step 4 — it *searches* for
whichever local partial alignment is thermodynamically most favorable,
tolerating large unpaired dangling ends. At low mismatch counts (REDUCED,
3–5 mm) the full intended site stays the optimum, so the numbers are sound.
At high mismatch counts (LIKELY FAILURE, 6+ mm) primer3 instead finds a
short, mostly-unrelated stable sub-stretch and reports *its* Tm — verified
by hand on WEIDMANN2004_ENZ_P vs one BDBV genome: 8 mismatches produced a
duplex Tm of **−25°C** (an "89°C shift"), which is a search artifact, not
biology. primer3 has no mode to force the known fixed-frame alignment.
Mismatch count/3′-position is already decisive at the LIKELY FAILURE tier
regardless — that's exactly why Step 5 calls it LIKELY FAILURE — so nothing
is lost by reporting mismatch counts alone there instead of a misleading Tm
number. Full technical detail in the script's docstring.

Degenerate IUPAC bases in several primers (Ogawa, Greninger, the Panning/
Weidmann probes) are resolved per-position against the actual genome base
using Step 4's own match/mismatch/skip rule before Tm modeling (see
`realize_duplex()`), so mismatch counts always agree with Step 4/5.

**REDUCED-tier results, rerun 2026-08-25 after excluding the 8 patent-division
records (worst/classification-driving component, Tm shift in °C):**
- Panning panfilo / Sanchez FiloAB (share FiloB reverse primer): 4 mismatches
  → 19.67°C shift, **all 108 REDUCED rows identical**, now spanning all three
  arms uniformly (2007 n=12, 2012 n=17, 2026 n=25) — the same divergence
  pattern in the founder 2007 outbreak as in 2012 and 2026, which is exactly
  the "historical, not new" evidence this project needs.

**The Greninger PfiloL_F "2007 arm, n=3, 36.35°C shift" result from the
2026-08-24 run no longer exists and must not appear in the manuscript.**
It was never real biology: it was three copies of a single patent-office
sequence (GenBank Patent division, prefixes HC/HV/HZ — see Step 1). With
those records excluded, Greninger PfiloL scores 0 mismatches against every
BDBV genome in the dataset (all three arms) and is 100% LIKELY DETECTED
wherever evaluable. This is also independently supported by the Greninger
2024 paper's own wet-lab validation: their tested BDBV isolate was
NGS-identified (Supp. Table S3) as best-matching NC_014373.1, which is
genome `PP_006X6ZU.1` in this dataset — and that genome scores 0 mismatches
against PfiloL_F/R in this pipeline too. The assay's own validation isolate
and this project's real BDBV genomes agree; only the patent artifacts disagreed.

All REDUCED-tier rows land in the >5°C "THERMO CONFIRMS" band — i.e. even
the borderline mismatch-count calls correspond to a substantial predicted
duplex destabilization, not a marginal one.

Outputs: `step6_output/thermodynamics_per_component.csv` (REDUCED tier, all
components), `step6_output/thermodynamics_summary.csv` (worst-component
summary), `step6_output/likely_failure_mismatch_summary.csv` (LIKELY FAILURE
mismatch counts, no Tm), `step6_output/step6_report.txt`.

**Step 6b — pinned-alignment reanalysis (added 2026-08-26, in response to
review):** `calc_heterodimer()`'s free-search behavior (above) was flagged
as a gap even at the REDUCED tier — never independently hand-verified there
the way it was at LIKELY FAILURE. `step6_output/step6b_nn_reanalysis.py`
closes it using `Bio.SeqUtils.MeltingTemp.Tm_NN` with an explicit `c_seq`
(no new dependency — Biopython is already required). Given an explicit
`c_seq`, `Tm_NN` cannot search; it scores nearest-neighbor thermodynamics,
with internal-mismatch parameters (Allawi & SantaLucia 1997), for exactly
the alignment it's handed. Reused the same realized primer/template pair
Step 6 already builds, so both methods score the identical fixed alignment.
**Result: 22.80°C shift (ref 55.62°C, duplex 32.82°C), vs primer3's
19.67°C — larger, not smaller, uniform across all 108 REDUCED rows for the
same reason primer3's number was (one primer, one fixed mismatch pattern).**
Absolute Tm differs between methods (different default salt/dNTP/oligo-
concentration parameters and NN tables) but the shift, which is what Step 5/6
actually uses, agrees in direction and magnitude. This is real independent
confirmation, not circular: `Tm_NN` structurally cannot reproduce the
local-optimum-search artifact that broke `calc_heterodimer()` at high
mismatch counts, because it has no search step at all. c_seq orientation
convention (non-obvious, verified empirically before trusting it): `Tm_NN`
wants `c_seq[i]` = the base pairing with `seq[i]`, in the SAME left-to-right
order as `seq` — i.e. `c_seq = complement(seq)`, NOT `reverse_complement(seq)`
(the latter raises `ValueError: no thermodynamic data for neighbors`). Since
`build_template()` already returns the template 5'->3' (= `reverse_complement(primer)`
in the perfect-match case), the correct `c_seq` is `template[::-1]` (reversal
is a pure string op, holds in the mismatch case too).
Outputs: `step6_output/nn_reanalysis_per_component.csv`,
`step6_output/nn_vs_primer3_comparison.csv`, `step6_output/step6b_report.txt`.

---

## What comes next (in order)

1. ~~Figures~~ — COMPLETE 2026-08-25, see Step 7 below.
2. Write up results

---

### Step 7 — Figures (COMPLETE 2026-08-25)
**Scripts:** `step7_output/fig1_heatmap.py`, `fig2_binding_schematics.py`,
`fig3_phylogeny.py`. **Outputs:** `figures/fig1_mismatch_heatmap.{png,pdf}`,
`fig2_binding_site_schematics.{png,pdf}`, `fig3_phylogeny.{png,pdf}`.

New dependencies installed into the `amr-multiomics` conda env: **mafft**
7.525 and **fasttree** 2.2.0 (via bioconda/conda-forge — needed for Fig 3
only; Figs 1-2 use only matplotlib/pandas/biopython, already present).

- **Fig 1 — mismatch heatmap.** 8 assay calls × 70 genomes (from
  `step5_output/assay_coverage_per_genome.csv`, worst-component mismatch
  count), columns grouped by outbreak arm, NOT EVALUABLE cells rendered as
  a distinct grey rather than folded into the colour scale. Confirms the
  whole-project story at a glance: Greninger and Ogawa mostly light,
  Panning/Sanchez a uniform mid-blue band across all three arms, Weidmann +
  MARV control uniformly dark (expected).
- **Fig 2 — binding-site schematics.** Per-position mismatch frequency
  (recomputed fresh from Step 4's own `best_alignment()`, not read from a
  CSV — no prior step kept per-position detail) for three primers chosen to
  contrast the coverage story: Panning/Sanchez shared FiloB (the real
  REDUCED-tier gap — 3 internal hotspots, none in the 3′ window, consistent
  with REDUCED rather than LIKELY FAILURE), Greninger PfiloL_F (clean),
  Weidmann ENZ_F (intentional off-target control, hotspots include the 3′
  window).
- **Fig 3 — phylogeny.** `step1_output/curated_sequences.fasta` (70
  genomes) aligned with MAFFT (`--auto`, ~1 min), tree built with FastTree
  (`-nt -gtr`, ~40s). Rendered with Bio.Phylo + matplotlib (headless, no
  ete3/Xvfb needed). Two tip strips: outbreak arm (categorical) and Ogawa
  2010 EBOV NP worst-mismatch count (the one assay with real per-genome
  texture, 0-11 mm — Panning/Sanchez FiloB was considered but is a flat 4 mm
  across every evaluable genome, so it has no texture to show on a tree).
  Tree cleanly separates into three arm-pure clades (2012 / 2007 / 2026), as
  expected for genomes sampled at different times — but the Ogawa mismatch
  colouring does **not** track those clades; it is scattered uniformly
  across all three, with the single MS811591 outlier (11 mm) visible
  sitting inside the 2007 clade rather than on some separate branch. This
  is the visual counterpart of the "species-level gap, not new-lineage
  escape" claim.
- Regenerate any figure after a Step 4/5/6 rerun by just re-running its
  script — none of the three modify pipeline outputs, they only read them
  (Fig 3 additionally reads/regenerates `step7_output/aligned.fasta` and
  `tree.nwk`, which do not need to change unless the genome set itself changes).

---

## Key interpretation notes (do not lose these)

- The temporal signal: **the identified primer incompatibilities (Panning/Sanchez
  FiloB, 4 mm, 19.67°C shift) were not unique to the 2026 lineage — they were
  already present in the founder 2007 outbreak genomes and in 2012, uniformly.**
  This is a species-level gap (BDBV vs EBOV/SUDV), not a new 2026 variant escape.
  Say it this way, not as a "not new variant escape" double-negative — reviewers
  will trust you more for the positive, specific claim. (Confirmed stronger than
  before: pre-correction, the 2007 arm's own signal was partly a patent-sequence
  artifact; post-correction, the real 2007 genomes show the identical uniform
  4 mm pattern seen in 2012/2026 — see Step 6.)
- Weidmann ENZ (Zaire, mean ~5.6 in 2026) sits closer to BDBV than ENS (Sudan, mean ~6.3). Fits known filovirus phylogeny. Worth a sentence.
- All Weidmann high mismatches against BDBV are expected — species-typing primers for a third species. Report as specificity confirmation, not failure, exactly as already framed above.
- In silico mismatch predicts risk, it does not prove assay failure. State this plainly every time you make a claim, using language like "predicted reduced amplification compatibility," "predicted thermodynamic destabilisation," or "potential reduction in amplification efficiency" rather than "PCR failure." Wet-lab confirmation is the final word.
- primer3's `calc_heterodimer()` cannot be trusted for Tm at high mismatch counts — it searches for the best local partial alignment rather than scoring a fixed known site, and finds nonsense (e.g. negative Tm) once mismatches are numerous. Only use it in the REDUCED-mismatch regime. Full explanation in `step6_thermodynamics.py`'s docstring.
- **GenBank Patent-division records are a real contamination risk in outbreak
  sequence pulls, not just short/fragment length.** Check 2-letter accession
  prefixes against https://www.ncbi.nlm.nih.gov/genbank/acc_prefix/ for every
  new download, not just N-content/length. See the 2026-08-25 Step 1 correction.
- Co-author needed: someone with diagnostics-lab experience, ideally through GHIDI. Also send courtesy email to INRB Kinshasa (submitting group on the Pathoplexus records) even though the data is open — they may want to be involved.

---

## Scripts in project folder

| Script | Purpose | Status |
|---|---|---|
| `step1_build_dataset.py` | Load, verify, curate genomes | Complete |
| `step2_validate_panel.py` | Validate assay panel CSV | Complete |
| `step4_oligo_alignment.py` | Mismatch-tolerant oligo vs genome | Complete |
| `step5_combine_assays.py` | Mixed-primer per-assay coverage | Complete |
| `step6_thermodynamics.py` | primer3-py Tm-shift modeling, REDUCED tier | Complete |
| `step7_output/fig1_heatmap.py` | assay x genome mismatch heatmap | Complete |
| `step7_output/fig2_binding_schematics.py` | per-position mismatch schematics | Complete |
| `step7_output/fig3_phylogeny.py` | MAFFT+FastTree phylogeny, mismatch-annotated | Complete |
| `step4b_target_coverage.py` | Target-site coverage detection (absent-site → NOT EVALUABLE) | Complete (2026-08-28) |
| `step6_output/step6b_nn_reanalysis.py` | Pinned-alignment Tm reanalysis (Bio.SeqUtils Tm_NN) | Complete (2026-08-26) |
| `step5_output/skip_fraction_sensitivity.py` | 0.10/0.15/0.20 NOT-EVALUABLE cutoff sensitivity | Complete (2026-08-26) |
| `step1_output/supp_excluded_records.py` | Refetch + classify the 8 excluded records for Supp. Table 1 | Complete (2026-08-26) |
| `step7_output/miqe_oligo_table.py` | MIQE-required oligo table (sequence/gene/source/position) | Complete (2026-08-26) |

`diagnose_phantom10.py` was never created/committed — the skip-fraction fix
was applied directly to `step5_combine_assays.py` without a separate
diagnostic run. Not needed now that Step 5 is fixed and verified.

---

## Install dependencies if needed
```
pip install pandas biopython primer3-py matplotlib seaborn --break-system-packages
conda install -n amr-multiomics -c bioconda -c conda-forge mafft fasttree -y
```
