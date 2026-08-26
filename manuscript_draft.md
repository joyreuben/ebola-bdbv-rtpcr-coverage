# In silico inclusivity and exclusivity analysis of five filovirus RT-PCR assays against Bundibugyo virus genomes from the 2007, 2012, and 2026 outbreaks

*Draft — Methods and Results only. Format not yet fixed to a specific journal.*

---

## Methods

### Study design

This is an in silico inclusivity and exclusivity analysis: a computational comparison of published RT-PCR primer/probe sequences against a curated panel of Bundibugyo virus (BDBV) genomes, with no accompanying wet-laboratory work. The analysis follows the general logic of the AOAC International Appendix Q framework and the FDA's molecular diagnostic assay evaluation template, and reports assay design information consistent with MIQE guidelines. In silico mismatch analysis predicts amplification risk; it does not, by itself, establish PCR failure or success, and all findings below are reported as predictions requiring wet-laboratory confirmation.

### Genome dataset

BDBV genome sequences with an open data-use designation were obtained from Pathoplexus (`ebola-bdbv` collection, downloaded 2026-08-22), spanning three outbreaks: the founder 2007 outbreak (Bundibugyo District, Uganda), the 2012 outbreak (Democratic Republic of the Congo), and the ongoing 2026 outbreak (DRC, Uganda, and one Germany-labeled record; the epidemiological link, e.g. travel history, was not available in the metadata and is not assumed here).

Of 78 downloaded records, 8 were excluded prior to analysis:

- One 133 bp fragment (accession LQ278392.1) flagged at initial curation as an automated INSDC ingest rather than a patient-derived sequence.
- Seven further records identified during downstream verification as GenBank Patent-division deposits rather than clinical or outbreak isolates: six 2,031 bp fragments (accessions HC918458.1, HC918461.1, HV753272.1, HV753275.1, HZ164395.1, HZ164398.1; two identical triplets) and one 173 bp fragment (LQ278396.1), all submitted under patent-agent names rather than a research group, with no collection date or country of origin. Patent-division status was confirmed against NCBI's GenBank accession-prefix registry (HC = EBI Patent; HV, HZ = DDBJ Patent/JPO; LQ = EBI Patent). All 7 records mismatched at least one assay primer at a level inconsistent with every genuine BDBV genome in the dataset (see Results), which is what prompted this check.

This second round of exclusions is reported explicitly because it materially changed a downstream result (see Results) and because it illustrates a general curation risk for public outbreak-genome pulls: short, undated, single-submitter records with a non-standard two-letter accession prefix should be checked against the GenBank patent-division prefix list before inclusion, not screened on length or ambiguity alone.

The final curated dataset comprises **70 genomes**: 12 from 2007, 23 from 2012, and 35 from 2026. Reference genomes for the other three orthoebolavirus/orthomarburgvirus species used for primer sanity-checking were taken from NCBI RefSeq: Zaire ebolavirus (NC_002549.1), Sudan ebolavirus (NC_006432.1), Marburg marburgvirus (NC_001608.3), and Bundibugyo ebolavirus itself (NC_014373.1, the species reference, also present in the 70-genome working set as a 2007 sample).

### Assay panel

Five published RT-PCR assays targeting the L or NP genes were assembled (22 oligonucleotides total): Panning et al. 2007 (pan-filovirus TaqMan, L gene, a primer mix of 3 forward and 2 reverse primers with 3 species-specific probes); Ogawa et al. 2010 (EBOV-specific conventional PCR, NP gene; this assay was previously mislabeled "Nakayama 2010" in earlier internal drafts of this analysis under a non-resolving DOI — the correct source is Ogawa H, Miyamoto H, Ebihara H, Ito K, Morikawa S, Feldmann H, Takada A. *J Virol Methods* 2010;171(1):310–313, doi:10.1016/j.jviromet.2010.11.010); Sanchez et al. 1999 (pan-filovirus conventional PCR, L gene, as re-reported in Ogawa et al. 2010; its reverse primer is identical to, and shared with, the Panning FiloB primer and is entered once); Weidmann et al. 2004 (three species-typing TaqMan assays on the NP gene: ENZ for Zaire ebolavirus, ENS for Sudan ebolavirus, MN for Marburg marburgvirus); and Greninger et al. 2024 (broad pan-filovirus RT-qPCR, L gene).

Panning's assay is used as a primer mix in the original protocol (3 forward + 2 reverse primers combined in one reaction); it is scored here by best-of-mix logic (the best-matching forward primer and the best-matching reverse primer each contribute, not a fixed single pair), never as individual-primer results. Inosine positions in four oligonucleotides (Panning FAMEBOSu, FAMEBOg probes; Greninger PfiloL_F, PfiloL_R) were encoded as N and excluded from match/mismatch scoring (neither counted as a match nor a mismatch). The Weidmann MN assay and the Panning FAMMBG probe are Marburg-specific and are included deliberately as negative/specificity controls: they are expected to mismatch BDBV heavily, and a high mismatch count for these two components is a correctness check, not a coverage concern.

As a transcription sanity check, every oligonucleotide was aligned against its own design-species reference genome (or, for pan-filovirus oligos, against whichever of the four references fit best). 21 of 22 fit their design species at ≤2 mismatches. The one exception, Panning FiloB-Ra (3 mismatches against its best reference, Marburg), was confirmed correct against the original publication: it is a deliberate mixed-primer-pool member not designed to match any single species' reference alone, and its behavior in the pool is expected.

### Mismatch-tolerant alignment

Each oligonucleotide was slid against both strands of every genome using a custom aligner with full IUPAC degenerate-base awareness (`step4_oligo_alignment.py`). For each oligonucleotide/genome pair, the best-scoring placement (strand and position minimizing total mismatches, with ties broken by fewest 3′-terminal mismatches) was retained, recording: total mismatches, mismatches within the terminal 5 nucleotides of the primer's 3′ end (the region most critical to polymerase extension), and the number of non-scoring positions (inosine-derived N in the oligo, or an ambiguous/N base in the genome — genome-side ambiguity is treated as not-evaluable at that position, never as a mismatch). The aligner was validated against a synthetic-sequence self-test suite covering exact matches, single/degenerate mismatches, inosine handling, reverse-complement detection, and 3′-end mismatch placement on both strands before being run on real data.

### Per-assay coverage classification

Each multi-component assay (primer mix and/or probe) was scored by its weakest required component per genome (`step5_combine_assays.py`), using the following decision rule, applied in order:

```
worst_skip_fraction > 0.15        → NOT EVALUABLE  (binding site too ambiguous to call)
worst 3′ mismatches ≥ 2            → LIKELY FAILURE
worst 3′ mismatches = 1 and
  worst total mismatches ≥ 3       → AT RISK
worst total mismatches ≥ 6         → LIKELY FAILURE
worst total mismatches ≥ 3         → REDUCED
otherwise                          → LIKELY DETECTED
```

These thresholds are decision aids for grouping results, not proof of wet-laboratory outcome; the underlying mismatch counts are reported alongside every classification.

### Thermodynamic modeling

For genomes classified REDUCED, predicted duplex melting temperature (Tm) was modeled with `primer3-py`'s `calc_heterodimer()` (`step6_thermodynamics.py`), comparing the intact primer/reference duplex Tm against the mismatched primer/genome duplex Tm at the same fixed binding-site alignment identified in the mismatch-tolerant aligner (degenerate IUPAC positions resolved per-position against the actual genome base before modeling). Genomes classified LIKELY FAILURE were **not** Tm-modeled: at high mismatch counts, `calc_heterodimer()` stops scoring the known, fixed binding site and instead searches for whichever short, unrelated local sub-sequence is thermodynamically most favorable, which can return physically nonsensical results (verified directly: an 8-mismatch case returned a duplex Tm of −25°C, an implied "89°C shift," an artifact of the search rather than real biology). Mismatch count and 3′-position are already the basis for a LIKELY FAILURE call, so nothing is lost by reporting mismatch counts alone at that tier.

### Phylogenetic analysis

The 70 curated genomes were aligned with MAFFT v7.525 (`--auto`) and a maximum-likelihood tree was built with FastTree v2.2.0 (`-nt -gtr`). The tree was rendered with Biopython's `Phylo` module and annotated per-tip with outbreak arm and per-genome mismatch count for the one assay component showing genuine per-genome variation (Ogawa 2010 EBOV, 0–11 mismatches), to visually assess whether elevated mismatch clusters within a particular clade or outbreak arm.

---

## Results

### Dataset composition and panel integrity

The final dataset (70 genomes: 12/2007, 23/2012, 35/2026) is described in Methods above, including the 8 exclusions. Panel transcription was confirmed correct for all 22 oligonucleotides (21/22 fit their design-species reference at ≤2 mismatches; the one flagged exception was confirmed as a deliberate mixed-pool design, not a transcription error), so all mismatches reported below against BDBV reflect real sequence divergence rather than a data-entry artifact.

### Coverage by assay

**Table 1. Per-assay coverage classification by outbreak arm** (n = total genomes in that arm; mean mm = mean worst-component mismatch count among evaluable genomes only — NOT EVALUABLE genomes have no real mismatch count and are excluded from this mean, not averaged in as zero)

| Assay | Arm | n | mean mm (evaluable only) | LIKELY DETECTED | REDUCED | LIKELY FAILURE | NOT EVALUABLE |
|---|---|---|---|---|---|---|---|
| Greninger 2024 (pan-filovirus L) | 2007 | 12 | 0.00 | 12 | 0 | 0 | 0 |
| | 2012 | 23 | 0.00 | 17 | 0 | 0 | 6 |
| | 2026 | 35 | 0.00 | 25 | 0 | 0 | 10 |
| Ogawa 2010 (EBOV NP) | 2007 | 12 | 2.75 | 11 | 0 | 1 | 0 |
| | 2012 | 23 | 2.00 | 17 | 0 | 0 | 6 |
| | 2026 | 35 | 2.00 | 25 | 0 | 0 | 10 |
| Panning 2007 (pan-filovirus L) | 2007 | 12 | 4.00 | 0 | 12 | 0 | 0 |
| | 2012 | 23 | 4.00 | 0 | 17 | 0 | 6 |
| | 2026 | 35 | 4.00 | 0 | 25 | 0 | 10 |
| Sanchez 1999 (pan-filovirus L) | 2007 | 12 | 4.00 | 0 | 12 | 0 | 0 |
| | 2012 | 23 | 4.00 | 0 | 17 | 0 | 6 |
| | 2026 | 35 | 4.00 | 0 | 25 | 0 | 10 |
| Weidmann ENZ (EBOV-Zaire NP)* | 2007 | 12 | 8.00 | 0 | 0 | 12 | 0 |
| | 2012 | 23 | 8.00 | 0 | 0 | 17 | 6 |
| | 2026 | 35 | 8.00 | 0 | 0 | 23 | 12 |
| Weidmann ENS (EBOV-Sudan NP)* | 2007 | 12 | 9.08 | 0 | 0 | 12 | 0 |
| | 2012 | 23 | 9.00 | 0 | 0 | 17 | 6 |
| | 2026 | 35 | 9.00 | 0 | 0 | 22 | 13 |
| Weidmann MN (Marburg NP, control)* | 2007 | 12 | 11.00 | 0 | 0 | 12 | 0 |
| | 2012 | 23 | 11.00 | 0 | 0 | 17 | 6 |
| | 2026 | 35 | 11.00 | 0 | 0 | 22 | 13 |
| Panning MARV probe (control)* | 2007 | 12 | 9.00 | 0 | 0 | 12 | 0 |
| | 2012 | 23 | 9.00 | 0 | 0 | 17 | 6 |
| | 2026 | 35 | 9.00 | 0 | 0 | 23 | 12 |

*Species-typing/specificity controls not designed to detect BDBV; high mismatch here is the expected, correct result (see below).

An earlier version of this table reported a lower, arm-varying mean (e.g. 2.86–4.00 for Panning/Sanchez across arms) by averaging worst_mm over *all* genomes in the arm, including NOT EVALUABLE ones — whose worst_mm is a placeholder 0 from skipped positions, not a real measurement. That average was mixing a real quantity with a not-a-number placeholder and understated the true, and completely uniform, per-genome mismatch count. The values above are corrected to exclude NOT EVALUABLE genomes from the mean entirely.

Greninger 2024 shows perfect predicted amplification compatibility with essentially every evaluable BDBV genome in all three outbreaks (0 mismatches, uniformly, in every evaluable genome in every arm). Ogawa 2010 also performs well overall — a uniform 2 mismatches in every evaluable genome in 2012 and 2026, and in 11 of 12 evaluable 2007 genomes (mean 2.75 in 2007 only because of the single MS811591 outlier discussed below).

Panning 2007 and Sanchez 1999 — which share a reverse primer (FiloB) — show a **consistent, moderate reduction in predicted compatibility across all three outbreaks**: every evaluable genome in every arm falls in the REDUCED tier, at **exactly 4 mismatches with zero variance in all three arms** (2007, 2012, and 2026 alike). This is the study's central positive finding and is addressed thermodynamically below.

The three species-typing/control components (Weidmann ENZ, ENS, MN, and the Panning MARV probe) show high mismatch counts against BDBV in every arm, exactly as expected: ENZ and ENS are Zaire- and Sudan-ebolavirus-specific typing primers, and MN plus the Panning MARV probe are Marburg-specific controls. None of these were designed to detect BDBV, so their LIKELY FAILURE classification here is a specificity confirmation, not evidence of assay defect. Consistent with known filovirus phylogeny, ENZ (Zaire) sits numerically closer to BDBV — uniformly 8 mismatches in every evaluable genome in every arm — than ENS (Sudan, 9 mismatches in every evaluable genome in 2012/2026, 9–10 in 2007), a small but expected and strikingly stable signal of relative evolutionary distance.

No 2007-arm genome was classified NOT EVALUABLE for any assay, but 26–37% of genomes per assay were in the 2012 and 2026 arms, reflecting ambiguous (N) bases within that assay's binding site in the underlying sequence data rather than any property of the assay itself; this is a data-completeness limitation of the available genome records, not evidence for or against detection, and is reported separately from the four detection-outcome tiers throughout.

### Thermodynamic confirmation of the Panning/Sanchez reduction

All 108 REDUCED-tier calls for Panning 2007 and Sanchez 1999 (12 in 2007, 17 in 2012, 25 in 2026, per assay) correspond to the identical 4-mismatch pattern on the shared FiloB reverse primer — the classification-driving component in every one of these rows — with zero of those mismatches in the 3′-terminal 5 nucleotides. Thermodynamic modeling of this fixed mismatch pattern gives a reference (fully matched) duplex Tm of 61.66°C and a mismatched-duplex Tm of 41.99°C, a uniform predicted shift of **19.67°C** across all 108 rows, placing every one in the >5°C "substantial predicted destabilization" band. The other two components in the same primer/probe mix are far less affected in the same genomes (forward primer FiloA2.2: 1 mismatch, 6.88°C shift; probe FAMEBOSu: 1 mismatch, 3.41°C shift, "THERMO MARGINAL"), indicating the reduction is specifically attributable to the shared FiloB reverse primer rather than a problem with the assay's other components. This is a large, consistent predicted thermodynamic effect, not a borderline mismatch-counting artifact — though, as throughout, this is a computational prediction of reduced amplification compatibility, not a demonstration of PCR failure; wet-laboratory confirmation would be needed to establish an actual effect on assay sensitivity.

For the LIKELY FAILURE tier (Weidmann ENZ/ENS/MN, the Panning MARV probe, and the single Ogawa outlier below), mismatch counts are reported without Tm modeling, for the reasons given in Methods: mismatch count and 3′-position are already the classification-driving evidence at that tier, and primer3's duplex search becomes physically unreliable at high mismatch counts. Mismatch counts here are themselves stable and uniform across arms (e.g., Weidmann ENZ at exactly 8 mismatches in all three arms, MN at exactly 11), indicating a consistent, real species-level sequence difference rather than noise.

### A single real outlier

One genome, MS811591.1 (2007 arm, a legitimate partial BDBV sequence, GenBank direct submission, not a patent record), shows 11 mismatches against Ogawa 2010 — its only LIKELY FAILURE call for this assay in the entire dataset. This is treated as a genuine biological/sequencing outlier worth a footnote rather than a coverage concern for the assay generally.

### The reduction is a historical species-level feature, not a 2026-specific change

The identified Panning/Sanchez primer incompatibility was not unique to the 2026 lineage: it is present, at essentially the same magnitude, in the founder 2007 outbreak genomes and in the 2012 genomes (Table 1; Figure 1). A maximum-likelihood phylogeny of all 70 genomes (Figure 3) separates cleanly into three outbreak-pure clades, as expected for samples collected roughly a decade apart each; the per-genome Ogawa mismatch count, plotted alongside, does not track this clade structure — it is scattered at a similarly low level across all three clades, with the single MS811591 outlier sitting inside, not outside, the 2007 clade. Taken together, these results indicate that the primer incompatibilities identified here reflect a longstanding species-level divergence between BDBV and the EBOV/SUDV sequences these assays were designed against, already present in historical BDBV sequences, rather than a mutation newly arising in the 2026 outbreak.

### Figures

- **Figure 1** — assay × genome mismatch heatmap, all eight assay components across all 70 genomes, grouped by outbreak arm.
- **Figure 2** — per-position mismatch-frequency schematics for three representative primers (Panning/Sanchez shared FiloB; Greninger PfiloL_F; Weidmann ENZ_F), showing where along each primer divergence concentrates relative to the 3′ end.
- **Figure 3** — MAFFT/FastTree phylogeny of the 70 genomes, annotated with outbreak arm and per-genome Ogawa mismatch count.
