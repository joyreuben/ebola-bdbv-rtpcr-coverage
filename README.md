# In silico inclusivity/exclusivity analysis of five filovirus RT-PCR assays vs. Bundibugyo virus

An in silico inclusivity and exclusivity analysis asking whether five published filovirus
RT-PCR assays still match Bundibugyo virus (BDBV) genomes from the 2007, 2012, and 2026
outbreaks — including the ongoing 2026 DRC/Uganda outbreak. No laboratory work; purely
computational. Follows the AOAC International Appendix Q framework, the FDA molecular
diagnostic assay evaluation template, and reports assay design information consistent with
MIQE guidelines.


## Headline finding

The Panning/Sanchez shared FiloB primer mismatch against BDBV (4 mismatches, ~19.7°C
predicted duplex Tm shift) is **not unique to the 2026 lineage** — it is already present,
uniformly, in the founder 2007 outbreak genomes and in 2012. This points to a species-level
gap between the assay design (EBOV/SUDV-based) and BDBV, rather than new-lineage escape in
the current outbreak. All findings are in silico predictions of reduced amplification
compatibility, not proof of assay failure — wet-lab confirmation is the necessary final step.

## Dataset

70 curated BDBV genomes (Pathoplexus `ebola-bdbv` collection, OPEN data-use terms,
downloaded 2026-08-22): 12 from the 2007 outbreak (Uganda), 23 from 2012 (DRC), 35 from
2026 (DRC, Uganda, Germany). 8 records were excluded as GenBank Patent-division deposits
or automated ingest artifacts rather than clinical isolates — see [CLAUDE.md](CLAUDE.md#step-1--dataset-complete-patent-record-exclusion-applied-2026-08-25)
for the full exclusion audit trail.

## Assay panel

Five published RT-PCR assays (22 oligonucleotides total), targeting the L or NP gene:

| Assay | Gene | Type | DOI |
|---|---|---|---|
| Panning 2007 | L | pan-filovirus TaqMan (primer mix) | 10.1086/520539 |
| Ogawa 2010 | NP | EBOV-specific conventional | 10.1016/j.jviromet.2010.11.010 |
| Sanchez 1999 | L | pan-filovirus conventional | via Ogawa 2010 |
| Weidmann 2004 | NP | species-typing TaqMan (Zaire/Sudan/Marburg) | 10.1016/j.jcv.2003.09.004 |
| Greninger 2024 | L | pan-filovirus broad RT-qPCR | 10.1128/spectrum.01010-24 |

See [assay_panel.csv](assay_panel.csv) for the full oligo table.

## Pipeline

| Step | Script | Purpose |
|---|---|---|
| 1 | `step1_build_dataset.py` | Load, verify, curate genomes |
| 2 | `step2_validate_panel.py` | Validate assay panel CSV |
| 4 | `step4_oligo_alignment.py` | Mismatch-tolerant oligo-vs-genome alignment, IUPAC-aware |
| 5 | `step5_combine_assays.py` | Mixed-primer per-assay coverage classification |
| 6 | `step6_thermodynamics.py` | primer3-py Tm-shift modeling (REDUCED-tier calls) |
| 7 | `step7_output/fig1_heatmap.py`, `fig2_binding_schematics.py`, `fig3_phylogeny.py` | Figures |

(There is no Step 3; numbering follows the project's internal step log in [CLAUDE.md](CLAUDE.md).)

Each script's output lives in a matching `stepN_output/` directory. Final figures are in
[figures/](figures/) as matched PNG/PDF pairs.

## Reproducing

```bash
conda activate amr-multiomics
pip install pandas biopython primer3-py matplotlib seaborn --break-system-packages
conda install -n amr-multiomics -c bioconda -c conda-forge mafft fasttree -y

python step1_build_dataset.py
python step2_validate_panel.py
python step4_oligo_alignment.py
python step5_combine_assays.py
python step6_thermodynamics.py
python step7_output/fig1_heatmap.py
python step7_output/fig2_binding_schematics.py
python step7_output/fig3_phylogeny.py
```

## Repository layout

```
step1_output/   curated_sequences.fasta, curated_metadata.csv
step2_output/   validated_panel.csv
step4_output/   oligo_vs_genome_full.csv (per-oligo mismatch alignment)
step5_output/   assay_coverage_per_genome.csv, assay_coverage_summary.csv, step5_report.txt
step6_output/   thermodynamics_per_component.csv, thermodynamics_summary.csv, step6_report.txt
step7_output/   figure-generating scripts + MAFFT/FastTree intermediates
figures/        final figures (fig1 heatmap, fig2 binding schematics, fig3 phylogeny)
references/     NCBI RefSeq genomes for EBOV/SUDV/MARV/BDBV
```

## Interpreting the results

In silico mismatch and predicted thermodynamic destabilization indicate *risk* of reduced
amplification efficiency — they do not establish PCR failure. All results in this repository
should be read as candidates for wet-lab confirmation, not as standalone diagnostic claims.

## Data source

Genome sequences: [Pathoplexus](https://pathoplexus.org/) `ebola-bdbv` collection, records
with OPEN data-use terms only. Reference genomes: NCBI RefSeq.

## Contact

Joy Okebugwu, MSc Molecular Biology and Genomics.
