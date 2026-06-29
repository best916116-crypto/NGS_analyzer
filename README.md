# Amplicon Analyzer

Amplicon Analyzer is a static web app for amplicon FASTQ analysis. It runs read QC, paired-end overlap joining, reference alignment, mutation calling, and figure export from a browser page.

Open the app:

```text
https://best916116-crypto.github.io/NGS_analyzer/
```

## Features

- FASTQ, FASTQ.GZ, and joined FASTQ text input
- File or folder selection for numbered sample batches
- Automatic R1/R2 sample grouping and overlap-consensus joining
- Custom, Cas / CRISPR, TALE / TALEN, and Prime Editor target annotation
- Optional multi-site analysis
- Setting CSV import and template export
- QC table, mutation table, allele spectrum, selected-sample allele reports, summary JSON, and heatmap export
- Zoomable heatmap with full-figure and visible-viewport PNG export

## Basic Use

1. Open the app.
2. Enter the site name and full reference amplicon sequence.
3. Select the assay type and exact target positions.
4. Add FASTQ files or choose a folder.
5. Adjust QC and alignment settings if needed. Leave `Read limit per sample` blank, or enter `0`, to process all reads.
6. Run the analysis.
7. Export tables, summary files, and figures.

## Inputs

| Input | Notes |
| --- | --- |
| `.fastq`, `.fq` | Plain FASTQ |
| `.fastq.gz`, `.fq.gz` | Gzip-compressed FASTQ |
| `.fastjoin`, `.fastqjoin`, `.fqjoin`, `.join`, `.txt` | Joined FASTQ-style text |
| Site name | Label used in figures, reports, and site-specific filenames |
| Reference amplicon | Full expected PCR amplicon sequence |
| Setting CSV | Optional run configuration |

For raw paired-end files, filenames with common `R1` / `R2` markers are grouped by sample. R2 reads are reverse-complemented and joined to R1 when a high-confidence overlap is found. Unjoined pairs are reported in the QC export.

If joined FASTQ files and raw R1/R2 files are provided for the same sample, joined files are used and raw files are skipped for that sample.

`Read limit per sample` is optional. When it is blank or `0`, every parsed read is processed. Use a numeric limit only for quick validation runs because downstream read counts and allele frequencies will be based on that subset.

## Setting CSV

The app can export a setting template with one row per amplicon site:

```text
site_name,reference_amplicon_sequence,target_positions,assay_type,spacer_sequence,pam_sequence,spacer_window,tale_left_sequence,tale_right_sequence,tale_spacer_sequence,tale_padding,prime_edited_amplicon_sequence,prime_intended_edits,prime_ignore_regions,expected_edit,sample_start,sample_end,expected_sample_count,preferred_input,notes
```

Blank fields are ignored. `assay_type` may be `custom`, `cas`, `tale`, or `prime`. If `assay_type` is blank, the app infers the assay from filled spacer, PAM, TALE binding-site, or Prime Editor fields. For TALE / TALEN rows, the matched spacer is annotated as a targetable region; exact target bases come from `target_positions`.

Prime Editor mode separates intended and unwanted edits. Define intended edits with a full `prime_edited_amplicon_sequence` or a compact `prime_intended_edits` list such as `C45T,52insA,del61`. Use `prime_ignore_regions` for PBS, RTT, primer, or other regions that should remain in raw tables but be excluded from unwanted-edit signal.

## Outputs

| Output | Contents |
| --- | --- |
| `processed_mutation_table.csv` | Position-level mutation counts, intended/unwanted edit counts, and edit rates |
| `target_window_table.csv` | Target-window subset of the mutation table, including Prime Editor intended/unwanted columns |
| `substitution_matrix.csv` | Substitution counts by reference and observed base |
| `allele_spectrum_table.csv` | Unique allele sequences and edit signatures |
| `allele_report_<site>_<sample>.svg` | Publication-style allele spectrum report for the selected sample |
| `allele_spectrum_<site>_<sample>.csv` | Allele spectrum CSV for the selected sample |
| `summary.mutations.csv` | Legacy-compatible C/G mutation outcome summary by sample |
| `summary.nonX_per_pos.mutations.csv` | Legacy-compatible per-position non-reference mutation counts |
| `qc_read_counts.csv` | Read count, filtering, alignment, and paired-end join summary |
| `amplicon_run_summary.json` | Run settings, site summaries, QC totals, and top edits |
| `amplicon_mutation_heatmap.svg` | Editable vector heatmap |
| `amplicon_mutation_heatmap.png` | Full heatmap PNG |
| `amplicon_mutation_heatmap_visible.png` | Current heatmap viewport PNG |
| `amplicon_run_report.md` | Markdown run report |
| `amplicon_run_report.html` | Self-contained HTML report |

## Validation

Run the public FASTQ benchmark:

```text
node tools/amplicon-analyzer-benchmark.mjs benchmarks/public/crispresso2_base_editor/config.json
```

Run the paired-end truth test:

```text
node tools/validate-paired-end.mjs
```

The paired-end test generates matching raw R1/R2 and joined FASTQ fixtures, then checks that both analysis paths report the same aligned read count, target edit position, and edit percentage.

## Repository Layout

```text
.
|-- index.html
|-- docs/
|   `-- user_guide.html
|-- benchmarks/
|   `-- public/
|       `-- crispresso2_base_editor/
|-- tools/
|   |-- amplicon-analyzer-benchmark.mjs
|   `-- validate-paired-end.mjs
`-- README.md
```
