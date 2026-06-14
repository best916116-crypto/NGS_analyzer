# AmpliconScope NGS Analyzer

AmpliconScope is a browser-only viewer for preprocessed amplicon NGS CSV files.

Open the app:

```text
https://best916116-crypto.github.io/NGS_analyzer/
```

The app runs in the browser. CSV files are read locally on your device and are not uploaded to a server.

## What You Can Do

- Load amplicon mutation CSV files by drag and drop
- Review read depth and basic QC checks
- View a mutation heatmap across amplicon positions
- Inspect allele-level reads for each sample
- Edit displayed sample labels in the app without changing CSV files
- Export SVG, PNG, CSV, and Markdown summary files
- Use compact multi-panel heatmaps when many samples are loaded

## Input Files

Required:

| File | What it contains |
| --- | --- |
| `all_mutation_raw.csv` | Allele-level rows with reference, aligned read, read count, ratio, and sample filename |
| `read_counts.csv` | Per-sample aligned read counts |

Optional:

| File | What it adds |
| --- | --- |
| `summary.nonX_per_pos.mutations.csv` | Position-level mutation counts for cross-checking and faster heatmap generation |
| `summary.mutations.csv` | Summary mutation counts by intended base or category |
| `sample_sheet.csv` | Cleaner sample labels, target names, condition names, and annotation ranges |

## Expected Columns

`all_mutation_raw.csv`

```text
ref,alignment,read,n_reads,ratio,filename
```

`read_counts.csv`

```text
filename,n_aligned_noindel,n_aligned_total
```

`summary.nonX_per_pos.mutations.csv`

```text
nonX_per_pos,A1,C2,G3,...,n_mut_any,n_aligned_noindel,n_aligned_total
```

The position columns are parsed as reference-base plus 1-based position, such as `G16` or `C19`.

## Basic Workflow

1. Open the app in a browser.
2. Drag in `all_mutation_raw.csv` and `read_counts.csv`.
3. Add `summary.nonX_per_pos.mutations.csv` if available.
4. Check sample count, read depth, mutation heatmap, and allele spectrum.
5. Choose the read count basis:
   - `Indel 제외 정렬 reads`: uses `n_aligned_noindel`
   - `전체 정렬 reads`: uses `n_aligned_total`
6. Set target positions, minimum reads, and figure layout.
7. Open the `Samples` tab if the displayed sample names need to be changed.
8. Export the figure or tables you need.

## Outputs

| Output | Use |
| --- | --- |
| `Figure SVG` | Editable vector heatmap |
| `Figure PNG` | Image export of the current heatmap |
| `Heatmap CSV` | Position-level percentages for downstream plotting |
| `Figure Report` | Markdown record of settings and summary values |

## Calculation Notes

Mutation percentage is calculated from read counts.

For each sample and position:

```text
mutation percentage = mutation count at position / selected read count basis * 100
```

The default read count basis is `n_aligned_noindel`.

The app does not use the input `ratio` column as the only source of truth. It recalculates displayed percentages from the available count columns.

Coordinates are 1-based amplicon positions.

## Privacy

AmpliconScope is static HTML, CSS, and JavaScript.

- No account login
- No server upload
- No backend database
- No installation step
- No raw sequencing files stored in this repository

Use private experiment files directly in your browser and do not commit them to this repository.

## Repository Layout

```text
.
├── index.html
├── demo/
│   └── ngs_amplicon_analyzer_ko.html
├── docs/
│   └── ngs_amplicon_analyzer_guide_ko.html
└── README.md
```

## Notes

This tool expects preprocessed CSV files. It is not a FASTQ aligner and does not perform primary sequencing alignment.

It is intended for visualization, QC review, and export of amplicon editing results.
