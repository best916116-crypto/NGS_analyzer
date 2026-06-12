# AmpliconScope NGS Analyzer

AmpliconScope is a browser-only amplicon editing NGS viewer and publication-grade figure exporter.

Open the app:

```text
https://best916116-crypto.github.io/NGS_analyzer/
```

## What It Does

- Loads preprocessed amplicon allele CSV files in the browser
- Recomputes mutation events from allele rows
- Shows QC summaries, mutation heatmaps, allele tables, and validation checks
- Exports SVG, PNG, CSV, and Markdown figure reports
- Keeps scientific labels in exported figures in English

## Input Files

Required:

- `all_mutation_raw.csv`
- `read_counts.csv`

Optional:

- `summary.nonX_per_pos.mutations.csv`
- `summary.mutations.csv`
- `sample_sheet.csv`

## Privacy

Files are parsed locally in your browser. AmpliconScope does not upload sequencing data to a server.

Do not upload private experiment data to this repository. Use your local files directly in the browser.

## Scientific Notes

The app does not treat the `ratio` column as the only source of truth. It recalculates both:

- total-aligned percentage
- no-indel percentage

Default heatmap denominator is no-indel aligned reads. Coordinates are 1-based.
