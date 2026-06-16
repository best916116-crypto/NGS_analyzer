# Amplicon Analyzer

Amplicon Analyzer is a single-page amplicon sequencing analyzer for FASTQ files, focused on rapid genome-editing QC, processed tables, and report-ready heatmap exports.

Open the app:

```text
https://best916116-crypto.github.io/NGS_analyzer/
```

## What You Can Do

- Drag and drop `.fastq`, `.fq`, `.fastq.gz`, `.fq.gz`, joined FASTQ text files, or choose a folder
- Import Amplicon Analyzer setting CSV files, or legacy setting CSV files that define `aseq`, `site`, `opts`, and `run_align_mutations`
- Paste a reference amplicon sequence and target positions
- Run read QC, unique-read aggregation, reference alignment, and mutation calling
- Track expected base-editor conversions such as `C>T` or `A>G`
- Review sample-level QC, read depth, alignment-pass counts, alignment identity, mutation heatmap, reference sequence, and allele spectrum
- Pan, zoom, fit, and save the currently visible heatmap viewport as a PNG
- Edit displayed sample labels without changing the original FASTQ files
- Export processed data, run summary, and figures as CSV, JSON, SVG, PNG, Markdown, and self-contained HTML report files

## Recommended Workflow

1. Open the app URL in a browser.
2. Paste or confirm the reference amplicon sequence.
3. Choose the assay type: `Custom window`, `Cas / CRISPR`, or `TALE / TALEN`.
4. Optional: add extra site rows in the multi-site manifest.
5. Optional: download the setting CSV template, fill it, then import it to configure the reference, assay, target window, and sample range automatically.
6. Drag in FASTQ files, choose files, or choose a folder.
7. Optional: set a first/last sample number to process only numbered sample groups. For 96-well batches, keep `Expected sample count` at `96`.
8. Review the preflight panel.
9. Click `Run analysis`.
10. Review QC, heatmap, processed table, and allele spectrum.
11. Export the files needed for downstream analysis or slides.

## Input Files

Required:

| File | Notes |
| --- | --- |
| `.fastq` / `.fq` | Plain FASTQ files |
| `.fastq.gz` / `.fq.gz` | Gzip-compressed FASTQ files in browsers that support `DecompressionStream` |
| `.fastjoin` / `.fastqjoin` / `.fqjoin` / `.join` / `.txt` | Joined or concatenated FASTQ text files |
| reference amplicon sequence | A/C/G/T/N bases pasted into the app |
| setting CSV | Optional Amplicon Analyzer template CSV, or legacy CSV with `arg`, `opts`, and `run_align_mutations` rows |
| additional site manifest | Optional rows for multi-site analysis |

The file picker accepts individual files, and the folder picker accepts a directory of many FASTQ-like files. Folder upload uses browser-local file access; files are still processed locally and are not uploaded to a server.

The reference should be the full expected PCR amplicon sequence used for alignment, not only the edited bases. The targeting window is a smaller region inside that reference that is highlighted in the figure and used as the expected editing region.

Paired-end files are grouped by sample name when filenames contain common `R1` / `R2` markers. If raw R1/R2 FASTQ files are supplied without joined FASTQ files, Amplicon Analyzer reverse-complements R2, finds a high-confidence overlap, and builds a consensus read before QC and alignment. Pairs that do not meet the overlap criteria are counted in QC as `unjoined_pairs`.

For large numbered batches, use `First sample number` and `Last sample number` to process only matching inferred sample groups. Files without a number are skipped while this filter is active.

For Illumina-style filenames such as `97_S97_L001_R1_001.fastq.gz`, the displayed sample is normalized to `97`. Joined files such as `97.fastqjoin` use the same sample name. If joined and raw R1/R2 files are mixed for the same sample, the joined files are used and raw files are skipped to avoid double counting.

## 96-Sample Batch Workflow

Amplicon Analyzer is optimized for numbered 96-sample batches:

1. Import the legacy setting CSV.
2. Choose the `fastqjoin` folder or the raw FASTQ.GZ folder. Raw R1/R2 files are joined automatically when a high-confidence overlap is found.
3. Set `First sample number` and `Last sample number` if needed.
4. Confirm that preflight reports `96` selected sample groups.

The recommended setting CSV can be downloaded from the app with `Download template`. Legacy setting CSV import still maps `--user_region_beg_offset` and `--user_region_length` to 1-based amplicon target positions. For example, if the site starts at position 61, `--user_region_beg_offset 32 --user_region_length 55` becomes target positions `29-83`.

## Setting CSV Template

The downloadable template uses one row per amplicon site:

```text
site_name,reference_amplicon_sequence,target_positions,assay_type,spacer_sequence,pam_sequence,spacer_window,tale_left_sequence,tale_right_sequence,tale_spacer_sequence,tale_padding,expected_edit,sample_start,sample_end,expected_sample_count,preferred_input,notes
```

Only filled fields are used. Blank fields are ignored.

Assay inference:

| Input columns | Inferred mode |
| --- | --- |
| `assay_type` set to `custom`, `cas`, or `tale` | Uses the selected mode |
| `pam_sequence` filled | Cas / CRISPR; `spacer_sequence` is required to map the guide |
| `tale_left_sequence` and `tale_right_sequence` | TALE / TALEN |
| `spacer_sequence` without `pam_sequence` | TALE spacer-only fallback |
| `target_positions` only | Custom window |

For Cas rows, `spacer_window` uses 1-based spacer coordinates and is mapped to amplicon coordinates after spacer/PAM matching. For TALE/TALEN rows, left/right binding sites are used to infer the intervening spacer. Additional rows are imported into the multi-site manifest with their inferred target positions.

## Multi-Site Analysis

Use the primary reference fields for the first site. To analyze additional amplicon sites in the same run, add rows to `Additional sites`:

```text
site_name, reference_amplicon_sequence, target_positions
site_B, ACGTACGT..., 16-23
site_C, ACGTACGT..., 12,18,22
```

FASTQ files are parsed and QC-filtered once per sample group, then each site reference is aligned and summarized separately. The heatmap displays one selected site at a time; the site selector switches the displayed heatmap and tables. CSV and report exports include a `site` column for all analyzed sites.

The `Minimum alignment identity (%)` setting is applied independently for every site. Reads that pass FASTQ QC but align below this identity threshold are excluded from mutation calling and counted as `low_identity_reads` in QC exports.

## Assay Type

| Mode | How the target window is set |
| --- | --- |
| `Custom window` | Uses the manual 1-based target positions entered in the app |
| `Cas / CRISPR` | Finds the spacer/guide sequence in the reference, checks the PAM pattern when possible, then maps the spacer window such as `4-8` to amplicon coordinates |
| `TALE / TALEN` | Finds left and right TALE binding sites and uses the spacer between them as the default targeting window; spacer-only matching is available as a fallback |

For TALE/TALEN assays, left/right binding sites are more reliable than spacer-only input because they define orientation and the exact intervening region.

## Outputs

| Output | Use |
| --- | --- |
| `processed_mutation_table.csv` | Site-aware position-level substitution, deletion, insertion, and edit-rate table |
| `target_window_table.csv` | Target-window-only subset of the processed mutation table |
| `substitution_matrix.csv` | Position-level substitution counts by reference and observed base |
| `allele_spectrum_table.csv` | Site-aware unique read sequences with counts, ratios, orientation, identity, and edit signature |
| `qc_read_counts.csv` | Site-aware sample read count and filtering summary |
| `amplicon_run_summary.json` | Run settings, sample-number filter, site summaries, QC totals, and top edited positions |
| `amplicon_mutation_heatmap.svg` | Editable vector heatmap |
| `amplicon_mutation_heatmap.png` | Full slide-ready heatmap image |
| `amplicon_mutation_heatmap_visible.png` | PNG crop matching the currently visible heatmap viewport |
| `amplicon_run_report.md` | Markdown record of run settings and sample summary |
| `amplicon_run_report.html` | Self-contained report with heatmap, QC, site summary, expected-edit summary, and methods notes |

## Calculation Notes

The edit signal shown in the heatmap can be selected in the app:

- `Substitution + indel edit rate`
- `Substitution only`
- `Indel only`

For each sample and position:

```text
edit percentage = selected signal reads at position / position-level covered reads * 100
```

Coordinates are 1-based amplicon positions. The heatmap uses a black grid for amplicon bases, red outlines and labels for target positions, and a reference sequence row below the figure. Terminal no-coverage gaps are not counted as deletion edits.

Raw paired-end R1/R2 files are first converted to overlap-consensus reads when possible. Reads are then filtered by FASTQ quality and length, then aligned in both orientations. Only reads at or above the configured minimum alignment identity are used for allele tables, position-level coverage, and mutation percentages.

The QC export includes `paired_reads`, `joined_pairs`, `unjoined_pairs`, and `unpaired_reads` so raw FASTQ-only runs can be audited against upstream join results.

When an expected base edit is selected, Amplicon Analyzer also reports target-window expected-edit reads:

```text
expected edit percentage = expected conversion reads at editable target bases / covered reads at editable target bases * 100
```

For example, with `C>T`, only target-window positions where the reference base is `C` are included in the expected-edit denominator.

## Browser Limits

This is a static browser app, so large FASTQ files are limited by browser memory and CPU. Start with the built-in sample data or a read limit such as `50,000` reads per sample, then increase the limit after confirming the configuration.

For production-scale primary analysis, keep the QC exports with the run and validate final calls against an established command-line workflow when changing assay design, read length, trimming, or join assumptions.

## Paired-End Validation

The repository includes a synthetic paired-end truth test:

```text
node tools/validate-paired-end.mjs
```

This creates raw R1/R2 FASTQ files and a matching joined FASTQ fixture, runs both paths through the analyzer, and checks that joined read count, aligned read count, expected edit position, and edit percentage match.

## Public FASTQ Benchmark

Amplicon Analyzer includes a reproducible public FASTQ fixture based on the CRISPResso2 EMX1 base-editor example dataset.

Source:

- CRISPResso2 examples: <https://docs.crispresso.com/latest/core/examples.html>
- FASTQ file: <https://crispresso.com/static/demo/base_editor.fastq.gz>

Run the benchmark:

```text
node tools/amplicon-analyzer-benchmark.mjs benchmarks/public/crispresso2_base_editor/config.json
```

Current benchmark result:

| Metric | Value |
| --- | ---: |
| QC-passed reads | 25,000 |
| Alignment-passed reads | 24,976 |
| Low-identity reads excluded | 24 |
| Alignment-passed unique reads | 3,915 |
| Mean alignment identity | 99.678% |
| Top target-window edit | C136, 27.709% substitution signal |
| Dominant substitution class | C>T, 12,877 substitution reads |

Generated outputs are stored in `benchmarks/public/crispresso2_base_editor/results/`, including position-level CSV, allele-spectrum CSV, QC CSV, substitution matrix, target-window table, SVG heatmap, summary JSON, and Markdown report. The app can also export a self-contained HTML report from the browser UI.

The app also includes a `Load public benchmark` button. It loads the same FASTQ fixture, reference amplicon, Cas spacer, PAM, and base-editor window into the browser UI so the benchmark can be reproduced without preparing files manually.

## Local Batch Runner

The CLI benchmark runner also supports local batch configs with a setting CSV and a folder of numbered FASTQ-like files:

```json
{
  "runName": "mito_ND1 96-sample local run",
  "settingCsv": "path/to/mito_ND1_setting.csv",
  "dataset": {
    "folder": "path/to/fastqjoin",
    "preferredInput": "joined",
    "sampleStart": 97,
    "sampleEnd": 192,
    "expectedSampleCount": 96
  },
  "settings": {
    "readLimit": 50000,
    "minQ": 20,
    "minLen": 30,
    "minIdentity": 0.8,
    "dropN": true,
    "signalMode": "all"
  },
  "outputs": {
    "resultDir": "examples/local/mito_ND1_96/results"
  }
}
```

`examples/local/` is ignored by git so local raw-data-derived outputs are not committed accidentally.

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
