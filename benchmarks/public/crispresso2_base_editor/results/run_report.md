# CRISPResso2 EMX1 Base Editor Public FASTQ Benchmark

- Created: 2026-06-15T10:13:03.766Z
- Dataset: CRISPResso2 base-editor example
- Source: https://docs.crispresso.com/latest/core/examples.html
- FASTQ: https://crispresso.com/static/demo/base_editor.fastq.gz
- Reference length: 238 bp
- Samples: 1
- QC-passed reads: 25000
- Alignment-passed reads: 24976
- Edit signal shown in heatmap: Substitution-only edit rate
- Target source: Cas forward spacer match - PAM matched - spacer window 4-13
- Target positions: 135, 136, 137, 138, 139, 140, 141, 142, 143, 144
- Minimum alignment identity: 80.0%
- Max edit percentage: 27.709%

## QC Summary

| sample | QC-passed reads | alignment-passed reads | low-identity excluded | unique reads | mean Q | mean identity | max edit % | warnings |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| emx1_base_editor | 25000 | 24976 | 24 | 3915 | 37.29 | 99.68% | 27.709% | low-identity reads excluded |

## Target Window

| position | ref | signal % | substitution reads | expected edit reads | expected edit % |
| ---: | --- | ---: | ---: | ---: | ---: |
| 136 | C | 27.709% | 6920 | 6708 | 26.860% |
| 137 | C | 22.916% | 5723 | 5650 | 22.624% |
| 135 | T | 0.364% | 91 | 0 | 0.000% |
| 141 | C | 0.144% | 36 | 16 | 0.064% |
| 138 | G | 0.092% | 23 | 0 | 0.000% |
| 139 | A | 0.088% | 22 | 0 | 0.000% |
| 142 | A | 0.072% | 18 | 0 | 0.000% |
| 144 | A | 0.072% | 18 | 0 | 0.000% |

## Substitution Spectrum

| conversion | reads | percent |
| --- | ---: | ---: |
| C>T | 12877 | 0.793% |
| T>C | 684 | 0.074% |
| C>A | 634 | 0.039% |
| G>A | 517 | 0.029% |
| T>G | 428 | 0.061% |
| A>G | 388 | 0.030% |
| A>C | 354 | 0.028% |
| C>G | 334 | 0.033% |

## Interpretation Notes

- Expected edit model: cytosine base editing (C>T).
- Reads were QC-filtered, unique-read aggregated, aligned to the amplicon in both orientations, and summarized at 1-based amplicon coordinates.
- Edit percentages use position-level covered reads as the denominator; terminal no-coverage gaps are not counted as deletion edits.
- This benchmark is intended as a reproducible validation fixture, not as a claim of equivalence to every CRISPResso2 output table.
