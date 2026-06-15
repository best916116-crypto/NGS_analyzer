# CRISPResso2 EMX1 Base Editor Benchmark

This fixture uses the public CRISPResso2 base-editor example FASTQ file and analyzes it with the local Amplicon Analyzer benchmark runner.

Source:

- CRISPResso2 examples: <https://docs.crispresso.com/latest/core/examples.html>
- FASTQ: <https://crispresso.com/static/demo/base_editor.fastq.gz>

Run:

```bash
node tools/amplicon-analyzer-benchmark.mjs benchmarks/public/crispresso2_base_editor/config.json
```

Outputs are written to `benchmarks/public/crispresso2_base_editor/results/`.

Current reproduced result with an 80% minimum alignment identity threshold:

| Metric | Value |
| --- | ---: |
| QC-passed reads | 25,000 |
| Alignment-passed reads | 24,976 |
| Low-identity reads excluded | 24 |
| Alignment-passed unique reads | 3,915 |
| Mean alignment identity | 99.678% |
| Top target-window edit | C136, 27.709% substitution signal |
| Dominant substitution class | C>T, 12,877 substitution reads |
