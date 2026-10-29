# Changelog

## 0.1.2 — 2026-10-29

- Standardize report and plot wording to P; use P column headers in primary CSV files.
- Remove repeated explanatory subtitles and replace ambiguous wording with explicit evidence definitions.
- Add an unmarked P volcano plot alongside the existing nonreference-marked version.
- Preserve count analysis, P/FDR calculations and selection thresholds.
- Rename candidate MDS figure files from rawP to P.


## 0.1.1 — 2026-10-28

- Add gene-level protein homology counts and English/Chinese summaries to Markdown, JSON and HTML outputs.
- Report expression and P-value differential candidate sets separately, distinguishing unannotated genes from searched genes without matches.
- Replace ambiguous completeness wording with explicit identity, target coverage, alignment length and E-value criteria.

## 0.1.0 — 2026-10-28

- Package the paired-end RNA workflow as a Linux CLI with portable configuration and sample sheets.
- Add observed nonreference assembly, fixed-model TPM/counts, evidence and all/longest-isoform sequence exports.
- Add local protein-homology evidence, optional R analyses and transparent figure export.
- Make P-value/FC differential candidates the readable primary output, retaining complete FDR statistics separately.
- Add strict import provenance, reference-derived sequence checks, CPU/disk controls and safe resume behavior.
- Research alpha; see validation scope and license-selection notes before public release.
