# nonref-rna Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Deliver an installable Linux CLI repository and ZIP implementing the approved nonreference RNA candidate workflow.

**Architecture:** Package the existing tested RNA engine with portable configuration and sample ingestion. Add independent evidence/sequence export and optional homology, R analysis and HTML reports. Validate core algorithms locally and external tools on the existing server.

**Tech Stack:** Python >=3.10, PyYAML, standard library, R, existing external bioinformatics tools.

**Spec:** ../../design.md (approved by user on 2026-10-28).

## Global Constraints

- Linux CLI, paired-end FASTQ, raw or cleaned; reuse prepared reference bundle.
- TPM expression and raw integer counts remain distinct. No unsupported claim of novel functional genes.
- No sequence uploads by default. No data, credentials or personal hostnames in public artifacts.
- Original datasets and outputs are read-only during validation.
- All new code lives in the new project directory; no changes to unrelated existing source.

## Review Focus

- Empty candidate set and no optional database must produce valid explicit outputs.
- FASTQ names, duplicate samples and relative paths must be checked before work.
- Contig selection must not include guide-only predictions or original reference genes.
- Resume after modified inputs/reference must fail; resource-only changes remain allowed.
- Longest isoform export must preserve every selected gene and explicitly identify isoforms.

## Task 1: Portable configuration and ingestion

Files: pyproject.toml, src/nonref_rna/config.py, samples.py, io.py, tests/test_config_samples.py.
Interfaces: load_config(path), load_samples(path), resolve_engine_config(config, work), build_resources(threads,jobs).

- [ ] Write failing tests for resource bounds, duplicate samples, plain/gzip paired FASTQ and relative paths.
- [ ] Implement configuration defaults, PATH tool resolution and sample schema.
- [ ] Run tests to green.

## Task 2: Workflow engine and CLI

Files: src/nonref_rna/engine/*, workflow.py, cli.py, tests/test_workflow.py.
Interfaces: run(config_path,samples_path,work,dry_run=False), validate(...), import_results(...).

- [ ] Add failing tests for changed inputs, nonreference contig identity and raw/clean read handling.
- [ ] Adapt existing tested engine with package imports, explicit contig selection, resource allocation and robust checkpoints.
- [ ] Add init, validate, run --dry-run, summarize and status commands; include per-run lock and provenance.
- [ ] Run unit and CLI tests to green.

## Task 3: Candidate evidence and FASTA/homology exports

Files: evidence.py, homology.py, report.py, tests/test_evidence.py.
Interfaces: summarize(work, output, config), annotate(query, db, ...).

- [ ] Add failing fixture tests for any-sample versus tissue mean, empty results, multiple isoforms and absent database.
- [ ] Implement observed RNA-only selection, deterministic per-gene longest sequence, model/FASTA consistency checks and counts/TPM reporting.
- [ ] Add BLASTX evidence tiers and gene-family ambiguity without definitive novelty claims.
- [ ] Reproduce current buffalo expression statistics and 48 selected FASTA mapping.

## Task 4: Optional R analyses, documentation and release

Files: scripts/analysis.R, examples/*, environment.yml, README.md, README.zh-CN.md, docs/methods.md, .github/workflows/tests.yml, tests/test_integration.py.

- [ ] Write test for group comparison validation and empty DEG outputs.
- [ ] Add raw counts edgeR/TMM with separate FDR/raw-P tables, TPM PCA/MDS and transparent scientific exports.
- [ ] Run a real external-tool synthetic smoke test; record exactly which stages were exercised.
- [ ] Run package install, CLI, tests and static checks, review for privacy, package GitHub-ready ZIP excluding all run data.

## Execution ledger

- Baseline: 13 tests in existing RNA code pass locally.
- Ruling: implement inline in a dedicated new directory rather than create a worktree: no existing tracked source is edited.
- Ruling: include a code-use notice and a license selection note, not an invented author identity or license grant.
- Task 1: complete. Configuration/ingestion tests passed; relative paths, duplicate IDs, paired records and CPU budgets covered.
- Task 2: complete. Real six-library raw-read pipeline completed. Empty splice-site indexing fixed after reproducing the HISAT2 error; resource-only resume preserved completed alignment/count artifacts.
- Task 3: complete. Source counting GTF, transcript models, quantification and reference-derived sequences are cross-checked. Real73 and real48 homology regression passed.
- Task 4: complete. R differential module reproduced 1,794/48 raw-P candidates; P-primary CSV and 48/87 sequence exports passed. Source/wheel packaging and privacy checks executed; see VALIDATION.md for exact scope.
- User steering: primary differential output uses raw P<0.05 and |log2FC|>1, with P-sorted CSV/FASTA/transparent figures under primary_Pvalue. FDR remains a separate reference result.
- Review: two Important findings resolved and focused re-review confirmed closure. Added regression for a reversed exported model despite a plus-strand counting annotation.
- Ruling: deliver a source ZIP as requested; do not create/push a repository, merge unrelated workspace history, or choose the author's license.
