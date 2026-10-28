# nonref-rna

**RNA-supported nonreference gene candidates from a prepared extended genome.**

[Scientific interpretation](docs/methods.md) · [Design](docs/design.md)

`nonref-rna` packages a buffalo RNA-seq workflow into a Linux command-line tool. Supply paired-end reads, a samples TSV and the reference bundle you already constructed. It aligns reads, assembles observed transcripts, quantifies a fixed annotation, exports nonreference candidate FASTA/GTF files and reports expression, prediction and optional protein-homology evidence separately.

 It reports RNA expression, transcript structure and protein homology separately; these analyses do not establish locus-specific protein function. See the evidence columns and validation report.

## Install

```bash
conda env create -f environment.yml
conda activate nonref-rna
python -m pip install .
nonref-rna --help
```

The environment file is a dependency specification, not a solved lockfile. Record actual tool versions from `logs/versions.txt` and the homology manifest. Linux and Python >=3.10 are required for full runs. Summaries and unit tests can also run on Windows.

## Start with an existing reference bundle

```bash
nonref-rna init --reference-bundle /path/to/buffalo_annotation --out config.local.yaml
# Edit config.local.yaml and create samples.tsv (see examples/).
nonref-rna validate --config config.local.yaml --samples samples.tsv
nonref-rna run --config config.local.yaml --samples samples.tsv --out results/run1 --dry-run
nonref-rna run --config config.local.yaml --samples samples.tsv --out results/run1
```

Reference layout:

```text
buffalo_annotation/
  rnaseq_reference/extended.fa
  rnaseq_reference/guide.gtf
  rnaseq_reference/nonref.guide.gtf
  reference/reference.normalized.gtf
  run_config.json                         # optional legacy source metadata
  blastp/gene_ids.80_80.txt                 # optional prior protein support
```

The GTF files must use IDs present in the extended FASTA. Original reference and nonreference contig sets must be disjoint. By default nonreference IDs begin with `NR_`; use an explicit contig list for other naming schemes.

## Samples

Tab-separated, with no duplicate sample IDs or reused FASTQ files:

```text
sample  r1                       r2                       tissue   breed    group  individual  strandedness
Early1  /reads/Early1_R1.fastq.gz /reads/Early1_R2.fastq.gz Mammary  Buffalo  Early  animal1     auto
Late1   /reads/Late1_R1.fastq.gz  /reads/Late1_R2.fastq.gz  Mammary  Buffalo  Late   animal2     auto
```

Only `sample`, `r1`, `r2` are required. The example above is visual; use real TAB characters, as in `examples/samples.tsv`. Paths are relative to the samples file. `.fastq` and `.fastq.gz` are accepted. `strandedness`: `auto`, `rf`, `fr`, `unstranded`. `individual` records metadata; the current differential model is an **unpaired two-group** model. Paired or repeated-measure designs must be analyzed with a suitable model outside this release.

Set `reads.cleaned: true` to reuse cleaned reads, or `false` to run fastp. Sample exclusions must be performed explicitly in the sample sheet; no breed or tissue is automatically excluded.

## Evidence and outputs

Default candidates: observed nonreference RNA loci with **TPM >=1 in at least one library**. Change `expression.tpm` and `expression.min_samples` to require more support. Tissue summaries separately use **arithmetic tissue mean TPM >1**, with unique gene unions. A sum of per-tissue sizes is not the number of unique genes.

```text
results/run1/
  run_manifest.json                        # input signatures, code, parameters
  status.json
  logs/                                   # commands and actual tool versions
  annotation/                             # observed merge and final counting GTF
  results/                                # engine matrices and evidence
  matrices/gene_TPM.tsv
  matrices/gene_counts.raw.tsv
  candidates/nonreference_candidates.tsv
  candidates/all_transcripts.fa
  candidates/longest_transcript_per_gene.fa
  candidates/gene_transcript_mapping.tsv
  candidates/selected_models.gtf
  homology/gene_homology.tsv               # if enabled
  differential_expression/                # if enabled
  figures/                                # if enabled
  report/summary.json
  report/report.html
```

FASTA files contain spliced nucleotide transcripts. The representative FASTA has one longest isoform per selected gene, with deterministic transcript-ID tie breaking; the complete isoform set is retained. AUGUSTUS comparison is structural support relative to the supplied predicted models.

### Protein annotation

Set `homology.enabled: true` and `homology.database: /path/to/combined_db.dmnd`. The workflow runs local DIAMOND BLASTX and retains raw matches, protein descriptions, identity, query and subject coverage, alignment lengths and E-values. It distinguishes broad protein matches from partial/weak matches and flags repeated matches to the same accession. Gene symbols are extracted only when present in database headers; otherwise use the protein description/accession. No online gene-symbol resolution or sequence upload is performed. 

### Differential expression and plots

`homology/homology_summary.md` and `.json` count matches for the expression candidate set. `primary_Pvalue/homology_summary.md` and `.json` independently count the raw-P differential candidate set. Both include reusable English/Chinese text, also displayed in the HTML report. Counts are per gene, not per transcript or protein hit. A broad high-identity match requires identity >=80%, target protein coverage >=80%, alignment >=50 aa and E-value <=1e-10. Protein coverage does not measure transcript completeness. Unannotated candidates are reported separately from searched candidates without matches; homology counts do not measure AUGUSTUS/assembly dual support.

Set `analysis.differential: true`, `control`, `treatment`, and use at least two independent biological samples per group. Counts are integer featureCounts fragment counts; edgeR uses CPM filtering, TMM, robust dispersion estimation and a QL F-test. Raw-P (`abs(log2FC)>1, P<0.05`) and BH-FDR (`abs(log2FC)>1, FDR<0.05`) results are **separate files**. Raw-P results are exploratory. The software cannot establish that samples are independent from their filenames.

Set `analysis.plots: true` for transparent R PNG/PDF expression, tissue-intersection, PCA and (with differential analysis) volcano/MDS outputs. PCA/MDS use log2(TPM+1), no gene variance scaling. 

**Primary differential results use raw P<0.05 and |log2FC|>1.** Open `primary_Pvalue/nonreference_candidates_P05.csv` for the P-sorted nonreference candidates, their available protein descriptions and representative transcript IDs. The same directory contains `nonreference_candidates_P05.longest.fa`, all-isoform FASTA, GTF, and (when plots are enabled) `Volcano_Pvalue.png` and `MDS_TPM.png`. `all_genes_P05.csv` includes reference and nonreference raw-P candidates. Primary CSVs emphasize P values; FDR results remain in `differential_expression/` for reference. 

## Resume, resources, existing results

Repeat the same `run` command to resume. Inputs/references/tools/biological parameters are signed with paths, sizes and modification times; code is SHA256-hashed. FASTQ contents are not fully hashed by default. Changed biological inputs require a **new output directory**. Resource-only changes are allowed. Index and sample artifacts are checked before reuse; modified/missing artifacts are rebuilt. A Linux advisory lock prevents concurrent writers to one run.

`resources.threads` is the planned tool CPU budget; worker allocation accounts for concurrent aligner and sorter threads. External libraries and OS I/O can add small overhead. Compressed alignments, indexes and scratch files can require hundreds of GB for real datasets. Change `minimum_free_gib` deliberately; disk reservations stop new workers before exhausting the configured headroom.

Export a compatible historical workspace without mapping again:

```bash
nonref-rna summarize --work /path/to/completed/bam-new --out results/imported
# Optional: --config config.local.yaml for thresholds / homology / R analyses
```

This requires its `samples.tsv`, `results/gene_TPM.tsv`, `results/gene_metadata.tsv`, `annotation/final_annotation.gtf`, `annotation/nonref.rna_merged.gtf`, `annotation/nonref.rna_transcripts.fa`, source `quant/<sample>/gene_abundance.tsv`, source engine config and input manifest, and the unchanged extended FASTA referenced by that manifest. Counts are additionally required for differential analysis. IDs, sample columns, counting annotation, source TPM values and transcript sequences reconstructed from reference/exons are checked. Do not mix `BFRNA` IDs from independent assemblies.

Validation checks the first FASTQ pair before alignment, and additionally checks the pilot pairs when automatic strand inference runs. It does not perform a separate full-file read-ID pairing audit; use correctly paired input files and inspect the mapping/QC logs.

## Tests and reproducibility

```bash
python -m unittest discover -s tests -v
# Real external-tool synthetic test, after installation:
python tests/make_smoke_fixture.py /tmp/nonref-smoke
nonref-rna run --config /tmp/nonref-smoke/config.yaml --samples /tmp/nonref-smoke/samples.tsv --out /tmp/nonref-smoke/run
```

Synthetic data are explicitly marked and are not biological discoveries. See `docs/VALIDATION.md` for tests actually executed; CI unit tests alone do not certify real-data biological accuracy. There is no automatic KEGG enrichment or proof of locus uniqueness in this version.

## Upstream documentation

- [HISAT2 manual](https://daehwankimlab.github.io/hisat2/manual/)
- [StringTie manual](https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual)
- [featureCounts](https://subread.sourceforge.net/featureCounts.html)
- [DIAMOND command reference](https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options)
