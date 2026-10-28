# Validation of 0.1.0 alpha

Date: 2026-10-28. These are executed tests, not a claim of universal biological validity.

## Executed checks

- 21 unit tests passed on Windows/Python 3.13.5. Tests cover plain/gzip paired-read ingestion, unsafe or repeated sample IDs, resource bounds, empty candidates, isoforms, explicit differential gene selection, reference/counting-model/source-TPM consistency, wrong same-length sequences, unpaired-model individual IDs, signature changes and preservation of completed status on rejected resume.
- Source built and installed as a Python wheel locally. CLI version/help and configuration generation were exercised.
- Six synthetic paired-end libraries were processed through real fastp, HISAT2, SAMtools, StringTie, gffcompare, gffread, featureCounts and R analysis on Linux/Python 3.10.22 with PyYAML 6.0.2. Two supported nonreference loci were exported; an unobserved guide model was excluded.
- A resource-only resume changed neither completed alignment nor count-file modification times.
- Final reference-aware summary matched existing 73-library results: 51,165 annotated loci, 10,850 nonreference loci meeting TPM>=1 in any library, 24,693 unique genes with tissue mean TPM>1, and 5,283 nonreference loci under that tissue criterion. 15,319 transcript isoforms were exported for the 10,850 candidates.
- The original eight-library differential module was replayed using the existing integer count matrix: 13,944 tested genes; 1,794 raw-P/FC candidates, including 48 nonreference candidates; zero candidates under FDR<0.05 and |log2FC|>1.
- The P-value primary exporter reproduced those 48 candidates (6 Up, 42 Down), a 48-record representative FASTA and an 87-record all-transcript FASTA. A synthetic end-to-end postprocessing test also exercised the primary P CSV/FASTA/figure output.
- Local DIAMOND BLASTX for the existing 48 sequences reproduced 42 with matches: 25 high-identity broad-protein-coverage, 12 high-identity partial-protein matches, 5 limited/divergent matches, and 6 without a detected match.
- Current-server configuration validation recognized all 73 sample pairs, 328,623 nonreference contigs, and the configured tools, without remapping the real data.
- Read-only review identified import-provenance and rejected-resume-state problems; both were fixed with failing-then-passing regression tests and rechecked.

## Tested tools

The existing server tools were used without replacing the user's environments. Version logs include HISAT2 2.2.1, SAMtools 1.19.2, featureCounts 2.1.1, gffcompare 0.12.10 and gffread 0.12.7. Protein regression used DIAMOND 2.2.2; R analysis used the existing R 4.3.1/edgeR 4.0.16 installation. Full version output and private paths remain in the user's test workspace, not this public archive.

`environment.yml` specifies acceptable dependencies; a fresh conda solver run for that exact file was not executed. It is not a locked environment or a guarantee of identical numeric output with future dependency versions.

## Scope limits

The new software did not remap all 73 real libraries during development; real-result regression used their existing assemblies and quantification with strict source checks. Synthetic reads test execution and logic, not sensitivity, false-discovery rates or biological gene novelty. Multi-day resource scaling, every organism/library layout, single-end/long-read/single-cell data and paired or batch-adjusted differential models were not validated.

Protein names are homolog labels from the configured local library. Online NCBI gene-symbol lookup performed in a separate research analysis is not an automatic feature of this release. No sequence is uploaded by the software. P-value primary outputs are exploratory raw-P selections, not BH-FDR-significant discoveries.

## 0.1.1 homology summary update (2026-10-28)

- 22 unit tests passed, including candidate subset counts, unannotated versus no-match states, primary output and HTML summaries.
- Replayed the existing 48-gene annotation archive: 42 matched, 25 broad/high-identity, 12 partial/high-identity, 5 limited/divergent and 6 without detected matches.
- This update reuses existing annotations; it does not rerun read alignment or protein database searches.

