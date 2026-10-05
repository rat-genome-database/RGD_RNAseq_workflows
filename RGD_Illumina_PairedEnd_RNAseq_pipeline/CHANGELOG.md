# Changelog

All notable changes to this pipeline are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [Unreleased]

---

## [2.1.1] — 2026-10-05

### Changed
- `ConflictedSampleReport_v8.sh` — the conflict report columns `Strand` and `LibraryPrep` are renamed
  `ComputedStrand` and `ComputedLibraryPrep`, matching `ComputedSex`. The study summary
  (`<BIOProjectID>_library_prep_summary.txt`) now reads the `ComputedLibraryPrep` column; its format is unchanged.
  Reports written by 2.1.0 keep the old column names until the script is re-run for that project.
- `BWjson_v7.sh` — each BigWig track's metadata (shown in JBrowse) now includes two fields after "Computed Sex",
  read from the sample's `<GSM>_library_prep.log`:
  - **Computed Strandedness:** Stranded (forward), Stranded (reverse) or Unstranded.
  - **Computed RNA Selection:** Poly(A) selection, rRNA depletion, or "Ambiguous (sample data doesn't
    definitively fit either group)".
  - Either field is "Unknown (sample data couldn't be measured)" when the log is missing, as for samples processed
    before 2.1.0. RNA Selection is also Unknown below 500,000 gene-assigned reads.
  - Track colors in the JBrowse session are unchanged.

---

## [2.1.0] — 2026-10-02

Library strandedness and library preparation (poly(A) selection vs rRNA depletion) are now inferred from the
data for every sample, and the sex call uses sex-linked gene expression.

### Added
- `STAR_bigwig3.sh` — replaces `STAR_bigwig2.sh`. STAR now runs with `--quantMode TranscriptomeSAM GeneCounts`.
  After alignment:
  - **Strand detection** from STAR gene counts: stranded if one strand has more than 3× the reads of the
    other, stranded with a weak-strand warning above 1.5× up to 3×, otherwise unstranded; fewer than 10,000
    stranded gene-assigned reads → unstranded (`AMBIGUOUS_LOW_COUNTS`).
  - **Picard CollectRnaSeqMetrics** with the detected strand (requires the `picard` module, `RRNA_INTERVALS`, `REF_FLAT`).
  - **Library prep call** through `libprep_call.sh`. Wall time 8 h.
- `libprep_call.sh` — per-sample library prep call, rule **LP2**: NonPolyA_perM (reads on non-polyadenylated
  small-RNA genes — snoRNA, snRNA, scaRNA, SRP, RNase MRP, telomerase RNA — per million gene-assigned reads)
  ≤ 350 → POLYA, ≥ 700 → RRNA_DEPLETED, between → AMBIGUOUS, < 500,000 gene-assigned reads → UNKNOWN.
  Calibrated on 617 rat samples (24 study sets, 23 GEO series) with a documented prep: 100% correct.
  Writes `<GSM>_library_prep.log` with the call, a plain-language reason, rule version, thresholds, gene-list
  MD5, Picard metrics as supporting evidence, and a QC flag (LOW_MAPPING: STAR unique < 50% and multimapped > 30%).
- `libprep_recall.sh` — re-calls library prep for already-processed projects without realignment.
- `build_nonpolyA_genes.sh` — one-time build of the read-only non-polyA gene list from the GTF.
- `ConflictedSampleReport_v8.sh` — replaces `ConflictedSampleReport_v4.sh`.
  - **Sex call:** the final call now comes from sex-linked gene TPMs, R = (Xist + 1) / (median of Uty, Ddx3y,
    Kdm5d, Eif2s3y + 1): R ≥ 10 → F, R ≤ 1 → M, otherwise Undetermined / Conflict. Sry is reported but not used.
    The original chrX/Y calls are kept in `<BIOProjectID>_sex_result.xy.txt`.
  - **New columns:** Strand, LibraryPrep and the library prep metrics, rule, QC flag and reason, read from each
    sample's log.
  - **Study summary:** writes `<BIOProjectID>_library_prep_summary.txt` with a StudyCall (POLYA or RRNA_DEPLETED
    when ≥ 80% of definite sample calls agree, otherwise MIXED).
- `JBrowseSession_v2.sh` — replaces `JBrowseSession_v1.sh`; builds the session from the STARQC PASS
  accession list, so QC-failed samples are never included. BigWigs of failed samples stay in the download folder.

### Changed
- `run_RNApipeline_pairedG8_diskGuard.bash` — calls `STAR_bigwig3.sh` and `JBrowseSession_v2.sh`; wall time 72 h.
- `RSEMmatrix_v5.sh` — calls `ConflictedSampleReport_v8.sh`; memory 16 GB.
- `RSEM_noBW.bash` — adds `--no-bam-output` (RSEM's own BAM is not used downstream); wall time 24 h.
- `make_jbrowse_session_for_bioproject.py` — includes only samples in the STARQC PASS accession list
  (optional 4th argument; aborts if the list is missing unless `ALLOW_MISSING_QC=1`).
- `BWjson_v7.sh` — loads a Python module before JSON validation.

### Fixed
- **RSEM reference reuse:** the existing-reference check in `run_RNApipeline_pairedG8_diskGuard.bash` and
  `RSEMref_v4.sh` now looks for `<scratch_dir>/rsemref.transcripts.fa`, where `rsem-prepare-reference` writes it.
  The previous checks looked in other locations and never matched, so the reference was rebuilt on every run.
- **Scratch paths in the public scripts:** several scripts used `scratch_dir` without setting it. It is now
  built as `SCRATCH_BASE/<BIOProjectID>` in `run_RNApipeline_pairedG8_diskGuard.bash`, `STAR_bigwig3.sh`,
  `RSEMref_v4.sh`, `starRef_v4.sh`, `ComputeSex_v5.sh` and `pSTARQC_v1.sh`. `run_SRA2QC_diskGuard.bash` sets
  `scratchDir` and `SCRATCH_MOUNT` the same way. `SCRATCH_BASE` now means the same thing in every script.
- `RSEMmatrix_v5.sh` — looks for the custom `rsem-generate-data-matrix` scripts in `../dependencies`
  (the repository's `dependencies/` folder, beside `scripts/`). It previously looked in `scripts/dependencies`,
  so the matrix step failed when the repository was used as laid out.
- `utilities/README_utilities.md` — corrected the link to `combined_project_processing/`.

### Removed (superseded)
- `STAR_bigwig2.sh` → `STAR_bigwig3.sh`
- `ConflictedSampleReport_v4.sh` → `ConflictedSampleReport_v8.sh`
- `JBrowseSession_v1.sh` → `JBrowseSession_v2.sh`

### Removed (not distributed)
- `sex_json_regen_v2.sh` (`scripts/` and `utilities/`) — RGD-specific utility for regenerating sex reports and
  JBrowse files after manual curation.
- `standalone_SRA2QC.sh` — download retry helper outside the pipeline; download SRA data with your own tools.

---

## [2.0.1] — 2026-04-29

### Changed
- `RSEM_noBW.bash` — removed `--sort-bam-by-coordinate`; RSEM's coordinate-sorted BAM was not used by any
  downstream step (sorting is done on the STAR genome BAM).
- `SRA2QC_production.sh` — FastQC runs with `-t 2` and an 8 GB Java heap to prevent out-of-memory failures.

---

## [2.0.0] — 2026-03-17

### Initial public release (GRCr8 workflow)

**New scripts**
- `SRA2QC_production.sh` — replaces `SRA2QC.sh`; adds exponential-backoff prefetch retry (up to 8 attempts), `.sralite` support, `vdb-config` cache redirection, and exit code 2 for single-end layout detection
- `bulk_orchestrator_production_diskGuard.bash` — top-level batch orchestrator with hybrid small/large project scheduling and scratch disk guard

**Updated scripts**
- `STAR_bigwig2.sh` (Feb 2026) — migrated to GRCr8; multi-run support per GSM; BigWig generation via bamCoverage (BPM, unique mappers)
- `ComputeSex_v5.sh` (Feb 2026) — fixed BAM paths to match STAR output structure; processes PASS-only samples
- `RSEMmatrix_v5.sh` — calls custom `rsem-generate-data-matrix` scripts from `dependencies/` via `SCRIPT_DIR`
- `make_jbrowse_session_for_bioproject.py` (Mar 2026) — updated BigWig URL path to `Genome-wide_read_coverage_BigWig_files/`
- `sex_json_regen_v2.sh` — fixed `LIB` and `ConflictedSampleReport` paths to use `SCRIPT_DIR`
- `JBrowseSession_v1.sh` — fixed hardcoded Python script path to use `SCRIPT_DIR`
- `run_SRA2QC_diskGuard.bash` — calls `SRA2QC_production.sh` instead of `SRA2QC.sh`; `screenconfig` set via `SCRIPT_DIR/dependencies/`

**Dependencies added**
- `dependencies/rsem-generate-data-matrix` — custom Perl script (TPM output); adds `File::Basename` for clean column headers (Aug 2024)
- `dependencies/rsem-generate-data-matrix-counts` — custom Perl script (expected counts output)
- `dependencies/fastq_screen.conf` — FastQ-Screen config template with all database paths as placeholders

**Removed scripts** (superseded by diskGuard versions)
- `SRA2QC.sh`
- `run_SRA2QC.bash`
- `run_RNApipeline_pairedG8.bash`
