# Changelog

All notable changes to this pipeline are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [Unreleased]

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
