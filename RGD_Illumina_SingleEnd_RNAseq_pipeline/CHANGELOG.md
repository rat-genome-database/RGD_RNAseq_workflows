# Changelog

All notable changes to this pipeline are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [Unreleased]

---

## [2.1.0] — 2026-10-09

First versioned release of the single-end workflow. Library strandedness and library preparation
(poly(A) selection vs rRNA depletion) are now inferred from the data for every sample, the sex call uses
sex-specific gene expression, and JBrowse sessions include only samples that pass STAR alignment QC.

### Added
- `STAR_SE_v2.sh` — replaces `STAR_SE_v1.sh`. STAR now runs with `--quantMode TranscriptomeSAM GeneCounts`.
  After alignment:
  - **Strand detection** from STAR gene counts: stranded if one strand has more than 3× the reads of the
    other, stranded with a weak-strand warning above 1.5× up to 3×, otherwise unstranded; fewer than 10,000
    stranded gene-assigned reads → unstranded (`AMBIGUOUS_LOW_COUNTS`).
  - **Picard CollectRnaSeqMetrics** with the detected strand (requires the `picard` module and the
    `RRNA_INTERVALS` and `REF_FLAT` reference files; see the README).
  - **Library prep call** through `libprep_call.sh`. Wall time 8 h.
- `libprep_call.sh` — per-sample library prep call, rule **LP2**: NonPolyA_perM (reads on small-RNA genes
  lacking stable poly(A) tails — snoRNA, snRNA, scaRNA, SRP RNA, RNase MRP RNA, telomerase RNA — per million
  gene-assigned reads) ≤ 350 → POLYA, ≥ 700 → RRNA_DEPLETED, between → AMBIGUOUS, < 500,000 gene-assigned
  reads → UNKNOWN. Writes `<GSM>_library_prep.log` with the call, a plain-language reason, rule version,
  thresholds, gene-list MD5, Picard metrics as supporting evidence, and a QC flag (LOW_MAPPING: STAR unique
  < 50% and multimapped > 30%).
- `libprep_recall.sh` — re-calls library prep for projects already processed with `STAR_SE_v2.sh`, without
  realignment.
- `build_nonpolyA_genes.sh` — one-time build of the read-only gene list from the GTF.
- `ConflictedSampleReport_v7.sh` — replaces `ConflictedSampleReport_v5.sh`.
  - **Sex call:** the final call now comes from sex-specific gene TPMs, R = (Xist + 1) / (median of Uty,
    Ddx3y, Kdm5d, Eif2s3y + 1): R ≥ 10 → F, R ≤ 1 → M, otherwise Undetermined / Conflict. Sry is reported but
    not used. The original chrX/chrY calls are kept in `<BIOProjectID>_sex_result.xy.txt`.
  - **New columns:** ComputedStrand, ComputedLibraryPrep and the library prep metrics, rule, QC flag and
    reason, read from each sample's log.
  - **Line 1** of the report gives the workflow version and the sex-call rule.
  - **Study summary:** writes `<BIOProjectID>_library_prep_summary.txt` with a StudyCall (POLYA or
    RRNA_DEPLETED when ≥ 80% of definite sample calls agree, otherwise MIXED).
- `RSEMmatrix_v7.sh` — replaces `RSEMmatrix_v6.sh`; calls `ConflictedSampleReport_v7.sh`.
- `run_RNApipeline_SE_diskGuard_v2.bash` — replaces `run_RNApipeline_SE_diskGuard_v1.bash`; calls
  `STAR_SE_v2.sh`, `RSEMmatrix_v7.sh` and `JBrowseSession_v2.sh`.
- `JBrowseSession_v2.sh` — replaces `JBrowseSession_v1.sh`; builds the session from the STAR QC PASS
  accession list, so QC-failed samples are never included. BigWigs of failed samples stay in the sample folder.
- `CHANGELOG.md` (this file).

### Changed
- `BWjson_v7.sh` — track metadata gains "Computed Strandedness" (Stranded (forward), Stranded (reverse),
  Unstranded) and "Computed Library Capture" (Poly(A) selection, rRNA depletion, Ambiguous); either is
  "Unknown (sample data couldn't be measured)" when no call is available. "Data Processing" now reads
  "HPC RGD single-end workflow 2.1.0" (set by `WORKFLOW_VERSION` at the top of the script).
- `make_jbrowse_session_for_bioproject.py` — includes only samples in the STAR QC PASS accession list
  (optional 4th argument; aborts if the list is missing unless `ALLOW_MISSING_QC=1`).
- `bulk_orchestrator_production_diskGuard.bash` — Step 2 runs `run_RNApipeline_SE_diskGuard_v2.bash`.
- `SRA2QC_SE_v1.sh` — FastQC runs with a 2 GB Java heap (`_JAVA_OPTIONS=-Xmx2g`).
- `README.md` — documents the new steps, reference files and outputs; corrects the usage lines for both
  steps, the AccList description, and the project-list format.

### Validation
The library prep thresholds were calibrated on paired-end data and checked on single-end studies whose
library preparation is documented in GEO:
- rRNA depletion: GSE53960 (rat BodyMap; Ribo-Zero, unstranded, 50 bp; 320 samples, 11 organs) — all 320
  called RRNA_DEPLETED, NonPolyA_perM 1,028–19,255.
- Poly(A) selection: GSE100355 (cerebral cortex, 6 samples) and GSE314094 (left ventricle, 8 samples) — all
  called POLYA, NonPolyA_perM 50–170.
- Blind test: GSE310748 (dorsal root ganglion, 14 samples; GEO description ambiguous) — all called POLYA
  (NonPolyA_perM 20–75), consistent with its "polyA RNA" molecule annotation.

### Removed (superseded)
- `STAR_SE_v1.sh` → `STAR_SE_v2.sh`
- `ConflictedSampleReport_v5.sh` → `ConflictedSampleReport_v7.sh`
- `RSEMmatrix_v6.sh` → `RSEMmatrix_v7.sh`
- `run_RNApipeline_SE_diskGuard_v1.bash` → `run_RNApipeline_SE_diskGuard_v2.bash`
- `JBrowseSession_v1.sh` → `JBrowseSession_v2.sh`

### Removed (not distributed)
- `sex_json_regen_v3.sh` — RGD-specific utility for regenerating sex reports and JBrowse files after manual
  curation.
