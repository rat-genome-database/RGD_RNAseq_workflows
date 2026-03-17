# RGD RNAseq Workflows

Bioinformatics workflows developed at the [Rat Genome Database (RGD)](https://rgd.mcw.edu), Medical College of Wisconsin, for processing bulk RNA-seq data from NCBI GEO/SRA on SLURM-based HPC clusters.

All workflows align reads to rat reference genomes, quantify expression with RSEM, and generate outputs compatible with [JBrowse2](https://jbrowse.org) for visualization on the RGD genome browser.

---

## Current Workflow

### [RGD_Illumina_PairedEnd_RNAseq_pipeline](./RGD_Illumina_PairedEnd_RNAseq_pipeline/)

The current production pipeline for **paired-end Illumina** RNA-seq data.

**Reference genome:** GCF_036323735.1 GRCr8  
**Key tools:** STAR 2.7.10b · RSEM 1.3.3 · samtools 1.20 · deeptools 3.5.1 · FastQC · FastQ-Screen · MultiQC

**What it does:**
- Downloads SRA data with robust retry logic (up to 8 prefetch attempts)
- Runs FastQC and FastQ-Screen QC on raw reads
- Aligns to GRCr8 with STAR; generates BPM-normalized BigWig coverage tracks
- Filters samples by alignment rate; estimates biological sex from chrX/Y read depth
- Quantifies expression with RSEM; produces TPM and count matrices
- Generates JBrowse2 track configs and session JSON for visualization
- Supports batch processing of multiple projects via an orchestrator with disk space guards

See the [pipeline README](./RGD_Illumina_PairedEnd_RNAseq_pipeline/README.md) for full documentation, configuration instructions, and usage.

---

## Archive

Older workflow versions are preserved in the [`archive/`](./archive/) directory for reference. These were developed for earlier reference genomes and pipeline architectures and are no longer actively maintained.

| Directory | Description |
|-----------|-------------|
| `GRC8_bulk_illuminaPaired_oneSRA` | GRCr8 paired-end, one SRA run per sample |
| `illuminaPaired_multiplePairs_perSample` | Paired-end, multiple read pairs per sample |
| `illuminaPaired_multipleSRA_IDperSample_v2` | Paired-end, multiple SRA IDs per sample (v2) |
| `illuminaPaired_multipleSRA_IDperSample_v3` | Paired-end, multiple SRA IDs per sample (v3) |
| `illuminaPaired_oneSRA_IDperSample` | Paired-end, one SRA ID per sample |
| `illuminaPaired_oneSRA_IDperSample_v2` | Paired-end, one SRA ID per sample (v2) |
| `illuminaSingle_1readPerSample` | Single-end, one read file per sample |
| `illuminaSingle_1readPerSample_updateAug2025` | Single-end, August 2025 update |
| `illuminaSingle_1readPerSample_v3` | Single-end (v3) |
| `illuminaSingle_multRunsPerSample` | Single-end, multiple runs per sample |

---

## Planned Workflows

- **RGD_Illumina_SingleEnd_RNAseq_pipeline** — in development

---

## License

Pipeline scripts: MIT — see [LICENSE](./RGD_Illumina_PairedEnd_RNAseq_pipeline/LICENSE)  
RSEM-derived scripts in `dependencies/`: GPL-3.0 — see [dependencies/LICENSE](./RGD_Illumina_PairedEnd_RNAseq_pipeline/dependencies/LICENSE)
