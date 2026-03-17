# RGD Bulk RNA-seq Pipeline (GRCr8)

A SLURM-based pipeline for processing paired-end bulk RNA-seq data from NCBI SRA. Designed for HPC clusters with SLURM job scheduling. Produces STAR-aligned BAM files, BigWig coverage tracks, RSEM expression matrices, sex estimation, and JBrowse2 visualization sessions.

**Genome:** GCF_036323735.1 GRCr8 (rat)  
**Developed at:** Rat Genome Database (RGD), Medical College of Wisconsin

---

## Pipeline Overview

```
Input: Accession list (TSV) + BioProject ID
       │
       ▼
STEP 1: SRA Download & QC
  ├── prefetch (up to 8 attempts with exponential backoff)
  ├── vdb-validate
  ├── fasterq-dump → paired FASTQ (3 retries; exits code 2 if single-end detected)
  ├── FastQC + FastQ-Screen (per sample)
  └── MultiQC report
       │
       ▼
STEP 2: Alignment, Quantification & Visualization
  ├── STAR genome index (GRCr8)
  ├── STAR alignment → sorted BAM + BigWig (BPM normalized)
  ├── STARQC alignment rate summary → PASS/FAIL filter (threshold: <50% unmapped)
  ├── Sex estimation (chrX/Y read depth ratio)
  ├── RSEM reference + per-sample quantification (PASS samples only)
  ├── RSEM expression matrices (TPM + counts, genes + transcripts)
  ├── MultiQC final report
  ├── Sex conflict report (Xist, Uty, Sry, Ddx3y, Kdm5d, Eif2s3y TPMs)
  ├── BigWig JSON track configs (for JBrowse2)
  └── JBrowse2 session JSON
```

---

## Requirements

### Software (via `module load` on your cluster)

| Tool | Version used |
|------|-------------|
| STAR | 2.7.10b |
| RSEM | 1.3.3 |
| samtools | 1.20 |
| deeptools | 3.5.1 |
| SRA Toolkit | 3.1.1 |
| FastQC | 0.11.9 |
| FastQ-Screen | 0.15.2 |
| Bowtie2 | 2.5 |
| MultiQC | 1.18 |
| Python | 3.9+ |
| Perl | 5.x |
| pigz | any |

### Reference Files (GRCr8)
Download from NCBI accession `GCF_036323735.1`:
- Genome FASTA: `GCF_036323735.1_GRCr8_genomic.fna`
- Genome GTF: `GCF_036323735.1_GRCr8_genomic.gtf`

### FastQ-Screen Genomes
Download the standard genome set using FastQ-Screen's built-in tool:
```bash
fastq_screen --get_genomes --outdir /path/to/your/FastQ_Screen_Genomes
```
Then update the paths in `dependencies/fastq_screen.conf`.

---

## Repository Structure

```
RGD_RNAseq_pipeline/
├── scripts/                  Pipeline scripts (see table below)
├── dependencies/
│   ├── rsem-generate-data-matrix        Custom Perl script (outputs TPM matrix)
│   ├── rsem-generate-data-matrix-counts Custom Perl script (outputs counts matrix)
│   └── fastq_screen.conf               FastQ-Screen config template
├── docs/
│   ├── example_AccList.txt             Example accession list format
│   └── example_project_list.txt        Example orchestrator project list
├── README.md
├── CHANGELOG.md
├── CONTRIBUTING.md
└── LICENSE
```

---

## Configuration

Before running, update the `USER CONFIGURATION` block at the top of each script. The variables to set are:

| Variable | Description |
|----------|-------------|
| `SCRIPT_DIR` | Full path to the `scripts/` directory |
| `myDir` | Your home/base directory |
| `SCRATCH_BASE` | Scratch filesystem mount point |
| `REF_GTF` | Path to GRCr8 GTF file |
| `GENOME_FASTA` | Path to GRCr8 genome FASTA |
| `--account` | Your SLURM account name (in `#SBATCH` headers) |
| `--mail-user` | Your email for SLURM failure notifications |

Also update `dependencies/fastq_screen.conf` with paths to your local FastQ-Screen genome indices.

---

## Input Format

Tab-delimited accession list with a header row:

```
Run	geo_accession	Tissue	Strain	Sex	PMID	GEOpath	Title	Sample_characteristics	StrainInfo
SRR12345678	GSM1234567	Liver	BN/NHsdMcwi	M	12345678	https://...	Study Title	age: 12 weeks	https://...
```

Multiple SRA runs per GSM sample are supported — they are automatically grouped and merged for alignment.

---

## Usage

### Recommended: Batch Mode (Orchestrator)

```bash
# project_list.txt: <acclist_path>  <BioProjectID>  <read_length>
sbatch scripts/bulk_orchestrator_production_diskGuard.bash project_list.txt
```

The orchestrator enforces concurrent project limits, disk space guards, and handles large vs. small project scheduling automatically.

### Single Project

```bash
# Step 1: Download and QC
sbatch scripts/run_SRA2QC_diskGuard.bash /path/to/AccList.txt PRJNA123456

# Step 2: Alignment, quantification, visualization
sbatch scripts/run_RNApipeline_pairedG8_diskGuard.bash /path/to/AccList.txt PRJNA123456 150
```

### Utility Scripts

```bash
# Retry a failed SRA download
sbatch scripts/standalone_SRA2QC.sh PRJNA123456 SRR12345678 GSM1234567

# Regenerate sex conflict report and BigWig JSONs after correcting sex metadata
sbatch scripts/sex_json_regen_v2.sh /path/to/PASS_AccList.txt PRJNA123456

# Count unique samples in an accession list
bash scripts/sample_counting.sh /path/to/AccList.txt
```

---

## Script Reference

| Script | Role |
|--------|------|
| `bulk_orchestrator_production_diskGuard.bash` | Batch orchestrator with hybrid scheduling and disk guard |
| `run_SRA2QC_diskGuard.bash` | Step 1 controller: SRA download and QC |
| `run_RNApipeline_pairedG8_diskGuard.bash` | Step 2 controller: alignment through visualization |
| `SRA2QC_production.sh` | Per-sample SLURM job: SRA → FASTQ → FastQC/FastQ-Screen |
| `standalone_SRA2QC.sh` | Standalone retry for failed SRA downloads |
| `starRef_v4.sh` | Generates STAR genome index |
| `STAR_bigwig2.sh` | STAR alignment → sorted BAM → BigWig (BPM) |
| `pSTARQC_v1.sh` | Parses STAR logs; generates PASS/FAIL alignment summary |
| `ComputeSex_v5.sh` | Estimates sex from chrX/Y read depth ratio |
| `ConflictedSampleReport_v4.sh` | Cross-checks input vs. computed sex using sex-linked gene TPMs |
| `RSEMref_v4.sh` | Generates RSEM reference |
| `RSEM_noBW.bash` | Per-sample RSEM expression quantification |
| `RSEMmatrix_v5.sh` | Combines RSEM results into project-level matrices + final MultiQC |
| `BWjson_v7.sh` | Generates per-sample JBrowse2 BigWig track JSON |
| `JBrowseSession_v1.sh` | SLURM wrapper for the JBrowse session builder |
| `make_jbrowse_session_for_bioproject.py` | Builds combined JBrowse2 session JSON with color grouping |
| `sex_json_regen_v2.sh` | Utility: regenerate sex report + JSONs after metadata correction |
| `sample_counting.sh` | Utility: count unique samples in an accession list |
| `lib_v10.sh` | Shared helper functions (logging, job submission, status checking) |

---

## Output Files

All permanent outputs go to `<myDir>/data/expression/GEO/<BIOProjectID>/`:

```
<BIOProjectID>/
├── reads_fastq/<GSM_ID>/
│   ├── RNAseq_<unique_name>.bigwig
│   ├── RNAseq_<unique_name>.json
│   └── log_files/
├── log_files/
│   ├── STARQC/<BIOProjectID>_STAR_Align_sum.txt
│   └── ...
├── <BIOProjectID>.genes.TPM.matrix
├── <BIOProjectID>.genes.counts.matrix
├── <BIOProjectID>.transcripts.TPM.matrix
├── <BIOProjectID>.transcripts.counts.matrix
├── <BIOProjectID>_sex_result.txt
├── <BIOProjectID>_sex_conflict_report.txt
├── <BIOProjectID>_fastq_multiQC_report.html
├── <BIOProjectID>_final_multiQC_report.html
└── <BIOProjectID>_jbrowse_session_GRCr8.json
```

---

## Custom Dependencies

The `dependencies/` directory contains two custom Perl scripts derived from RSEM that output **TPM values** and **expected counts** respectively, with filenames extracted using `File::Basename` (added August 2024) for cleaner matrix column headers.

---

## Citation

If you use this pipeline, please cite:
- Dobin et al. (2013) STAR. *Bioinformatics* 29(1):15–21
- Li & Dewey (2011) RSEM. *BMC Bioinformatics* 12:323
- Andrews S. (2010) FastQC. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- Ewels et al. (2016) MultiQC. *Bioinformatics* 32(19):3047–3048
- Ramírez et al. (2016) deepTools2. *Nucleic Acids Research* 44(W1):W160–W165

---

## License

Pipeline scripts: MIT — see [LICENSE](LICENSE)  
RSEM-derived scripts in `dependencies/`: GPL-3.0 — see [dependencies/LICENSE](dependencies/LICENSE)
