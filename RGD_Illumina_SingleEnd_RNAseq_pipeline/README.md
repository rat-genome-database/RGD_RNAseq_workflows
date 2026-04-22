# RGD Illumina Single-End RNAseq Pipeline

Workflow for processing GEO single-end Illumina RNA-seq data using STAR + RSEM on a SLURM HPC cluster.

## Overview

This pipeline downloads SRA data, performs QC, aligns reads with STAR, and quantifies expression with RSEM. It is designed for bulk RNA-seq datasets from NCBI GEO.

### Pipeline Steps

| Step | Script | Description |
|------|--------|-------------|
| 1 | `run_SRA2QC_SE_v1.bash` | Download SRA reads and run FastQC / FastQ Screen |
| 2 | `run_RNApipeline_SE_diskGuard_v1.bash` | STAR alignment + RSEM quantification |
| — | `bulk_orchestrator_production_diskGuard.bash` | Orchestrates both steps across multiple projects |

### Supporting Scripts

| Script | Description |
|--------|-------------|
| `SRA2QC_SE_v1.sh` | Core SRA download + QC logic |
| `STAR_SE_v1.sh` | STAR alignment (single-end) |
| `RSEM_SE_v1.sh` | RSEM quantification (single-end) |
| `RSEMref_v4.sh` | Build RSEM reference |
| `RSEMmatrix_v6.sh` | Generate count/TPM expression matrices |
| `starRef_v4.sh` | Build STAR genome index |
| `pSTARQC_v1.sh` | Parse STAR alignment QC logs |
| `ComputeSex_v6.sh` | Compute biological sex from expression data |
| `ConflictedSampleReport_v5.sh` | Flag samples with conflicting sex metadata |
| `sex_json_regen_v3.sh` | Regenerate sex inference JSON files |
| `BWjson_v7.sh` | Generate BigWig track JSON for JBrowse |
| `JBrowseSession_v1.sh` | Create JBrowse2 session files |
| `make_jbrowse_session_for_bioproject.py` | Python helper for JBrowse session generation |
| `lib_v10.sh` | Shared library functions |
| `sample_counting.sh` | Count samples per project |

## Requirements

- SLURM workload manager
- [SRA Toolkit](https://github.com/ncbi/sra-tools)
- [STAR](https://github.com/alexdobin/STAR)
- [RSEM](https://github.com/deweylab/RSEM)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/)
- [deepTools](https://deeptools.readthedocs.io/) (for BigWig generation)
- Python 3

## Configuration

Before running, update the `USER CONFIGURATION` block at the top of each script:

```bash
SCRIPT_DIR="/path/to/RGD_Illumina_SingleEnd_RNAseq_pipeline"   # Location of pipeline scripts
myDir="/path/to/home"                                            # Your home/base directory
SCRATCH_BASE="/path/to/scratch"                                  # Scratch filesystem for temp files
REF_GTF="/path/to/genome/GCF_036323735.1_GRCr8/..."            # GTF annotation file
GENOME_FASTA="/path/to/genome/GCF_036323735.1_GRCr8/..."       # Genome FASTA file
```

Also update SLURM directives at the top of each script:
```bash
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --mail-user=your@email.edu
```

## Input Files

Place the following in your working directory before running:

- **AccList file**: One SRA accession per line (see `docs/example_AccList.txt`)
- **Project list file**: Tab-separated `BIOProjectID <tab> AccList_filename` (see `docs/example_project_list.txt`)

## Usage

**Single project (two-step):**
```bash
# Step 1: Download and QC
sbatch run_SRA2QC_SE_v1.bash <BIOProjectID> <AccList_file>

# Step 2: Align and quantify
sbatch run_RNApipeline_SE_diskGuard_v1.bash <BIOProjectID> <AccList_file>
```

**Multiple projects (orchestrated):**
```bash
bash bulk_orchestrator_production_diskGuard.bash project_list.txt
```

## Reference Genome

Scripts are configured for **GRCr8** (rat genome GCF_036323735.1). Commented-out lines for **mRatBN7.2** are retained for reference.

## License

GPL-3.0 — see `dependencies/LICENSE`
