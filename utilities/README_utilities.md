# Utilities

Standalone scripts for re-running or correcting specific steps outside of the main RNAseq pipeline. These are used when a project needs targeted fixes after manual QC review.

---

## Scripts

| Script | Purpose |
|--------|---------|
| [`sex_json_regen_v2.sh`](#sex_json_regen_v2sh) | Run after manually correcting sex assignments in the sex result file. Regenerates the sex conflict report and per-sample BigWig JSON (JBrowse2 track) files for a BioProject. |

---

## `sex_json_regen_v2.sh`

**When to use:** Sex assignment is part of manual QC. After reviewing the pipeline output and editing `<BioProject>_sex_result.txt` directly to correct any misassignments, run this script to regenerate the downstream files that depend on it:

1. The sex conflict report (`<BioProject>_sex_conflict_report.txt`)
2. Per-sample BigWig JSON track definition files used by JBrowse2
3. The combined JBrowse2 session JSON

### Prerequisites

- SLURM cluster with `sbatch`, `squeue`, `sacct`
- The shared pipeline library (`lib_v<N>.sh`) and the following sub-scripts must be present in `SCRIPT_DIR`:
  - `ConflictedSampleReport_v<N>.sh`
  - `BWjson_v<N>.sh`
  - `JBrowseSession_v<N>.sh`
- `mail` utility (optional, for failure notifications)

### Configuration

Open the script and fill in the `CONFIGURATION` block near the top:

```bash
SCRIPT_DIR="/path/to/your/pipeline/scripts"   # directory containing lib and sub-scripts
LIB_NAME="lib_v10.sh"                          # shared library filename
DATA_BASE="/path/to/your/expression/data"      # root directory for all BioProject data
SCRATCH_BASE="/scratch/your/scratch/area"      # scratch filesystem for temporary files
REFdir="/path/to/your/reference/genome"        # reference genome directory
NOTIFY_EMAIL="YOUR_EMAIL@example.com"          # failure notification recipient

CONFLICT_SCRIPT="ConflictedSampleReport_v4.sh" # sex conflict report script name
BWJSON_SCRIPT="BWjson_v7.sh"                   # BigWig JSON script name
SESSION_SCRIPT="JBrowseSession_v1.sh"          # JBrowse session builder script name

BWJSON_REQUIRED=true  # set to false to warn-and-continue on BWjson failures
```

Also update the SLURM `#SBATCH` header at the top of the file:
```bash
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --mail-user=YOUR_EMAIL@example.com
```

### Usage

```bash
# Interactive
bash sex_json_regen_v2.sh <SRA_accession_list.txt> <BioProject_ID>

# SLURM submission
sbatch sex_json_regen_v2.sh <SRA_accession_list.txt> <BioProject_ID>
```

**Arguments:**

| Argument | Description |
|----------|-------------|
| `<SRA_accession_list.txt>` | Tab-delimited accession list used by the main pipeline. Must include a header row with columns: `Run`, `geo_accession`, `Tissue`, `Strain`, `Sex`, `PMID`, `GEOpath`, `Title`, `Sample_characteristics`, `StrainInfo` |
| `<BioProject_ID>` | GEO/BioProject accession (e.g. `GSE12345`). Used to locate the project directory and label output files. |

**Example:**
```bash
sbatch sex_json_regen_v2.sh GSE12345_accession_list.txt GSE12345
```

> **Important:** Pass the accession list that reflects the final corrected sample set so that regeneration does not fail due to unprocessed samples.

### What it does

| Step | Action |
|------|--------|
| **1** | Runs `ConflictedSampleReport` to regenerate `<BioProject>_sex_conflict_report.txt` |
| **2** | Submits one `BWjson` SLURM job per unique `geo_accession` in the accession list; waits for all to complete |
| **3** | Submits `JBrowseSession` (dependent on all BWjson jobs) to rebuild the session JSON |

### Output

All output is written relative to `DATA_BASE/<BioProject_ID>/`:

```
<BioProject_ID>/
├── <BioProject_ID>_sex_conflict_report.txt       # regenerated
├── <BioProject_ID>_jbrowse_session.json          # regenerated
├── log_files/
│   ├── <BioProject_ID>_regenSexConflictandjson.out
│   └── JBrowseSession/
│       └── JBrowseSession-<jobid>.out/.err
└── reads_fastq/
    └── <geo_accession>/
        └── log_files/
            └── BWjson/
                └── BWjson-<jobid>.out/.err
```

### BWjson failure behavior

The `BWJSON_REQUIRED` flag controls what happens when a BWjson job fails:

| Value | Behavior |
|-------|----------|
| `true` (default) | Pipeline exits with error; no JBrowse session submitted |
| `false` | Warning email sent; pipeline continues to submit JBrowse session for samples that succeeded |

### Orchestrated mode

If this script is called from a parent orchestration script, set `ORCHESTRATED_MODE=true` in the environment to suppress duplicate email notifications:

```bash
ORCHESTRATED_MODE=true sbatch sex_json_regen_v2.sh GSE12345_accession_list.txt GSE12345
```

---

## Related

- [`combine_workflow/`](../combine_workflow/README.md) — merges two separately processed runs of the same BioProject
- Main pipeline scripts — upstream alignment, quantification, and initial QC steps
