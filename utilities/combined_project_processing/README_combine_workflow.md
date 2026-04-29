# RNAseq Multi-Strategy / Multi-Batch Combine Workflow

These scripts merge two RNAseq workflow output directories that belong to the same BioProject but were processed in separate runs — for example, paired-end vs. single-end libraries, or libraries with different read lengths — into a single unified output directory.

## Scripts

| Script | Language | Purpose |
|--------|----------|---------|
| `combine_multStrategy_Directories.bash` | Bash | Orchestrates the full merge: result files, matrices, logs, sample data, QC reports, and JBrowse session JSON |
| `make_jbrowse_session_for_combined_bioproject_v2.py` | Python 3 | Builds a JBrowse2 session JSON from per-sample `RNAseq_*.json` track definition files found recursively under the combined project directory |

---

## Prerequisites

- **Bash ≥ 4.0** with GNU coreutils (`join`, `sort`, `cut`, `head`, `tail`, `sed`, `find`, `wc`)
- **Python ≥ 3.8** (standard library only — no additional packages required)
- **SLURM** (optional): the bash script includes `#SBATCH` directives for HPC submission; it can also be run interactively with `bash combine_multStrategy_Directories.bash ...`

---

## Configuration

Before running, open each script and set the variables marked with `# <-- replace`:

### `combine_multStrategy_Directories.bash`

```bash
# Root directory containing all BioProject subdirectories
BASE_PATH="/path/to/your/expression/data"

# Path to the companion JBrowse session builder script
JBROWSE_PY="/path/to/scripts/make_jbrowse_session_for_combined_bioproject_v2.py"
```

Also update the SLURM header if submitting to a cluster:
```bash
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --mail-user=YOUR_EMAIL@example.com
```

If your cluster uses Lmod, the script calls `module load python/3.12.10`. Adjust the module name to match your environment, or remove that line if Python is already on your `PATH`.

### `make_jbrowse_session_for_combined_bioproject_v2.py`

Three variables near the top of `main()` control deployment-specific settings:

```python
# Base URL from which BigWig files are publicly served
BIGWIG_BASE_URL = "https://your.download.server/expression"

# JBrowse2 configuration key for the reference gene/transcript annotation track
GENE_TRACK_CFG = "Your_Species_Assembly_Genes_and_Transcripts-AssemblyName"

# Default viewport (chromosome, coordinates, assembly name)
target_refname = "Chr1"
target_start_1 = 1
target_end_1   = 10_000_000
assembly_name  = "YourAssembly"
```

---

## Usage

### Bash script (interactive)

```bash
bash combine_multStrategy_Directories.bash <dir1> <dir2> <combined_dir>
```

**Arguments:**

| Argument | Description |
|----------|-------------|
| `<dir1>` | Subdirectory name of the first processed run (e.g. `GSE12345_PE`) |
| `<dir2>` | Subdirectory name of the second processed run (e.g. `GSE12345_SE`) |
| `<combined_dir>` | Subdirectory name for the merged output (e.g. `GSE12345`) |

All three names are resolved relative to `BASE_PATH`.

**Examples:**
```bash
# Paired-end vs. single-end
bash combine_multStrategy_Directories.bash GSE12345_PE GSE12345_SE GSE12345

# Different read lengths
bash combine_multStrategy_Directories.bash GSE12345_75bp GSE12345_150bp GSE12345
```

### SLURM submission

```bash
sbatch combine_multStrategy_Directories.bash GSE12345_PE GSE12345_SE GSE12345
```

### Python script (standalone)

The Python script is called automatically by the bash script (Step 7), but can also be run independently:

```bash
python make_jbrowse_session_for_combined_bioproject_v2.py COMBINED_ID PRJdir baseDir
```

**Example:**
```bash
python make_jbrowse_session_for_combined_bioproject_v2.py GSE12345 \
    /path/to/data/expression/GEO/GSE12345 \
    /path/to/data/expression/GEO/GSE12345
```

---

## What Each Step Does

| Step | Action |
|------|--------|
| **1** | Merges `<study>_sex_result.txt` files from both runs; header row written once |
| **2** | Merges `<study>_sex_conflict_report.txt` files; 2-line header written once |
| **3a** | Merges `genes.TPM.matrix` files by joining on gene name (column 1) |
| **3b** | Merges `genes.counts.matrix` files by joining on gene name (column 1) |
| **3c** | Merges `transcripts.TPM.matrix` files by joining on transcript name (column 1) |
| **3d** | Merges `transcripts.counts.matrix` files by joining on transcript name (column 1) |
| **4** | Copies `log_files/` from each run into `combined/log_files/log_files_<dir1>/` and `.../log_files_<dir2>/` |
| **5** | Copies per-sample subdirectories from `reads_fastq/` in each run into the combined `reads_fastq/`; warns on duplicate sample names |
| **6** | Copies the per-run final MultiQC HTML reports into the combined directory |
| **7** | Calls `make_jbrowse_session_for_combined_bioproject_v2.py` to regenerate the JBrowse2 session JSON |

### Matrix merging notes

- Windows-style line endings (`\r\n`) are stripped before merging.
- Trailing blank lines are removed.
- Headers are merged by taking the full header from `dir1` and appending all non-ID columns from `dir2`.
- Rows are joined on the feature name in column 1 (after sorting both files). A warning is issued if the output row count is less than either input, which indicates some feature IDs did not match between the two runs.

---

## Output Directory Structure

After a successful run, the combined directory will contain:

```
GSE12345/
├── GSE12345_sex_result.txt
├── GSE12345_sex_conflict_report.txt
├── GSE12345.genes.TPM.matrix
├── GSE12345.genes.counts.matrix
├── GSE12345.transcripts.TPM.matrix
├── GSE12345.transcripts.counts.matrix
├── GSE12345_jbrowse_session_<AssemblyName>.json
├── GSE12345_PE_final_multiQC_report.html
├── GSE12345_SE_final_multiQC_report.html
├── log_files/
│   ├── GSE12345_combine.log
│   ├── log_files_GSE12345_PE/
│   └── log_files_GSE12345_SE/
└── reads_fastq/
    ├── <sample_A>/
    ├── <sample_B>/
    └── ...
```

---

## JBrowse2 Session JSON

The Python script (`make_jbrowse_session_for_combined_bioproject_v2.py`) finds all `RNAseq_*.json` track definition files recursively under the combined project directory (excluding `*geneTPMbed.json` and `*TXTPMbed.json`) and:

- Updates public-facing metadata URLs (`RGD Metadata Report`, `Project Repository Link`) to use the combined accession ID
- Sets BigWig adapter URIs to the configured `BIGWIG_BASE_URL`
- Assigns a color from a 10-color palette to each unique biological group (Tissue / Strain / Sex / Sample Characteristic combination)
- Writes a complete JBrowse2 session JSON with a LinearGenomeView containing all tracks

Color assignments are printed to stderr for review.

---

## Related Workflows

This combine workflow is designed to run after the main RNAseq processing pipeline. See the other scripts in this repository for upstream alignment, quantification, and QC steps.
