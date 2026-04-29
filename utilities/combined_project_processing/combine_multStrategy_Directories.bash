#!/usr/bin/env bash
#SBATCH --job-name=CMBD
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=06:30:00
#SBATCH --account=YOUR_SLURM_ACCOUNT       # <-- replace with your SLURM account
#SBATCH --output=%x-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=YOUR_EMAIL@example.com  # <-- replace with your email
#SBATCH --error=combineDir.err

#####################################################################################################################
# combine_multStrategy_Directories.bash
#
# Combines two RNAseq workflow output directories that belong to the same BioProject
# but were processed in separate runs (e.g. paired-end vs. single-end, or different
# read-length batches) into a single merged output directory.
#
# Usage:
#   bash combine_multStrategy_Directories.bash <dir1> <dir2> <combined_dir>
#
# Arguments:
#   <dir1>          Subdirectory name of the first processed run  (e.g. GSE70012_PE)
#   <dir2>          Subdirectory name of the second processed run (e.g. GSE70012_SE)
#   <combined_dir>  Subdirectory name for the merged output       (e.g. GSE70012)
#
# All three names are resolved relative to BASE_PATH (set below).
#
# Examples:
#   bash combine_multStrategy_Directories.bash GSE70012_PE  GSE70012_SE  GSE70012
#   bash combine_multStrategy_Directories.bash GSE70012_75bp GSE70012_150bp GSE70012
#
# Steps performed:
#   1. Merge <study>_sex_result.txt files           (header preserved once)
#   2. Merge <study>_sex_conflict_report.txt files  (2-line header preserved once)
#   3. Merge gene/transcript TPM and counts matrices on matching feature IDs (col 1)
#      3a. genes.TPM.matrix
#      3b. genes.counts.matrix
#      3c. transcripts.TPM.matrix
#      3d. transcripts.counts.matrix
#   4. Copy log_files/ subdirectories from each run into combined log_files/
#   5. Copy reads_fastq/ sample subdirectories from each run into combined reads_fastq/
#   6. Copy per-run final MultiQC HTML reports into the combined directory
#   7. Regenerate the JBrowse2 session JSON with make_jbrowse_session_for_combined_bioproject_v2.py
#
# Dependencies:
#   - bash >= 4.0, GNU coreutils (join, sort, cut, head, tail, sed, find, wc)
#   - Python 3 with make_jbrowse_session_for_combined_bioproject_v2.py accessible
#   - SLURM environment (optional; can also be run interactively with: bash <script> ...)
#
# Output:
#   All merged files are written to <BASE_PATH>/<combined_dir>/
#   A timestamped combine log is written to <BASE_PATH>/<combined_dir>/log_files/
#####################################################################################################################
set -euo pipefail

# ---------------------------------------------------------------------------
# CONFIGURATION — edit these two variables before running
# ---------------------------------------------------------------------------
# Root directory containing all BioProject subdirectories
BASE_PATH="/path/to/your/expression/data"   # <-- replace with your base data path

# Path to the companion JBrowse session builder script
JBROWSE_PY="/path/to/scripts/make_jbrowse_session_for_combined_bioproject_v2.py"  # <-- replace with actual path
# ---------------------------------------------------------------------------

# --- Argument parsing ---
if [[ $# -ne 3 ]]; then
    echo "Usage: $0 <dir1> <dir2> <combined_dir>"
    echo "Example: $0 GSE70012_PE GSE70012_SE GSE70012"
    echo "Example: $0 GSE70012_75bp GSE70012_150bp GSE70012"
    exit 1
fi

DIR1_NAME="${1}"
DIR2_NAME="${2}"
COMBINED_DIR_NAME="${3}"

# --- Derive study ID from combined dir name (used for output filenames) ---
STUDY_ID="${COMBINED_DIR_NAME}"

# --- Full paths ---
DIR1="${BASE_PATH}/${DIR1_NAME}"
DIR2="${BASE_PATH}/${DIR2_NAME}"
COMBINED_DIR="${BASE_PATH}/${COMBINED_DIR_NAME}"

# --- Validate input directories exist ---
if [[ ! -d "$DIR1" ]]; then
    echo "ERROR: Directory 1 not found: $DIR1" >&2
    exit 1
fi
if [[ ! -d "$DIR2" ]]; then
    echo "ERROR: Directory 2 not found: $DIR2" >&2
    exit 1
fi

# --- Create combined output directory ---
mkdir -p "$COMBINED_DIR"
echo "Base path       : $BASE_PATH"
echo "Directory 1     : $DIR1"
echo "Directory 2     : $DIR2"
echo "Combined output : $COMBINED_DIR"

LOG_FILE="${COMBINED_DIR}/log_files/${STUDY_ID}_combine.log"
mkdir -p "$(dirname "$LOG_FILE")"
exec > >(tee -a "$LOG_FILE") 2>&1

# -------------------------------------------------------
# STEP 1: Merge sex_result.txt files (header only once)
# -------------------------------------------------------
# Fixes applied:
#   - sed 's/\r//g' : strips Windows-style CRLF line endings that would
#                     corrupt the header and data rows
#   - grep -v "^$"  : drops any trailing blank lines from input files
#                     so they don't propagate into the combined output
# -------------------------------------------------------
DIR1_SEX="${DIR1}/${DIR1_NAME}_sex_result.txt"
DIR2_SEX="${DIR2}/${DIR2_NAME}_sex_result.txt"
COMBINED_SEX="${COMBINED_DIR}/${STUDY_ID}_sex_result.txt"

echo ""
echo "=== Merging sex result files ==="
echo "  Dir1 source : $DIR1_SEX"
echo "  Dir2 source : $DIR2_SEX"
echo "  Output      : $COMBINED_SEX"

if [[ ! -f "$DIR1_SEX" ]]; then
    echo "  WARNING: Dir1 sex result file not found: $DIR1_SEX" >&2
fi
if [[ ! -f "$DIR2_SEX" ]]; then
    echo "  WARNING: Dir2 sex result file not found: $DIR2_SEX" >&2
fi

if [[ -f "$DIR1_SEX" && -f "$DIR2_SEX" ]]; then
    {
        head -n 1 "$DIR1_SEX"    | sed 's/\r//g'
        tail -n +2 "$DIR1_SEX"   | sed 's/\r//g' | grep -v "^$"
        tail -n +2 "$DIR2_SEX"   | sed 's/\r//g' | grep -v "^$"
    } > "$COMBINED_SEX"
    echo "  Dir1 rows added : $(tail -n +1 "$DIR1_SEX" | grep -v "^$" | wc -l | tr -d ' ')"
    echo "  Dir2 rows added : $(tail -n +2 "$DIR2_SEX" | grep -v "^$" | wc -l | tr -d ' ')"
    echo "  Total data rows : $(tail -n +2 "$COMBINED_SEX" | wc -l | tr -d ' ')"
elif [[ -f "$DIR1_SEX" ]]; then
    sed 's/\r//g' "$DIR1_SEX" | grep -v "^$" > "$COMBINED_SEX"
    echo "  WARNING: Only Dir1 file found - copied as-is."
elif [[ -f "$DIR2_SEX" ]]; then
    sed 's/\r//g' "$DIR2_SEX" | grep -v "^$" > "$COMBINED_SEX"
    echo "  WARNING: Only Dir2 file found - copied as-is."
else
    echo "  ERROR: Neither sex result file found. Skipping." >&2
fi

echo ""
echo "Step 1 complete. Ready for next step (sex_conflict_report)."

# -------------------------------------------------------
# STEP 2: Merge sex_conflict_report.txt files (2-line header only once)
# -------------------------------------------------------
# Fixes applied:
#   - sed 's/\r//g' : strips Windows-style CRLF line endings that would
#                     corrupt the header and data rows
#   - grep -v "^$"  : drops any trailing blank lines from input files
#                     so they don't propagate into the combined output
# -------------------------------------------------------
DIR1_SEX="${DIR1}/${DIR1_NAME}_sex_conflict_report.txt"
DIR2_SEX="${DIR2}/${DIR2_NAME}_sex_conflict_report.txt"
COMBINED_SEX="${COMBINED_DIR}/${STUDY_ID}_sex_conflict_report.txt"

echo ""
echo "=== Merging sex conflict report files ==="
echo "  Dir1 source : $DIR1_SEX"
echo "  Dir2 source : $DIR2_SEX"
echo "  Output      : $COMBINED_SEX"

if [[ ! -f "$DIR1_SEX" ]]; then
    echo "  WARNING: Dir1 sex conflict file not found: $DIR1_SEX" >&2
fi
if [[ ! -f "$DIR2_SEX" ]]; then
    echo "  WARNING: Dir2 sex conflict file not found: $DIR2_SEX" >&2
fi

if [[ -f "$DIR1_SEX" && -f "$DIR2_SEX" ]]; then
    {
        head -n 2 "$DIR1_SEX"    | sed 's/\r//g'
        tail -n +3 "$DIR1_SEX"   | sed 's/\r//g' | grep -v "^$"
        tail -n +3 "$DIR2_SEX"   | sed 's/\r//g' | grep -v "^$"
    } > "$COMBINED_SEX"
    echo "  Dir1 rows added : $(tail -n +2 "$DIR1_SEX" | grep -v "^$" | wc -l | tr -d ' ')"
    echo "  Dir2 rows added : $(tail -n +3 "$DIR2_SEX" | grep -v "^$" | wc -l | tr -d ' ')"
    echo "  Total data rows : $(tail -n +3 "$COMBINED_SEX" | wc -l | tr -d ' ')"
elif [[ -f "$DIR1_SEX" ]]; then
    sed 's/\r//g' "$DIR1_SEX" | grep -v "^$" > "$COMBINED_SEX"
    echo "  WARNING: Only Dir1 file found - copied as-is."
elif [[ -f "$DIR2_SEX" ]]; then
    sed 's/\r//g' "$DIR2_SEX" | grep -v "^$" > "$COMBINED_SEX"
    echo "  WARNING: Only Dir2 file found - copied as-is."
else
    echo "  ERROR: Neither sex conflict report file found. Skipping." >&2
fi

echo ""
echo "Step 2 complete. Ready for next step (merge matrices)."

# -------------------------------------------------------
# STEP 3a: Merge gene TPM matrix
# -------------------------------------------------------
# Fixes applied:
#   - sed 's/\r//g' : strips Windows-style CRLF line endings
#   - grep -v "^$"  : drops trailing blank lines
# Matrices are joined on gene/transcript name (column 1) after sorting.
# A count warning is issued if the output row count is less than either input.
# -------------------------------------------------------
DIR1_geneTPM="${DIR1}/${DIR1_NAME}.genes.TPM.matrix"
DIR2_geneTPM="${DIR2}/${DIR2_NAME}.genes.TPM.matrix"
COMBINED_geneTPM="${COMBINED_DIR}/${STUDY_ID}.genes.TPM.matrix"

echo ""
echo "=== Merging gene TPM matrix files ==="
echo "  Dir1 source : $DIR1_geneTPM"
echo "  Dir2 source : $DIR2_geneTPM"
echo "  Output      : $COMBINED_geneTPM"

if [[ ! -f "$DIR1_geneTPM" ]]; then
    echo "  WARNING: Dir1 gene TPM matrix file not found: $DIR1_geneTPM" >&2
fi
if [[ ! -f "$DIR2_geneTPM" ]]; then
    echo "  WARNING: Dir2 gene TPM matrix file not found: $DIR2_geneTPM" >&2
fi

if [[ -f "$DIR1_geneTPM" && -f "$DIR2_geneTPM" ]]; then

    TMP1=$(mktemp)
    TMP2=$(mktemp)
    sed 's/\r//g' "$DIR1_geneTPM" | grep -v "^$" > "$TMP1"
    sed 's/\r//g' "$DIR2_geneTPM" | grep -v "^$" > "$TMP2"

    HEADER1=$(head -n 1 "$TMP1")
    HEADER2=$(head -n 1 "$TMP2" | cut -f2-)
    echo -e "${HEADER1}\t${HEADER2}" > "$COMBINED_geneTPM"

    join -t $'\t' -1 1 -2 1 \
        <(tail -n +2 "$TMP1" | sort -k1,1) \
        <(tail -n +2 "$TMP2" | sort -k1,1) \
        >> "$COMBINED_geneTPM"

    DIR1_GENES=$(tail -n +2 "$TMP1" | wc -l | tr -d ' ')
    DIR2_GENES=$(tail -n +2 "$TMP2" | wc -l | tr -d ' ')
    COMBINED_GENES=$(tail -n +2 "$COMBINED_geneTPM" | wc -l | tr -d ' ')
    echo "  Dir1 genes      : $DIR1_GENES"
    echo "  Dir2 genes      : $DIR2_GENES"
    echo "  Genes in output : $COMBINED_GENES"
    if [[ "$COMBINED_GENES" -lt "$DIR1_GENES" || "$COMBINED_GENES" -lt "$DIR2_GENES" ]]; then
        echo "  WARNING: Gene count in output is less than one or both inputs." >&2
        echo "           Some genes may not have matched between DIR1 and DIR2." >&2
    fi

    rm -f "$TMP1" "$TMP2"

elif [[ -f "$DIR1_geneTPM" ]]; then
    sed 's/\r//g' "$DIR1_geneTPM" | grep -v "^$" > "$COMBINED_geneTPM"
    echo "  WARNING: Only Dir1 file found - copied as-is."
elif [[ -f "$DIR2_geneTPM" ]]; then
    sed 's/\r//g' "$DIR2_geneTPM" | grep -v "^$" > "$COMBINED_geneTPM"
    echo "  WARNING: Only Dir2 file found - copied as-is."
else
    echo "  ERROR: Neither gene TPM matrix file found. Skipping." >&2
fi

echo ""
echo "Step 3a (Combine Gene TPM matrices) complete."

# -------------------------------------------------------
# STEP 3b: Merge gene counts matrix
# -------------------------------------------------------
DIR1_geneCOUNTS="${DIR1}/${DIR1_NAME}.genes.counts.matrix"
DIR2_geneCOUNTS="${DIR2}/${DIR2_NAME}.genes.counts.matrix"
COMBINED_geneCOUNTS="${COMBINED_DIR}/${STUDY_ID}.genes.counts.matrix"

echo ""
echo "=== Merging gene COUNTS matrix files ==="
echo "  Dir1 source : $DIR1_geneCOUNTS"
echo "  Dir2 source : $DIR2_geneCOUNTS"
echo "  Output      : $COMBINED_geneCOUNTS"

if [[ ! -f "$DIR1_geneCOUNTS" ]]; then
    echo "  WARNING: Dir1 gene COUNTS matrix file not found: $DIR1_geneCOUNTS" >&2
fi
if [[ ! -f "$DIR2_geneCOUNTS" ]]; then
    echo "  WARNING: Dir2 gene COUNTS matrix file not found: $DIR2_geneCOUNTS" >&2
fi

if [[ -f "$DIR1_geneCOUNTS" && -f "$DIR2_geneCOUNTS" ]]; then

    TMP1=$(mktemp)
    TMP2=$(mktemp)
    sed 's/\r//g' "$DIR1_geneCOUNTS" | grep -v "^$" > "$TMP1"
    sed 's/\r//g' "$DIR2_geneCOUNTS" | grep -v "^$" > "$TMP2"

    HEADER1=$(head -n 1 "$TMP1")
    HEADER2=$(head -n 1 "$TMP2" | cut -f2-)
    echo -e "${HEADER1}\t${HEADER2}" > "$COMBINED_geneCOUNTS"

    join -t $'\t' -1 1 -2 1 \
        <(tail -n +2 "$TMP1" | sort -k1,1) \
        <(tail -n +2 "$TMP2" | sort -k1,1) \
        >> "$COMBINED_geneCOUNTS"

    DIR1_GENES=$(tail -n +2 "$TMP1" | wc -l | tr -d ' ')
    DIR2_GENES=$(tail -n +2 "$TMP2" | wc -l | tr -d ' ')
    COMBINED_GENES=$(tail -n +2 "$COMBINED_geneCOUNTS" | wc -l | tr -d ' ')
    echo "  Dir1 genes      : $DIR1_GENES"
    echo "  Dir2 genes      : $DIR2_GENES"
    echo "  Genes in output : $COMBINED_GENES"
    if [[ "$COMBINED_GENES" -lt "$DIR1_GENES" || "$COMBINED_GENES" -lt "$DIR2_GENES" ]]; then
        echo "  WARNING: Gene count in output is less than one or both inputs." >&2
        echo "           Some genes may not have matched between DIR1 and DIR2." >&2
    fi

    rm -f "$TMP1" "$TMP2"

elif [[ -f "$DIR1_geneCOUNTS" ]]; then
    sed 's/\r//g' "$DIR1_geneCOUNTS" | grep -v "^$" > "$COMBINED_geneCOUNTS"
    echo "  WARNING: Only Dir1 file found - copied as-is."
elif [[ -f "$DIR2_geneCOUNTS" ]]; then
    sed 's/\r//g' "$DIR2_geneCOUNTS" | grep -v "^$" > "$COMBINED_geneCOUNTS"
    echo "  WARNING: Only Dir2 file found - copied as-is."
else
    echo "  ERROR: Neither gene COUNTS matrix file found. Skipping." >&2
fi

echo ""
echo "Step 3b (Combine Gene COUNTS matrices) complete."

# -------------------------------------------------------
# STEP 3c: Merge transcript TPM matrix
# -------------------------------------------------------
DIR1_txTPM="${DIR1}/${DIR1_NAME}.transcripts.TPM.matrix"
DIR2_txTPM="${DIR2}/${DIR2_NAME}.transcripts.TPM.matrix"
COMBINED_txTPM="${COMBINED_DIR}/${STUDY_ID}.transcripts.TPM.matrix"

echo ""
echo "=== Merging transcript TPM matrix files ==="
echo "  Dir1 source : $DIR1_txTPM"
echo "  Dir2 source : $DIR2_txTPM"
echo "  Output      : $COMBINED_txTPM"

if [[ ! -f "$DIR1_txTPM" ]]; then
    echo "  WARNING: Dir1 transcript TPM matrix file not found: $DIR1_txTPM" >&2
fi
if [[ ! -f "$DIR2_txTPM" ]]; then
    echo "  WARNING: Dir2 transcript TPM matrix file not found: $DIR2_txTPM" >&2
fi

if [[ -f "$DIR1_txTPM" && -f "$DIR2_txTPM" ]]; then

    TMP1=$(mktemp)
    TMP2=$(mktemp)
    sed 's/\r//g' "$DIR1_txTPM" | grep -v "^$" > "$TMP1"
    sed 's/\r//g' "$DIR2_txTPM" | grep -v "^$" > "$TMP2"

    HEADER1=$(head -n 1 "$TMP1")
    HEADER2=$(head -n 1 "$TMP2" | cut -f2-)
    echo -e "${HEADER1}\t${HEADER2}" > "$COMBINED_txTPM"

    join -t $'\t' -1 1 -2 1 \
        <(tail -n +2 "$TMP1" | sort -k1,1) \
        <(tail -n +2 "$TMP2" | sort -k1,1) \
        >> "$COMBINED_txTPM"

    DIR1_TX=$(tail -n +2 "$TMP1" | wc -l | tr -d ' ')
    DIR2_TX=$(tail -n +2 "$TMP2" | wc -l | tr -d ' ')
    COMBINED_TX=$(tail -n +2 "$COMBINED_txTPM" | wc -l | tr -d ' ')
    echo "  Dir1 transcripts      : $DIR1_TX"
    echo "  Dir2 transcripts      : $DIR2_TX"
    echo "  Transcripts in output : $COMBINED_TX"
    if [[ "$COMBINED_TX" -lt "$DIR1_TX" || "$COMBINED_TX" -lt "$DIR2_TX" ]]; then
        echo "  WARNING: Transcript count in output is less than one or both inputs." >&2
        echo "           Some transcripts may not have matched between DIR1 and DIR2." >&2
    fi

    rm -f "$TMP1" "$TMP2"

elif [[ -f "$DIR1_txTPM" ]]; then
    sed 's/\r//g' "$DIR1_txTPM" | grep -v "^$" > "$COMBINED_txTPM"
    echo "  WARNING: Only Dir1 file found - copied as-is."
elif [[ -f "$DIR2_txTPM" ]]; then
    sed 's/\r//g' "$DIR2_txTPM" | grep -v "^$" > "$COMBINED_txTPM"
    echo "  WARNING: Only Dir2 file found - copied as-is."
else
    echo "  ERROR: Neither transcript TPM matrix file found. Skipping." >&2
fi

echo ""
echo "Step 3c (Combine Transcript TPM matrices) complete."

# -------------------------------------------------------
# STEP 3d: Merge transcript counts matrix
# -------------------------------------------------------
DIR1_txCOUNTS="${DIR1}/${DIR1_NAME}.transcripts.counts.matrix"
DIR2_txCOUNTS="${DIR2}/${DIR2_NAME}.transcripts.counts.matrix"
COMBINED_txCOUNTS="${COMBINED_DIR}/${STUDY_ID}.transcripts.counts.matrix"

echo ""
echo "=== Merging transcript counts matrix files ==="
echo "  Dir1 source : $DIR1_txCOUNTS"
echo "  Dir2 source : $DIR2_txCOUNTS"
echo "  Output      : $COMBINED_txCOUNTS"

if [[ ! -f "$DIR1_txCOUNTS" ]]; then
    echo "  WARNING: Dir1 transcript counts matrix file not found: $DIR1_txCOUNTS" >&2
fi
if [[ ! -f "$DIR2_txCOUNTS" ]]; then
    echo "  WARNING: Dir2 transcript counts matrix file not found: $DIR2_txCOUNTS" >&2
fi

if [[ -f "$DIR1_txCOUNTS" && -f "$DIR2_txCOUNTS" ]]; then

    TMP1=$(mktemp)
    TMP2=$(mktemp)
    sed 's/\r//g' "$DIR1_txCOUNTS" | grep -v "^$" > "$TMP1"
    sed 's/\r//g' "$DIR2_txCOUNTS" | grep -v "^$" > "$TMP2"

    HEADER1=$(head -n 1 "$TMP1")
    HEADER2=$(head -n 1 "$TMP2" | cut -f2-)
    echo -e "${HEADER1}\t${HEADER2}" > "$COMBINED_txCOUNTS"

    join -t $'\t' -1 1 -2 1 \
        <(tail -n +2 "$TMP1" | sort -k1,1) \
        <(tail -n +2 "$TMP2" | sort -k1,1) \
        >> "$COMBINED_txCOUNTS"

    DIR1_TX=$(tail -n +2 "$TMP1" | wc -l | tr -d ' ')
    DIR2_TX=$(tail -n +2 "$TMP2" | wc -l | tr -d ' ')
    COMBINED_TX=$(tail -n +2 "$COMBINED_txCOUNTS" | wc -l | tr -d ' ')
    echo "  Dir1 transcripts      : $DIR1_TX"
    echo "  Dir2 transcripts      : $DIR2_TX"
    echo "  Transcripts in output : $COMBINED_TX"
    if [[ "$COMBINED_TX" -lt "$DIR1_TX" || "$COMBINED_TX" -lt "$DIR2_TX" ]]; then
        echo "  WARNING: Transcript count in output is less than one or both inputs." >&2
        echo "           Some transcripts may not have matched between DIR1 and DIR2." >&2
    fi

    rm -f "$TMP1" "$TMP2"

elif [[ -f "$DIR1_txCOUNTS" ]]; then
    sed 's/\r//g' "$DIR1_txCOUNTS" | grep -v "^$" > "$COMBINED_txCOUNTS"
    echo "  WARNING: Only Dir1 file found - copied as-is."
elif [[ -f "$DIR2_txCOUNTS" ]]; then
    sed 's/\r//g' "$DIR2_txCOUNTS" | grep -v "^$" > "$COMBINED_txCOUNTS"
    echo "  WARNING: Only Dir2 file found - copied as-is."
else
    echo "  ERROR: Neither transcript counts matrix file found. Skipping." >&2
fi

echo ""
echo "Step 3d (Combine Transcript COUNTS matrices) complete."

# -------------------------------------------------------
# Step 4: Copy log_files/ subdirectories, namespaced by source dir
# -------------------------------------------------------
DIR1_logDIR="${DIR1}/log_files"
DIR2_logDIR="${DIR2}/log_files"
COMBINED_logDIR="${COMBINED_DIR}/log_files"
echo ""
echo "=== Transferring log_files directories ==="
echo "  Dir1 source : $DIR1_logDIR"
echo "  Dir2 source : $DIR2_logDIR"
echo "  Output      : $COMBINED_logDIR"

mkdir -p "$COMBINED_logDIR"

if [[ -d "$DIR1_logDIR" ]]; then
    cp -R "$DIR1_logDIR" "$COMBINED_logDIR/log_files_${DIR1_NAME}"
else
    echo "WARNING: Dir1 log directory not found: $DIR1_logDIR" >&2
fi

if [[ -d "$DIR2_logDIR" ]]; then
    cp -R "$DIR2_logDIR" "$COMBINED_logDIR/log_files_${DIR2_NAME}"
else
    echo "WARNING: Dir2 log directory not found: $DIR2_logDIR" >&2
fi
echo ""
echo "Step 4 (Transfer log_files directory) complete."

# -------------------------------------------------------
# Step 5: Copy reads_fastq/ sample subdirectories (sample logs, bigwigs, JSONs)
# Duplicate sample directory names across DIR1 and DIR2 are detected and reported.
# -------------------------------------------------------
DIR1_readsDIR="${DIR1}/reads_fastq"
DIR2_readsDIR="${DIR2}/reads_fastq"
COMBINED_readsDIR="${COMBINED_DIR}/reads_fastq"

echo ""
echo "=== Transferring reads_fastq directories ==="
echo "  Dir1 source : $DIR1_readsDIR"
echo "  Dir2 source : $DIR2_readsDIR"
echo "  Output      : $COMBINED_readsDIR"

mkdir -p "$COMBINED_readsDIR"

DIR1_count=0
DIR2_count=0
DUP_COUNT=0
TMP_DUP1=$(mktemp)
TMP_DUP2=$(mktemp)

if [[ -d "$DIR1_readsDIR" ]]; then
    DIR1_count=$(find "$DIR1_readsDIR" -mindepth 1 -maxdepth 1 -type d | wc -l | tr -d ' ')
    find "$DIR1_readsDIR" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort > "$TMP_DUP1"
    cp -R "$DIR1_readsDIR/." "$COMBINED_readsDIR"
else
    echo "WARNING: Dir1 reads_fastq directory not found: $DIR1_readsDIR" >&2
fi

if [[ -d "$DIR2_readsDIR" ]]; then
    DIR2_count=$(find "$DIR2_readsDIR" -mindepth 1 -maxdepth 1 -type d | wc -l | tr -d ' ')
    find "$DIR2_readsDIR" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort > "$TMP_DUP2"
    cp -R "$DIR2_readsDIR/." "$COMBINED_readsDIR"
else
    echo "WARNING: Dir2 reads_fastq directory not found: $DIR2_readsDIR" >&2
fi

if [[ -s "$TMP_DUP1" && -s "$TMP_DUP2" ]]; then
    DUPLICATES=$(comm -12 "$TMP_DUP1" "$TMP_DUP2" || true)
    if [[ -n "$DUPLICATES" ]]; then
        DUP_COUNT=$(echo "$DUPLICATES" | grep -c . || true)
        echo "WARNING: Duplicate sample directory names found in DIR1 and DIR2:" >&2
        echo "$DUPLICATES" >&2
    fi
fi

rm -f "$TMP_DUP1" "$TMP_DUP2"

echo ""
echo "Step 5 (Copy Sample data) complete."
echo "  Sample directories copied from $DIR1_NAME : $DIR1_count"
echo "  Sample directories copied from $DIR2_NAME : $DIR2_count"
echo "  Duplicate sample directory names detected : $DUP_COUNT"
echo "  Total input sample directories             : $((DIR1_count + DIR2_count))"

# -------------------------------------------------------
# Step 6: Copy final MultiQC HTML reports
# -------------------------------------------------------
DIR1_MultiQC="${DIR1}/${DIR1_NAME}_final_multiQC_report.html"
DIR2_MultiQC="${DIR2}/${DIR2_NAME}_final_multiQC_report.html"

echo ""
echo "=== Transferring final MultiQC html reports ==="
echo "  Dir1 source : $DIR1_MultiQC"
echo "  Dir2 source : $DIR2_MultiQC"
echo "  Output Directory    : $COMBINED_DIR"

if [[ -f "$DIR1_MultiQC" ]]; then
    cp -R "$DIR1_MultiQC" "$COMBINED_DIR"
else
    echo "WARNING: Dir1 final MultiQC report not found: $DIR1_MultiQC" >&2
fi

if [[ -f "$DIR2_MultiQC" ]]; then
    cp -R "$DIR2_MultiQC" "$COMBINED_DIR"
else
    echo "WARNING: Dir2 final MultiQC report not found: $DIR2_MultiQC" >&2
fi
echo ""
echo "Step 6 (Copy final multiQC reports) complete."

# -------------------------------------------------------
# Step 7: Regenerate JBrowse2 session JSON
# -------------------------------------------------------
COMBINED_SESSION_JSON="${COMBINED_DIR}/${STUDY_ID}_jbrowse_session_GRCr8.json"

echo ""
echo "=== Building combined JBrowse session JSON ==="
echo "  BIOProjectID : $STUDY_ID"
echo "  PRJdir       : $COMBINED_DIR"
echo "  baseDir      : $COMBINED_DIR"
echo "  Script       : $JBROWSE_PY"
echo "  Output       : $COMBINED_SESSION_JSON"

if [[ ! -f "$JBROWSE_PY" ]]; then
    echo "ERROR: JBrowse session builder script not found: $JBROWSE_PY" >&2
    exit 1
fi

# Load Python — adjust the module name or remove this line if not using Lmod
module load python/3.12.10 2>/dev/null || true

python "$JBROWSE_PY" \
    "$STUDY_ID" \
    "$COMBINED_DIR" \
    "$COMBINED_DIR"

if [[ -f "$COMBINED_SESSION_JSON" ]]; then
    TRACK_JSON_COUNT=$(find "$COMBINED_DIR" -type f -name 'RNAseq_*.json' \
        ! -name '*geneTPMbed.json' \
        ! -name '*TXTPMbed.json' | wc -l | tr -d ' ')
    echo "Step 7 (Build combined JBrowse session JSON) complete."
    echo "  RNAseq track JSON files detected : $TRACK_JSON_COUNT"
    echo "  Session JSON written to          : $COMBINED_SESSION_JSON"
else
    echo "ERROR: Combined JBrowse session JSON was not created: $COMBINED_SESSION_JSON" >&2
    exit 1
fi
