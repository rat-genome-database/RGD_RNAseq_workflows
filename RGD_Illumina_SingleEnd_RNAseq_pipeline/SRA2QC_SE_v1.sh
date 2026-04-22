#!/usr/bin/env bash
#SBATCH --job-name=SRA2QC_SE
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=1gb
#SBATCH --time=06:00:00
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=normal
#SBATCH --output=%x-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=your@email.edu
#SBATCH --error=SRA2QC_SE.err

###############################################################################
# USER CONFIGURATION
# Update the variables below for your system before running.
###############################################################################
# Directory containing all pipeline scripts
myDir="/path/to/home"
SCRIPT_DIR="$myDir/STAR_RSEM_pipeline/RGD_Illumina_SingleEnd_RNAseq_pipeline"
# Scratch filesystem root; scratch_dir is built as SCRATCH_BASE/BIOProjectID
SCRATCH_BASE="/path/to/scratch"
###############################################################################

# Script details
## Downloads SRA files, converts to single-end FASTQ, runs QC
## Outputs to scratch directory: SCRATCH_BASE/BIOProjectID/Run/

set -x

# ---------------------------------------------------------------------------
# Retry configuration
# ---------------------------------------------------------------------------
PREFETCH_MAX_ATTEMPTS=8
PREFETCH_RETRY_WAIT=120
PREFETCH_MAX_SIZE="100G"
FASTERQ_MAX_ATTEMPTS=3
FASTERQ_RETRY_WAIT=10

# Load required modules
module load sratoolkit/3.1.1 fastq-screen/0.15.2 fastqc/0.11.9 bowtie/2.5 python/3.9.1

# Retrieve environment variables
Run=${Run}
geo_accession=${geo_accession}
PRJdir=${PRJdir}
BIOProjectID=${BIOProjectID}

echo "============================================"
echo "SRA2QC_SE Processing Started: $(date "+%Y-%m-%d %H:%M:%S")"
echo "============================================"
echo "Run:           $Run"
echo "geo_accession: $geo_accession"
echo "PRJdir:        $PRJdir"
echo "screenconfig:  $screenconfig"
echo "BIOProjectID:  $BIOProjectID"

scratch_dir="${SCRATCH_BASE}/${BIOProjectID}"

mkdir -p "$scratch_dir/$Run"
echo "Scratch directory: $scratch_dir/$Run"
echo ""

# ---------------------------------------------------------------------------
# Helper: remove incomplete prefetch output
# ---------------------------------------------------------------------------
clean_prefetch_output() {
    local sra_dir="$scratch_dir/$Run"
    rm -f "$sra_dir/${Run}.sra" "$sra_dir/${Run}.sra.lite"
    rm -rf "$scratch_dir/${Run}.tmp" "$scratch_dir/${Run}.lock"
}

echo "============================================"
echo "Step 1: Prefetching SRA file (up to ${PREFETCH_MAX_ATTEMPTS} attempts)..."
echo "============================================"

vdb-config -s "/repository/user/main/public/root=$scratch_dir" 2>/dev/null || true

prefetch_attempt=1
prefetch_success=false
wait_time=$PREFETCH_RETRY_WAIT

while [ $prefetch_attempt -le $PREFETCH_MAX_ATTEMPTS ]; do
    echo "Prefetch attempt $prefetch_attempt of $PREFETCH_MAX_ATTEMPTS for $Run"

    if prefetch \
            --max-size "$PREFETCH_MAX_SIZE" \
            --progress \
            --resume yes \
            --output-directory "$scratch_dir" \
            "$Run"; then

        sra_path=$(find "$scratch_dir/$Run" \
            \( -name "${Run}.sra" -o -name "${Run}.sralite" \) \
            2>/dev/null | head -1)

        if [ -n "$sra_path" ]; then
            echo "Prefetch succeeded on attempt $prefetch_attempt"
            echo "  File: $sra_path"
            prefetch_success=true
            break
        else
            echo "WARNING: prefetch exited 0 but no .sra file found — treating as failure"
            clean_prefetch_output
        fi
    else
        echo "ERROR: prefetch failed on attempt $prefetch_attempt for $Run"
    fi

    prefetch_attempt=$((prefetch_attempt + 1))
    if [ $prefetch_attempt -le $PREFETCH_MAX_ATTEMPTS ]; then
        echo "Waiting ${wait_time}s before next attempt..."
        sleep "$wait_time"
        wait_time=$(( wait_time * 2 ))
        [ "$wait_time" -gt 600 ] && wait_time=600
    fi
done

if [ "$prefetch_success" = false ]; then
    echo "============================================"
    echo "FATAL ERROR: prefetch failed for $Run after $PREFETCH_MAX_ATTEMPTS attempts"
    echo "============================================"
    exit 1
fi

echo ""
echo "============================================"
echo "Step 2: Validating SRA file..."
echo "============================================"
sra_file_path=$(find "$scratch_dir/$Run" -name "${Run}.sra" -o -name "${Run}.sralite" 2>/dev/null | head -1)
vdb-validate "$sra_file_path" || {
    echo "ERROR: Validation failed for $Run"
    exit 1
}
echo "Validation completed successfully"
echo ""

echo "============================================"
echo "Step 3: Converting SRA to FASTQ (single-end)..."
echo "============================================"

fasterq_attempt=1
fasterq_success=false
new_fastq=""

while [ $fasterq_attempt -le $FASTERQ_MAX_ATTEMPTS ]; do
    echo "fasterq-dump attempt $fasterq_attempt of $FASTERQ_MAX_ATTEMPTS for $Run"

    if fasterq-dump.3.1.1 "$scratch_dir/$Run" \
            --threads 12 \
            --outdir "$scratch_dir/$Run" \
            -t "$scratch_dir/$Run"; then

        fastq_se="${scratch_dir}/${Run}/${Run}.fastq"
        fastq1="${scratch_dir}/${Run}/${Run}_1.fastq"

        if [ -f "$fastq_se" ]; then
            echo "Single-end FASTQ detected"
            new_fastq="${scratch_dir}/${Run}/${geo_accession}_${Run}.fastq"
            mv "$fastq_se" "$new_fastq"
            echo "  Renamed: $(basename "$new_fastq")"
            fasterq_success=true
            break

        elif [ -f "$fastq1" ]; then
            # Paired-end data encountered — exit with code 2 so caller can route to PE pipeline
            echo "============================================"
            echo "PAIRED-END LAYOUT DETECTED — NOT PROCESSED"
            echo "============================================"
            echo "Run:           $Run"
            echo "geo_accession: $geo_accession"
            echo "BIOProjectID:  $BIOProjectID"
            echo ""
            echo "fasterq-dump produced paired-end FASTQs for a single-end pipeline submission."
            echo "Please resubmit $Run through the paired-end pipeline."
            echo "Exiting with code 2 to flag as PE, not a pipeline error."
            echo "============================================"
            exit 2

        else
            echo "ERROR: No FASTQ output found for $Run on attempt $fasterq_attempt"
            ls -lh "$scratch_dir/${Run}/" || true
        fi
    else
        echo "ERROR: fasterq-dump failed for $Run on attempt $fasterq_attempt"
    fi

    fasterq_attempt=$((fasterq_attempt + 1))
    if [ $fasterq_attempt -le $FASTERQ_MAX_ATTEMPTS ]; then
        echo "Retrying in ${FASTERQ_RETRY_WAIT}s..."
        sleep "$FASTERQ_RETRY_WAIT"
    fi
done

if [ "$fasterq_success" = false ]; then
    echo "============================================"
    echo "FATAL ERROR: fasterq-dump failed for $Run after $FASTERQ_MAX_ATTEMPTS attempts"
    echo "============================================"
    exit 1
fi
echo ""

echo "============================================"
echo "Step 4: Running FastQC..."
echo "============================================"
fastqc "$new_fastq" --outdir "$scratch_dir/${Run}/" || {
    echo "ERROR: FastQC failed"
    exit 1
}
echo "FastQC completed successfully"
echo ""

echo "============================================"
echo "Step 5: Running FastQ-Screen..."
echo "============================================"
fastq_screen --outdir "$scratch_dir/${Run}/" --conf "$screenconfig" "$new_fastq" || {
    echo "ERROR: FastQ-Screen failed"
    exit 1
}
echo "FastQ-Screen completed successfully"
echo ""

echo "============================================"
echo "Step 6: Compressing FASTQ file..."
echo "============================================"
pigz "$new_fastq"
[ -f "${new_fastq}.gz" ] || {
    echo "ERROR: FASTQ compression failed"
    exit 1
}
echo "  Created: $(basename "${new_fastq}.gz")"
echo "Compression completed successfully"
echo ""

echo "============================================"
echo "Step 7: Cleaning up SRA file..."
echo "============================================"
for f in "${scratch_dir}/${Run}/${Run}.sra" "${scratch_dir}/${Run}/${Run}.sralite"; do
    if [ -f "$f" ]; then
        rm "$f"
        echo "Removed: $(basename "$f")"
    fi
done
echo ""

echo "============================================"
echo "SRA2QC_SE Processing Complete!"
echo "============================================"
echo "Sample:  $geo_accession  (Run: $Run)"
echo "Layout:  single-end"
echo "Out dir: $scratch_dir/$Run/"
echo "Files created:"
echo "  - $(basename "${new_fastq}.gz")"
echo "  - FastQC HTML report"
echo "  - FastQ-Screen HTML report"
echo "Completed: $(date "+%Y-%m-%d %H:%M:%S")"
echo "============================================"

exit 0
