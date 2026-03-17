#!/bin/bash
#SBATCH --job-name=SRA2QC
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=1gb
#SBATCH --time=06:00:00
#SBATCH --account=your-slurm-account
#SBATCH --partition=normal
#SBATCH --output=%x-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=your@email.edu
#SBATCH --error=SRA2QC.err

###############################################################################
# USER CONFIGURATION
# Update the variables below for your system before running.
###############################################################################
# Directory containing all pipeline scripts — screenconfig is set from this
SCRIPT_DIR="/path/to/your/pipeline/scripts"
# Scratch filesystem root; scratch_dir is built as SCRATCH_BASE/BIOProjectID
SCRATCH_BASE="/path/to/your/scratch/mount"
###############################################################################


# Script details
## Downloads SRA files, converts to FASTQ, runs QC
## Outputs to scratch directory: /path/to/your/scratch/$BIOProjectID/$Run/

set -x

# ---------------------------------------------------------------------------
# Retry configuration — tune these without touching the logic below
# ---------------------------------------------------------------------------
PREFETCH_MAX_ATTEMPTS=8       # number of retry attempts
PREFETCH_RETRY_WAIT=120       # base wait (seconds) between retries — doubles each failure
PREFETCH_MAX_SIZE="100G"      # covers large SRA files
FASTERQ_MAX_ATTEMPTS=3        # how many times to retry fasterq-dump
FASTERQ_RETRY_WAIT=10         # seconds to wait between fasterq-dump retries

# Load required modules
module load sratoolkit/3.1.1 fastq-screen/0.15.2 fastqc/0.11.9 bowtie/2.5 python/3.9.1

# Retrieve environment variables
Run=${Run}
geo_accession=${geo_accession}
PRJdir=${PRJdir}
screenconfig="$SCRIPT_DIR/dependencies/FastQ_Screen_Genomes/fastq_screen.conf"
BIOProjectID=${BIOProjectID}

echo "============================================"
echo "SRA2QC Processing Started: $(date "+%Y-%m-%d %H:%M:%S")"
echo "============================================"
echo "Run:           $Run"
echo "geo_accession: $geo_accession"
echo "PRJdir:        $PRJdir"
echo "screenconfig:  $screenconfig"
echo "BIOProjectID:  $BIOProjectID"


mkdir -p "$scratch_dir/$Run"
echo "Scratch directory: $scratch_dir/$Run"
echo ""

# ---------------------------------------------------------------------------
# Helper: remove incomplete/corrupt prefetch output so retries start clean
# ---------------------------------------------------------------------------
clean_prefetch_output() {
    local sra_dir="$scratch_dir/$Run"
    local sra_file="$sra_dir/${Run}.sra"
    local sra_file_lite="$sra_dir/${Run}.sra.lite"
    # Remove partial downloads
    rm -f "$sra_file" "$sra_file_lite"
    # sratoolkit sometimes leaves a lock or partial sub-directory
    rm -rf "$scratch_dir/${Run}.tmp" "$scratch_dir/${Run}.lock"
}

echo "============================================"
echo "Step 1: Prefetching SRA file (up to ${PREFETCH_MAX_ATTEMPTS} attempts)..."
echo "============================================"

# Set the SRA cache root to scratch so prefetch writes directly there.
# This also prevents it from defaulting to home directory and filling quota.
vdb-config -s "/repository/user/main/public/root=$scratch_dir" 2>/dev/null || true

prefetch_attempt=1
prefetch_success=false
wait_time=$PREFETCH_RETRY_WAIT

while [ $prefetch_attempt -le $PREFETCH_MAX_ATTEMPTS ]; do
    echo "Prefetch attempt $prefetch_attempt of $PREFETCH_MAX_ATTEMPTS for $Run"

    # --resume yes: picks up a partial download instead of restarting from zero.
    #               Critical when HTTPS drops mid-file after 20+ minutes.
    # --max-size:   raised to 100G to prevent silent failures on large files.
    # --progress:   surfaces download progress in the log.
    if prefetch \
            --max-size "$PREFETCH_MAX_SIZE" \
            --progress \
            --resume yes \
            --output-directory "$scratch_dir" \
            "$Run"; then

        # Confirm the file actually landed — prefetch can exit 0 spuriously
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
        # Exponential backoff: double the wait each time, cap at 10 minutes
        wait_time=$(( wait_time * 2 ))
        [ "$wait_time" -gt 600 ] && wait_time=600
    fi
done

if [ "$prefetch_success" = false ]; then
    echo "============================================"
    echo "FATAL ERROR: prefetch failed for $Run after $PREFETCH_MAX_ATTEMPTS attempts"
    echo "  Transports tried: ${TRANSPORT_CYCLE[*]:0:$PREFETCH_MAX_ATTEMPTS}"
    echo "============================================"
    exit 1
fi
echo ""

echo "============================================"
echo "Step 2: Validating SRA file..."
echo "============================================"
# Point vdb-validate at the actual file, not just the directory
sra_file_path=$(find "$scratch_dir/$Run" -name "${Run}.sra" -o -name "${Run}.sralite" 2>/dev/null | head -1)
vdb-validate "$sra_file_path" || {
    echo "ERROR: Validation failed for $Run";
    exit 1;
}
echo "Validation completed successfully"
echo ""

echo "============================================"
echo "Step 3: Converting SRA to FASTQ..."
echo "============================================"

fasterq_attempt=1
fasterq_success=false
new_fastq1=""
new_fastq2=""
new_fastq3=""

while [ $fasterq_attempt -le $FASTERQ_MAX_ATTEMPTS ]; do
    echo "fasterq-dump attempt $fasterq_attempt of $FASTERQ_MAX_ATTEMPTS for $Run"

    if fasterq-dump.3.1.1 "$scratch_dir/$Run" \
            --threads 12 \
            --outdir "$scratch_dir/$Run" \
            -t "$scratch_dir/$Run" \
            -x; then

        fastq1="${scratch_dir}/${Run}/${Run}_1.fastq"
        fastq2="${scratch_dir}/${Run}/${Run}_2.fastq"
        fastq_se="${scratch_dir}/${Run}/${Run}.fastq"

        # ---------------------------------------------------------------
        # Layout detection
        # ---------------------------------------------------------------
        if [ -f "$fastq1" ] && [ -f "$fastq2" ]; then
            echo "Paired-end FASTQ files detected — proceeding"

            new_fastq1="${scratch_dir}/${Run}/${geo_accession}_${Run}_1.fastq"
            new_fastq2="${scratch_dir}/${Run}/${geo_accession}_${Run}_2.fastq"
            mv "$fastq1" "$new_fastq1"
            mv "$fastq2" "$new_fastq2"
            echo "  Renamed: $(basename "$new_fastq1")"
            echo "  Renamed: $(basename "$new_fastq2")"

            # Unpaired remainder file (present in some paired runs)
            if [ -f "$fastq_se" ]; then
                new_fastq3="${scratch_dir}/${Run}/${geo_accession}_${Run}.fastq"
                mv "$fastq_se" "$new_fastq3"
                echo "  Renamed unpaired remainder: $(basename "$new_fastq3")"
            fi

            fasterq_success=true
            break

        elif [ -f "$fastq_se" ]; then
            # ---------------------------------------------------------------
            # Single-end detected — this pipeline does not process single-end.
            # Exit with code 2 so the caller can distinguish this from a true
            # failure (exit 1) and route the sample to the SE pipeline instead.
            # ---------------------------------------------------------------
            echo "============================================"
            echo "SINGLE-END LAYOUT DETECTED — NOT PROCESSED"
            echo "============================================"
            echo "Run:          $Run"
            echo "geo_accession: $geo_accession"
            echo "BIOProjectID: $BIOProjectID"
            echo ""
            echo "fasterq-dump produced a single-end FASTQ:"
            echo "  $fastq_se"
            echo ""
            echo "This script handles paired-end data only."
            echo "Please resubmit $Run through the single-end pipeline."
            echo "Exiting with code 2 to flag as SE, not a pipeline error."
            echo "============================================"
            exit 2

        else
            echo "ERROR: No FASTQ output found for $Run on attempt $fasterq_attempt"
            echo "Contents of output directory:"
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
fastqc "$new_fastq1" "$new_fastq2" --outdir "$scratch_dir/${Run}/" || {
    echo "ERROR: FastQC failed"; exit 1;
}
echo "FastQC completed successfully"
echo ""

echo "============================================"
echo "Step 5: Running FastQ-Screen..."
echo "============================================"
fastq_screen --outdir "$scratch_dir/${Run}/" --conf "$screenconfig" "$new_fastq1" "$new_fastq2" || {
    echo "ERROR: FastQ-Screen failed"; exit 1;
}
echo "FastQ-Screen completed successfully"
echo ""

echo "============================================"
echo "Step 6: Compressing FASTQ files..."
echo "============================================"
pigz "$new_fastq1" "$new_fastq2"
[ -f "${new_fastq1}.gz" ] && [ -f "${new_fastq2}.gz" ] || {
    echo "ERROR: FASTQ compression failed"; exit 1;
}
echo "  Created: $(basename "${new_fastq1}.gz")"
echo "  Created: $(basename "${new_fastq2}.gz")"
if [ -n "$new_fastq3" ] && [ -f "$new_fastq3" ]; then
    pigz "$new_fastq3"
    [ -f "${new_fastq3}.gz" ] || { echo "ERROR: Unpaired remainder compression failed"; exit 1; }
    echo "  Created: $(basename "${new_fastq3}.gz") (unpaired remainder)"
fi
echo "Compression completed successfully"
echo ""

echo "============================================"
echo "Step 7: Cleaning up SRA file..."
echo "============================================"
sra_file="${scratch_dir}/${Run}/${Run}.sra"
sra_file_lite="${scratch_dir}/${Run}/${Run}.sralite"
for f in "$sra_file" "$sra_file_lite"; do
    if [ -f "$f" ]; then
        rm "$f"
        echo "Removed: $(basename "$f")"
    fi
done
echo ""

echo "============================================"
echo "SRA2QC Processing Complete!"
echo "============================================"
echo "Sample:  $geo_accession  (Run: $Run)"
echo "Layout:  paired-end"
echo "Out dir: $scratch_dir/$Run/"
echo "Files created:"
echo "  - $(basename "${new_fastq1}.gz")"
echo "  - $(basename "${new_fastq2}.gz")"
[ -n "$new_fastq3" ] && [ -f "${new_fastq3}.gz" ] && \
    echo "  - $(basename "${new_fastq3}.gz") (unpaired remainder)"
echo "  - FastQC HTML reports"
echo "  - FastQ-Screen HTML reports"
echo "Completed: $(date "+%Y-%m-%d %H:%M:%S")"
echo "============================================"

exit 0