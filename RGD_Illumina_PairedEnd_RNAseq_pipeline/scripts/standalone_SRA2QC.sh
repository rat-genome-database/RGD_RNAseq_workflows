#!/bin/bash
#SBATCH --job-name=sSRA2QC
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=1gb
#SBATCH --time=01:35:00
#SBATCH --account=your-slurm-account
#SBATCH --partition=normal
#SBATCH --output=%x-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=your@email.edu
#SBATCH --error=ssSRA2QC.err

###############################################################################
# USER CONFIGURATION
# Update the variables below for your system before running.
###############################################################################
# Directory containing all pipeline scripts
SCRIPT_DIR="/path/to/your/pipeline/scripts"
# Scratch filesystem root; scratch_dir is built as SCRATCH_BASE/BIOProjectID
SCRATCH_BASE="/path/to/your/scratch/mount"
###############################################################################


# Script details
## Downloads SRA files, converts to FASTQ, runs QC
##  Use to retry sample retrieval and  conversion
## Usage sh standalone_SRA2QC.sh $BIOProjectID $Run $geo_accession

# Enable debugging by printing each command before execution
set -x

# Load required modules
module load sratoolkit/3.1.1 fastq-screen/0.15.2 fastqc/0.11.9 bowtie/2.5 python/3.9.1

# Retrieve environment variables

BIOProjectID=${1}
Run=${2}
geo_accession=${3}


echo "============================================"
echo "SRA2QC Processing Started: $(date "+%Y-%m-%d %H:%M:%S")"
echo "============================================"
echo "Run: $Run"
echo "geo_accession: $geo_accession"
echo "PRJdir: $PRJdir"
echo "screenconfig: $screenconfig"
echo "BIOProjectID: $BIOProjectID"


# Create scratch directory for this run
mkdir -p "$scratch_dir/$Run"
echo "Scratch directory: $scratch_dir/$Run"
echo ""

echo "============================================"
echo "Step 1: Prefetching SRA file..."
echo "============================================"
prefetch -O "$scratch_dir" "$Run" --max-size 30G || {
    echo "ERROR: Prefetch failed for $Run";
    exit 1;
}
echo "Prefetch completed successfully"
echo ""

echo "============================================"
echo "Step 2: Validating SRA file..."
echo "============================================"
vdb-validate "$scratch_dir/$Run" || {
    echo "ERROR: Validation failed for $Run";
    exit 1;
}
echo "Validation completed successfully"
echo ""

echo "============================================"
echo "Step 3: Converting SRA to FASTQ..."
echo "============================================"

# Attempt to run fasterq-dump up to 3 times if needed
attempt=1
max_attempts=3
success=false

while [ $attempt -le $max_attempts ]; do
    echo "Attempt $attempt of $max_attempts: Running fasterq-dump for $Run"

    if fasterq-dump.3.1.1 "$scratch_dir/$Run" \
        --threads 12 \
        --outdir "$scratch_dir/$Run" \
        -t "${SCRATCH_BASE}/${BIOProjectID}" \
        -x; then

        # Define expected FASTQ file paths
        fastq1="${scratch_dir}/${Run}/${Run}_1.fastq"
        fastq2="${scratch_dir}/${Run}/${Run}_2.fastq"
        fastq3="${scratch_dir}/${Run}/${Run}.fastq"

        # Check if paired-end files were created
        if [ -f "$fastq1" ] && [ -f "$fastq2" ]; then
            echo "Paired-end FASTQ files found"

            # Rename files to include geo_accession
            echo "Renaming FASTQ files to include GEO accession: $geo_accession"
            new_fastq1="${scratch_dir}/${Run}/${geo_accession}_${Run}_1.fastq"
            new_fastq2="${scratch_dir}/${Run}/${geo_accession}_${Run}_2.fastq"
            mv "$fastq1" "$new_fastq1"
            mv "$fastq2" "$new_fastq2"
            echo "  - Renamed to: $(basename "$new_fastq1")"
            echo "  - Renamed to: $(basename "$new_fastq2")"

            # Handle unpaired reads if they exist
            if [ -f "$fastq3" ]; then
                new_fastq3="${scratch_dir}/${Run}/${geo_accession}_${Run}.fastq"
                mv "$fastq3" "$new_fastq3"
                echo "  - Renamed unpaired reads to: $(basename "$new_fastq3")"
            fi
            echo ""

            echo "============================================"
            echo "Step 4: Running FastQC..."
            echo "============================================"
            fastqc "$new_fastq1" "$new_fastq2" --outdir "$scratch_dir/${Run}/" || {
                echo "ERROR: FastQC failed";
                exit 1;
            }
            echo "FastQC completed successfully"
            echo ""

            echo "============================================"
            echo "Step 5: Running FastQ-Screen..."
            echo "============================================"
            fastq_screen --outdir "$scratch_dir/${Run}/" --conf "$screenconfig" "$new_fastq1" "$new_fastq2" || {
                echo "ERROR: FastQ-Screen failed";
                exit 1;
            }
            echo "FastQ-Screen completed successfully"
            echo ""

            echo "============================================"
            echo "Step 6: Compressing FASTQ files..."
            echo "============================================"
            echo "Compressing paired-end FASTQ files..."
            pigz "$new_fastq1" "$new_fastq2"

            # Verify compression succeeded
            if [ ! -f "${new_fastq1}.gz" ] || [ ! -f "${new_fastq2}.gz" ]; then
                echo "ERROR: FASTQ compression failed"
                exit 1
            fi

            echo "  - Created: $(basename "${new_fastq1}.gz")"
            echo "  - Created: $(basename "${new_fastq2}.gz")"

            # Compress unpaired reads if they exist
            if [ -f "$new_fastq3" ]; then
                echo "Compressing unpaired FASTQ file..."
                pigz "$new_fastq3"

                if [ ! -f "${new_fastq3}.gz" ]; then
                    echo "ERROR: Unpaired FASTQ compression failed"
                    exit 1
                fi

                echo "  - Created: $(basename "${new_fastq3}.gz")"
            fi
            echo "Compression completed successfully"
            echo ""

            echo "============================================"
            echo "Step 7: Cleaning up SRA file..."
            echo "============================================"
            sra_file="${scratch_dir}/${Run}/${Run}.sra"
            if [ -f "$sra_file" ]; then
                rm "$sra_file"
                echo "Removed SRA file: $(basename "$sra_file")"
            else
                echo "No SRA file to remove (already cleaned up)"
            fi
            echo ""

            echo "============================================"
            echo "SRA2QC Processing Complete!"
            echo "============================================"
            echo "Sample: $geo_accession (Run: $Run)"
            echo "Output directory: $scratch_dir/$Run/"
            echo "Files created:"
            echo "  - ${geo_accession}_${Run}_1.fastq.gz"
            echo "  - ${geo_accession}_${Run}_2.fastq.gz"
            if [ -f "${new_fastq3}.gz" ]; then
                echo "  - ${geo_accession}_${Run}.fastq.gz (unpaired)"
            fi
            echo "  - FastQC HTML reports"
            echo "  - FastQ-Screen HTML reports"
            echo "Completed: $(date "+%Y-%m-%d %H:%M:%S")"
            echo "============================================"

            success=true
            exit 0

        else
            echo "ERROR: Paired-end FASTQ files not found for $Run on attempt $attempt"
            echo "Expected files:"
            echo "  - $fastq1"
            echo "  - $fastq2"
            ls -lh "$scratch_dir/${Run}/" || true
        fi
    else
        echo "ERROR: fasterq-dump failed for $Run on attempt $attempt"
    fi

    attempt=$((attempt + 1))

    if [ $attempt -le $max_attempts ]; then
        echo "Retrying in 10 seconds..."
        sleep 10
    fi
done

if [ "$success" = false ]; then
    echo "============================================"
    echo "FATAL ERROR: Failed to process $Run after $max_attempts attempts"
    echo "============================================"
    exit 1
fi
