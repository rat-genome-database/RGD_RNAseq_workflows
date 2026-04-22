#!/usr/bin/env bash
#SBATCH --job-name=computeSex
#SBATCH --ntasks=1
#SBATCH --mem=1gb
#SBATCH --time=00:10:00
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=normal
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=your@email.edu

###############################################################################
# USER CONFIGURATION
# Update the variables below for your system before running.
###############################################################################
# Your home/base directory
myDir="/path/to/home"
# Scratch filesystem root; scratch_dir is built as SCRATCH_BASE/BIOProjectID
SCRATCH_BASE="/path/to/scratch"
###############################################################################


# Script details
## Written by Wendy Demos
## Last Modified 31 March 2026 with Claude AI assistence by WMD to run both single read and multiple reads per sample
## Updated Feb 2026 (v5): Fixed BAM paths to match STAR output directory structure
##                        Now processes deduplicated PASS samples only
## Updated 7 April 2026 (v6): ensure final output is written to baseDir,
##                             keep temp work in scratch, and reconstruct
##                             scratch_dir if not exported by caller

# Enable debugging by printing each command before execution
set -x
set -euo pipefail

module load samtools/1.20 

# Retrieve environment variables
BIOProjectID=${BIOProjectID}
PRJdir=${PRJdir}
baseDir=${baseDir}
Logdir=${Logdir}
AccList=${AccList}
scratch_dir=${scratch_dir:-${SCRATCH_BASE}/${BIOProjectID}}

echo "============================================"
echo "ComputeSex Job Started: $(date "+%Y-%m-%d %H:%M:%S")"
echo "============================================"
echo "PRJdir: $PRJdir"
echo "baseDir: $baseDir"
echo "Logdir: $Logdir"
echo "AccList: $AccList"
echo "scratch_dir: $scratch_dir"

# Array to store samples with conflicts
conflicted_samples=()

# Verify required inputs
if [ -z "$BIOProjectID" ]; then
    echo "ERROR: BIOProjectID is not set."
    exit 1
fi

if [ -z "$baseDir" ]; then
    echo "ERROR: baseDir is not set."
    exit 1
fi

if [ -z "$Logdir" ]; then
    echo "ERROR: Logdir is not set."
    exit 1
fi

if [ -z "$AccList" ] || [ ! -f "$AccList" ]; then
    echo "ERROR: AccList is not set or file does not exist: $AccList"
    exit 1
fi

# Verify that the project directory exists
if [ ! -d "$PRJdir" ]; then
    echo "ERROR: Directory $PRJdir does not exist. Please check the project name."
    exit 1
fi

# Ensure output directories exist
mkdir -p "$baseDir"
mkdir -p "$Logdir/ComputeSex"
mkdir -p "$scratch_dir"

# Navigate to the scratch directory
cd "${scratch_dir}" || {
    echo "ERROR: Failed to change directory to $scratch_dir. Please check the project name.";
    exit 1;
}
echo "Working in directory: ${scratch_dir}"
echo ""

# Redirect output to the log file in the permanent log directory
log_file="${Logdir}/ComputeSex/${BIOProjectID}_computeSex.txt"
exec 1>"$log_file" 2>&1

# Capture start time
start_time=$(date +%s)

# Temporary file to store the results
temp_results="${scratch_dir}/${BIOProjectID}_temp_results.txt"

echo "============================================"
echo "Processing Samples for Sex Computation"
echo "============================================"

process_sample() {
    local run="$1"
    local sample="$2"
    local input_sex="$3"

    # BAM files are in ${scratch_dir}/${sample}/
    local BAM_FILE="${scratch_dir}/${sample}/${sample}_GENOME_SORT.bam"

    echo "Processing sample: $sample (Run: $run)"
    echo "  Expected BAM: $BAM_FILE"

    # Verify sample directory exists
    if [ ! -d "${scratch_dir}/${sample}" ]; then
        echo "  ERROR: Directory ${scratch_dir}/${sample} does not exist."
        return 1
    fi

    # Verify BAM file exists
    if [ ! -f "$BAM_FILE" ]; then
        echo "  ERROR: BAM file not found: $BAM_FILE"
        echo "  Contents of ${scratch_dir}/${sample}/:"
        ls -lh "${scratch_dir}/${sample}/" || echo "  Directory listing failed"
        return 1
    fi

    echo "  BAM file found: $(du -h "$BAM_FILE" | cut -f1)"

    # Extract and calculate read depth of X and Y chromosomes
    x_map=$(samtools idxstats "$BAM_FILE" | grep -P "chrX\s" | awk '{print $3}')
    x_len=$(samtools idxstats "$BAM_FILE" | grep -P "chrX\s" | awk '{print $2}')
    x_cov=$(echo "scale=6; $x_map/$x_len" | bc)

    y_map=$(samtools idxstats "$BAM_FILE" | grep -P "chrY\s" | awk '{print $3}')
    y_len=$(samtools idxstats "$BAM_FILE" | grep -P "chrY\s" | awk '{print $2}')
    y_cov=$(echo "scale=6; $y_map/$y_len" | bc)

    # Check if y_cov is zero to prevent divide by zero error
    if [ "$y_cov" == "0" ] || [ "$y_cov" == "0.000000" ]; then
        ratio="Inf"
        calculated_sex="F"
    else
        ratio=$(echo "scale=6; $x_cov/$y_cov" | bc -l)
        if (( $(echo "$ratio > 40.00" | bc -l) )); then
            calculated_sex="F"
        else
            calculated_sex="M"
        fi
    fi

    # Determine agreement
    if [ "$input_sex" == "$calculated_sex" ]; then
        agreement="Agree"
    else
        agreement="Conflict"
        conflicted_samples+=("$sample")
    fi

    echo "  X coverage: $x_cov, Y coverage: $y_cov, Ratio: $ratio"
    echo "  Input sex: $input_sex, Calculated sex: $calculated_sex, Agreement: $agreement"

    # Output results with five columns
    printf "%s\t%s\t%s\t%s\t%s\n" "$sample" "$input_sex" "$calculated_sex" "$ratio" "$agreement" >> "$temp_results"
}

# Create or clear the temporary results file with header
echo -e "SampleID\tInputSex\tComputedSex\tRatio\tAgreement" > "$temp_results"

# Process each sample from the input file
echo "Starting processing of samples..."
processed=0
failed=0

while IFS=$'\t' read -r Run geo_accession Tissue Strain Sex PMID GEOpath Title Sample_characteristics StrainInfo; do
    # Skip the header row
    if [ "$Run" = "Run" ]; then
        continue
    fi

    # Define sample ID and input sex
    sample="$geo_accession"
    input_sex="$Sex"

    # Process sample
    if process_sample "$Run" "$sample" "$input_sex"; then
        processed=$((processed + 1))
    else
        echo "ERROR: Processing failed for sample $sample"
        failed=$((failed + 1))
    fi
done < "$AccList"

echo ""
echo "============================================"
echo "Processing Summary"
echo "============================================"
echo "Successfully processed: $processed samples"
echo "Failed: $failed samples"

if [ ${#conflicted_samples[@]} -gt 0 ]; then
    echo ""
    echo "WARNING: Sex conflicts found for ${#conflicted_samples[@]} sample(s):"
    printf '  - %s\n' "${conflicted_samples[@]}"
fi

# Capture end time and calculate wall clock time
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
echo ""
echo "Processing completed in $elapsed_time seconds"
echo "============================================"

# Move temporary results to final permanent output file
final_output="${baseDir}/${BIOProjectID}_sex_result.txt"
mv "$temp_results" "$final_output"

echo ""
echo "Output files:"
echo "  - Sex results: $final_output"
echo "  - Log file: $log_file"
echo ""
echo "ComputeSex completed at: $(date "+%Y-%m-%d %H:%M:%S")"
echo "============================================"