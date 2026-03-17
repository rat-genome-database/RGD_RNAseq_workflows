#!/bin/bash
#SBATCH --job-name=computeSex
#SBATCH --ntasks=1
#SBATCH --mem=1gb
#SBATCH --time=00:10:00
#SBATCH --account=akwitek
#SBATCH --partition=normal
#SBATCH --output=%x-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wdemos@mcw.edu
#SBATCH --error=computeSex.err

# Script details
## Written by Wendy Demos
## Last Modified 28 Feb 2025 by WMD to run on RCC cluster
# Enable debugging by printing each command before execution
set -x

module load samtools/1.20 

# Retrieve environment variables
BIOProjectID=$1  # Pass only BIOProjectID, not geo_accession
PRJdir=${PRJdir}
WRKdir="/home/wdemos/data/expression/GEO/${BIOProjectID}"
Logdir=${Logdir}
AccList=${AccList}
scratch_dir="/scratch/g/akwitek/wdemos/${BIOProjectID}"
FinalOPdir=${Logdir}

# Ensure directories exist
if [ ! -d "$PRJdir" ] || [ ! -d "$scratch_dir" ]; then
    echo "Error: Directory $PRJdir or $scratch_dir does not exist. Please check the project name."
    exit 1
fi

# Navigate to the scratch directory
cd "${scratch_dir}" || { echo "Failed to change directory to $scratch_dir"; exit 1; }

# Redirect output to log file
log_file="${scratch_dir}/${BIOProjectID}_computeSex.log"
exec 1>"$log_file" 2>&1

# Capture start time
start_time=$(date +%s)

# Output file
#output_file="${scratch_dir}/${BIOProjectID}_sex_result.txt"
output_file="${WRKdir}/${BIOProjectID}_sex_result.txt"
echo -e "SampleID\tInputSex\tComputedSex\tRatio\tAgreement" > "$output_file"

# Function to process each unique sample
process_sample() {
    local sample="$1"
    local input_sex="$2"
    local BAM_FILE="${scratch_dir}/${sample}/${sample}_GENOME_SORT.bam"

    # Ensure BAM file exists
    if [ ! -f "$BAM_FILE" ]; then
        echo -e "$sample\t$input_sex\tERROR\tBAM not found\tNA" >> "$output_file"
        return
    fi

    # Compute read depth for X and Y chromosomes
    x_map=$(samtools idxstats "$BAM_FILE" | awk '$1 == "chrX" {print $3}')
    x_len=$(samtools idxstats "$BAM_FILE" | awk '$1 == "chrX" {print $2}')
    y_map=$(samtools idxstats "$BAM_FILE" | awk '$1 == "chrY" {print $3}')
    y_len=$(samtools idxstats "$BAM_FILE" | awk '$1 == "chrY" {print $2}')

    # Compute coverage ratios
    x_cov=$(echo "scale=6; $x_map/$x_len" | bc)
    y_cov=$(echo "scale=6; $y_map/$y_len" | bc)
    
    # Prevent division by zero
    if [ "$y_cov" == "0" ] || [ -z "$y_cov" ]; then
        ratio="$x_cov"
        computed_sex="F"
    else
        ratio=$(echo "scale=6; $x_cov/$y_cov" | bc)
        if (( $(echo "$ratio > 40.00" | bc -l) )); then
            computed_sex="F"
        else
            computed_sex="M"
        fi
    fi

    # Determine agreement based on the correct InputSex
    if [ "$input_sex" == "$computed_sex" ]; then
        agreement="Agree"
    else
        agreement="Conflict"
    fi

    # Write results
    echo -e "$sample\t$input_sex\t$computed_sex\t$ratio\t$agreement" >> "$output_file"
}

# Process each unique geo_accession, extracting the correct InputSex (column 5)
echo "Processing samples..."
awk -F'\t' 'NR>1 && !seen[$2]++ {print $2, $5}' "$AccList" | while read -r geo_accession input_sex; do
    process_sample "$geo_accession" "$input_sex"
done

# Capture end time and duration
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
echo "Processing completed in $elapsed_time seconds."

# Move results to log directory
#cp "$output_file" "$Logdir"
#echo "Results saved to $Logdir/${BIOProjectID}_sex_result.txt"
echo "Results saved to $WRKdir/${BIOProjectID}_sex_result.txt"
