#!/bin/bash
#SBATCH --job-name=computeSex
#SBATCH --ntasks=1
#SBATCH --mem=1gb
#SBATCH --time=00:10:00
#SBATCH --account=akwitek
#SBATCH --partition=normal
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wdemos@mcw.edu

# Script details
## Written by Wendy Demos
## Last Modified 29 Oct 2024 by WMD to run on RCC cluster
## 18 June 2025 removed from header: #SBATCH --output=%x-%j.out, #SBATCH --error=computeSex.err


# Enable debugging by printing each command before execution
set -x

module load samtools/1.20 

# Retrieve environment variables
Run=${Run}
geo_accession=${geo_accession}
PRJdir=${PRJdir}
BIOProjectID=${BIOProjectID}
echo "PRJdir is set to: $PRJdir"
scratch_dir="/scratch/g/akwitek/wdemos/${BIOProjectID}"
FinalOPdir=${Logdir}
echo "Logdir is set to: $Logdir"
 # Ensure this path is correct
conflicted_samples=()  # Array to store samples with conflicts

# Verify that the project directory exists
if [ ! -d "$PRJdir" ]; then
    echo "Error: Directory $PRJdir does not exist. Please check the project name."
    exit 1
fi

# Navigate to the scratch directory
cd "${scratch_dir}" || { echo "Failed to change directory to $scratch_dir. Please check the project name."; exit 1; }
echo "Working in directory: ${scratch_dir}"

# Redirect output to the log file in the project directory
log_file="${scratch_dir}/${BIOProjectID}_computeSex.txt"
exec 1>"$log_file" 2>&1

# Capture start time
start_time=$(date +%s)

# Temporary file to store the results
temp_results="${scratch_dir}/${BIOProjectID}_temp_results.txt"

process_sample() {
    local run="$1"
    local sample="$2"
    local input_sex="$3"
    local BAM_FILE="${scratch_dir}/${run}/${sample}_GENOME_SORT.bam"

    if [ ! -d "${scratch_dir}/${run}" ]; then
        echo "Error: Directory ${scratch_dir}/${run} does not exist."
        return 1
    fi
    if [ ! -f "$BAM_FILE" ]; then
        echo "Error: BAM file $BAM_FILE not found for sample $sample."
        return 1
    fi

    # Extract and calculate read depth of X and Y chromosomes
    x_map=$(samtools idxstats "$BAM_FILE" | grep -P "chrX\s" | awk '{print $3}')
    x_len=$(samtools idxstats "$BAM_FILE" | grep -P "chrX\s" | awk '{print $2}')
    x_cov=$(echo "scale=3; $x_map/$x_len" | bc)

    y_map=$(samtools idxstats "$BAM_FILE" | grep -P "chrY\s" | awk '{print $3}')
    y_len=$(samtools idxstats "$BAM_FILE" | grep -P "chrY\s" | awk '{print $2}')
    y_cov=$(echo "scale=6; $y_map/$y_len" | bc)

    # Check if y_cov is zero to prevent divide by zero error
    if [ "$y_cov" == "0" ]; then
        ratio=$(echo "scale=6; $x_cov" | bc)
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
        conflicted_samples+=("$sample")  # Add to the array of conflicted samples
    fi

    # Output results with four columns
#    printf "%s\t%s\t%s\t%s\n" "$sample" "$input_sex" "$calculated_sex" "$agreement" >> "$temp_results"
    # Output results with five columns
    printf "%s\t%s\t%s\t%s\t%s\n" "$sample" "$input_sex" "$calculated_sex" "$ratio" "$agreement" >> "$temp_results"
}

# Create or clear the temporary results file
#echo -e "SampleID\tInputSex\tComputedSex\tAgreement" > "$temp_results"
# 5 column output
echo -e "SampleID\tInputSex\tComputedSex\tRatio\tAgreement" > "$temp_results"

# Process each sample from the input file
echo "Starting processing of samples."
while IFS=$'\t' read -r Run geo_accession Tissue Strain Sex PMID GEOpath Title Sample_characteristics StrainInfo; do
    # Skip the header row
    if [ "$Run" = "Run" ]; then
        continue
    fi

    # Define sample ID and input sex
    sample="$geo_accession"
    input_sex="$Sex"

    # Process sample
    process_sample "$Run" "$sample" "$input_sex"
    if [ $? -ne 0 ]; then
        echo "Processing failed for sample $sample"
    fi
done < "$AccList"
echo "Finished processing samples."

# Capture end time and calculate wall clock time
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
echo "Processing completed in $elapsed_time seconds."

# Rename the temporary results file to final output file
mv "$temp_results" "${scratch_dir}/${BIOProjectID}_sex_result.txt"
cp "${scratch_dir}/${BIOProjectID}_sex_result.txt" $Logdir
cp "${scratch_dir}/${BIOProjectID}_computeSex.txt" $Logdir
echo "computation of sex is completed, a copy of the output is saved in the $Logdir"
