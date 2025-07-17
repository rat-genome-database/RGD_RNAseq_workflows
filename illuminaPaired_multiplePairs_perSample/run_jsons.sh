#!/bin/sh

# Check command line arguments
if [ $# -ne 2 ]; then
    echo "Usage: $0 <SRA_accession_list.txt> <BioProject_ID>"
    exit 1
fi

AccList=$1  # Accession list file
BIOProjectID=$2  # Project ID for output directory
length=$3
scratch_dir="/scratch/g/akwitek/wdemos/${BIOProjectID}"
log_file="${scratch_dir}/run_combination_log.txt"
INDEX_DIR="${scratch_dir}/RefIndex"
Logdir="${scratch_dir}/logs"
PRJdir="/home/wdemos/data/expression/GEO/${BIOProjectID}/reads_fastq"

# Ensure necessary directories exist
mkdir -p "$scratch_dir" "$Logdir"

# Clear the log file before starting
echo "Log of combined runs for each sample" > "$log_file"

###################################
# Generate JSON files #
###################################

declare -A json_job_ids  # Track job IDs for JSON jobs
declare -A seen_samples  # Track processed geo_accession values

while IFS=$'\t' read -r run geo_accession Tissue Strain Sex PMID GEOpath Title Sample_characteristics StrainInfo; do
    # Skip the header row
    if [ "$run" == "Run" ]; then
        continue
    fi

    # Process each geo_accession only once
    if [[ -n "${seen_samples[$geo_accession]}" ]]; then
        continue  # Skip if this geo_accession was already processed
    fi
    seen_samples["$geo_accession"]=1  # Mark as processed

    # Create a unique filename using Tissue, Sex, and Strain
    unique_name="${Tissue}_${Strain}_${Sex}_${geo_accession}"

    # Define output paths for this sample
    sample_output_dir="$PRJdir/$geo_accession"
    mkdir -p "$sample_output_dir"

    echo "Submitting JSON jobs for $geo_accession" >> "$log_file"

    # Submit jobs **without dependency**
    for job_type in "2bigWigjson" "geneJson" "txBedJson"; do
        script_name=""
        case "$job_type" in
            "2bigWigjson") script_name="BWjson.sh" ;;
            "geneJson") script_name="gJSON.sh" ;;
            "txBedJson") script_name="tJSON.sh" ;;
        esac

        json_job_ids["$geo_accession,$job_type"]=$(sbatch --job-name="${job_type}_$geo_accession" \
            --export=BIOProjectID="$BIOProjectID",Run="$run",geo_accession="$geo_accession",PRJdir="$PRJdir",tissue="$Tissue",strain="$Strain",sex="$Sex",title="$Title",Sample_characteristics="$Sample_characteristics",PMID="$PMID",GEOpath="$GEOpath",unique_name="$unique_name",Logdir="$Logdir",scratch_dir="$scratch_dir",StrainInfo="$StrainInfo" \
            --output="$sample_output_dir/${job_type}-%j.out" --error="$sample_output_dir/${job_type}-%j.err" "$script_name" | awk '{print $4}')

        if [ -z "${json_job_ids["$geo_accession,$job_type"]}" ]; then
            echo "Error submitting $job_type job for $geo_accession" >> "$log_file"
        else
            echo "Submitted $job_type for $geo_accession with job ID: ${json_job_ids["$geo_accession,$job_type"]}" >> "$log_file"
        fi
    done

done < "$AccList"

echo "All JSON generation jobs have been submitted." >> "$log_file"
