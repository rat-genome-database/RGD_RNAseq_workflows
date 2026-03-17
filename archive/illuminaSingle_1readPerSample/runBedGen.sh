#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1gb
#SBATCH --time=01:00:00
#SBATCH --account=akwitek
#SBATCH --partition=normal
#SBATCH --output=%x-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wdemos@mcw.edu
#SBATCH --error=bedNjson.err

# Enable debugging by printing each command before execution
#set -x

# Check command line arguments
if [ $# -ne 2 ]; then
    echo "Usage: $0 <SRA_accession_list.txt> <BioProject_ID>"
    exit 1
fi

AccList=$1 
BIOProjectID=$2
myDir="/home/wdemos"
PRJdir="$myDir/data/expression/GEO/$BIOProjectID/reads_fastq"
Logdir="$myDir/data/expression/GEO/$BIOProjectID"
scratch_dir="/scratch/g/akwitek/wdemos/$BIOProjectID"
INDEX_DIR="$scratch_dir/RefIndex"


# Create directories
mkdir -p "$PRJdir"
mkdir -p "$Logdir"

# Redirect output to log file
log_file="${Logdir}/${BIOProjectID}_bedsAndjsons.out"
exec 1>"$log_file" 2>&1

# Read accessions and prepare jobs for each sample
declare -A job_ids

while IFS=$'\t' read -r run geo_accession Tissue Strain Sex PMID GEOpath Title Sample_characteristics StrainInfo; do
    # Skip the header row
    if [ "$run" == "Run" ]; then
        continue
    fi

    # Create a unique filename using Tissue, Sex, and Strain
    unique_name="${Tissue}_${Strain}_${Sex}_${geo_accession}"

    # Define output paths for this sample
    sample_output_dir="$PRJdir/$geo_accession"
    mkdir -p "$sample_output_dir"

    # Debugging: Print the value of `Run` before submitting job
    echo "Debug: Submitting job for geo_accession $geo_accession, Run: $run" >> "$log_file"


    # Submit gene TPM bed job and capture its job ID
    job_ids["$geo_accession,GeneTPMbed"]=$(sbatch --job-name="gb_$geo_accession" \
        --export=scratch_dir="$scratch_dir",BIOProjectID="$BIOProjectID",Run="$run",geo_accession="$geo_accession",PRJdir="$PRJdir",tissue="$Tissue",strain="$Strain",sex="$Sex",unique_name="$unique_name",Logdir="$Logdir",StrainInfo="$StrainInfo" \
        --output="$sample_output_dir/Gbed-%j.out" --error="$sample_output_dir/Gbed-%j.err" GeneTPMbed.sh | awk '{print $4}')

    if [ $? -ne 0 ]; then
        echo "Error submitting gene TPM bed job for $geo_accession" >> "$log_file"
        continue
    fi

    echo "Gene TPM bed job submitted for $geo_accession with job ID: ${job_ids["$geo_accession,GeneTPMbed"]}" >> "$log_file"

    # Submit transcript TPM bed job and capture its job ID)
    job_ids["$geo_accession,txTPMbed"]=$(sbatch --job-name="txb_$geo_accession" \
        --export=scratch_dir="$scratch_dir",BIOProjectID="$BIOProjectID",Run="$run",geo_accession="$geo_accession",PRJdir="$PRJdir",tissue="$Tissue",strain="$Strain",sex="$Sex",unique_name="$unique_name",Logdir="$Logdir",StrainInfo="$StrainInfo" \
        --output="$sample_output_dir/TxBed-%j.out" --error="$sample_output_dir/TxBed-%j.err" txTPMbed.sh | awk '{print $4}')

    if [ $? -ne 0 ]; then
        echo "Error submitting transcript TPM bed job for $geo_accession" >> "$log_file"
        continue
    fi

done < "$AccList"

echo "All jobs for gene TPM bed JSON file generation are complete." >> "$log_file"

