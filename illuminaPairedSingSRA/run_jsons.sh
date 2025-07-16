#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1mb
#SBATCH --time=00:10:00
#SBATCH --account=akwitek
#SBATCH --partition=normal
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wdemos@mcw.edu

# 18 June 2025 removed #SBATCH --error=json.err and #SBATCH --output=%x-%j.out from header to reduce redundant log files
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
log_file="${Logdir}/${BIOProjectID}_jsons.out"
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

    # Submit bigWig JSON job
    job_ids["$geo_accession,2bigWigjson"]=$(sbatch --job-name="bgJ_$geo_accession" \
        --export=BIOProjectID="$BIOProjectID",Run="$run",geo_accession="$geo_accession",PRJdir="$PRJdir",tissue="$Tissue",strain="$Strain",sex="$Sex",title="$Title",Sample_characteristics="$Sample_characteristics",PMID="$PMID",GEOpath="$GEOpath",unique_name="$unique_name",Logdir="$Logdir",scratch_dir="$scratch_dir",StrainInfo="$StrainInfo" \
        --output="$sample_output_dir/bwjson-%j.log" --error="$sample_output_dir/bwjson-%j.log" BWjson.sh | awk '{print $4}')

    if [ $? -ne 0 ]; then
        echo "Error submitting bigWig json job for $geo_accession" >> "$log_file"
        continue
    fi

    echo "BigWig json generation job submitted for $geo_accession with job ID: ${job_ids["$geo_accession,2bigWigjson"]}" >> "$log_file"

 
    # Submit gene TPM bed JSON job and capture its job ID (dependent on transcript TPM bed job)
    job_ids["$geo_accession,geneJson"]=$(sbatch --job-name="gj_$geo_accession" \
        --export=BIOProjectID="$BIOProjectID",Run="$run",geo_accession="$geo_accession",PRJdir="$PRJdir",tissue="$Tissue",strain="$Strain",sex="$Sex",title="$Title",Sample_characteristics="$Sample_characteristics",PMID="$PMID",GEOpath="$GEOpath",unique_name="$unique_name",Logdir="$Logdir",scratch_dir="$scratch_dir",StrainInfo="$StrainInfo" \
        --output="$sample_output_dir/GbedJSON-%j.log" --error="$sample_output_dir/GbedJSON-%j.log" gJSON.sh | awk '{print $4}')

    if [ $? -ne 0 ]; then
        echo "Error submitting gene TPM bed JSON job for $geo_accession" >> "$log_file"
        continue
    fi

    echo "Gene TPM bed JSON job submitted for $geo_accession with job ID: ${job_ids["$geo_accession,geneJson"]}" >> "$log_file"

    # Submit transcript bed JSON job and capture its job ID (dependent on gene JSON job)
    job_ids["$geo_accession,txBedJson"]=$(sbatch --job-name="tj_$geo_accession" \
        --export=scratch_dir="$scratch_dir",BIOProjectID="$BIOProjectID",Run="$run",geo_accession="$geo_accession",PRJdir="$PRJdir",tissue="$Tissue",strain="$Strain",sex="$Sex",title="$Title",Sample_characteristics="$Sample_characteristics",PMID="$PMID",GEOpath="$GEOpath",unique_name="$unique_name",Logdir="$Logdir",StrainInfo="$StrainInfo" \
        --output="$sample_output_dir/TxBedJSON-%j.log" --error="$sample_output_dir/TxBedJSON-%j.log" tJSON.sh | awk '{print $4}')

    if [ $? -ne 0 ]; then
        echo "Error submitting transcript bed JSON job for $geo_accession" >> "$log_file"
        continue
    fi

    echo "Transcript bed JSON job submitted for $geo_accession with job ID: ${job_ids["$geo_accession,txBedJson"]}" >> "$log_file"

done < "$AccList"

echo "All jobs for gene TPM bed JSON file generation are complete." >> "$log_file"

