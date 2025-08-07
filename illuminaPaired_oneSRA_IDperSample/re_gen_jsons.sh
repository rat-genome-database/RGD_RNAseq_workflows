#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1mb
#SBATCH --time=00:10:00
#SBATCH --account=akwitek
#SBATCH --partition=normal
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wdemos@mcw.edu
#SBATCH --output=%x-%j.out

# Enable debugging by printing each command before execution
set -x

# Check command line arguments
if [ $# -ne 2 ]; then
    echo "Usage: $0 <SRA_accession_list.txt> <BioProject_ID>"
    exit 1
fi

AccList=$1 
BIOProjectID=$2
myDir="/home/wdemos"
PRJdir="$myDir/data/expression/GEO/$BIOProjectID/reads_fastq"
Logdir="$myDir/data/expression/GEO/$BIOProjectID/logfiles"
scratch_dir="/scratch/g/akwitek/wdemos/$BIOProjectID"
INDEX_DIR="$scratch_dir/RefIndex"
baseDir="$myDir/data/expression/GEO/$BIOProjectID"

# Create output and logging directories
mkdir -p "$PRJdir"
mkdir -p "$Logdir"

# Log file for this submission script (not individual JSON jobs)
log_file="${Logdir}/${BIOProjectID}_redo_jsons.out"
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
    sample_log_dir="$sample_output_dir/log_files"
    mkdir -p "$sample_log_dir"

    echo "Submitting JSON jobs for $geo_accession" >> "$log_file"

    # Submit bigWig JSON job
    job_ids["$geo_accession,2bigWigjson"]=$(sbatch --job-name="bgJ_$geo_accession" \
        --export=baseDir="$baseDir",BIOProjectID="$BIOProjectID",Run="$run",geo_accession="$geo_accession",PRJdir="$PRJdir",tissue="$Tissue",strain="$Strain",sex="$Sex",title="$Title",Sample_characteristics="$Sample_characteristics",PMID="$PMID",GEOpath="$GEOpath",unique_name="$unique_name",Logdir="$Logdir",scratch_dir="$scratch_dir",StrainInfo="$StrainInfo" \
        --output="$sample_log_dir/bwjson-%j.out" --error="$sample_log_dir/bwjson-%j.err" BWjson_v4.sh | awk '{print $4}')

    if [ $? -ne 0 ]; then
        echo "Error submitting bigWig json job for $geo_accession" >> "$log_file"
        continue
    fi

    echo -e "$geo_accession\t${job_ids["$geo_accession,2bigWigjson"]}\tbwjson" >> "$Logdir/${BIOProjectID}_jobids.tsv"

    # Submit gene TPM bed JSON job
    job_ids["$geo_accession,geneJson"]=$(sbatch --job-name="gj_$geo_accession" \
        --export=baseDir="$baseDir",BIOProjectID="$BIOProjectID",Run="$run",geo_accession="$geo_accession",PRJdir="$PRJdir",tissue="$Tissue",strain="$Strain",sex="$Sex",title="$Title",Sample_characteristics="$Sample_characteristics",PMID="$PMID",GEOpath="$GEOpath",unique_name="$unique_name",Logdir="$Logdir",scratch_dir="$scratch_dir",StrainInfo="$StrainInfo" \
        --output="$sample_log_dir/GbedJSON-%j.out" --error="$sample_log_dir/GbedJSON-%j.err" gJSON_v4.sh | awk '{print $4}')

    if [ $? -ne 0 ]; then
        echo "Error submitting gene TPM bed JSON job for $geo_accession" >> "$log_file"
        continue
    fi

    echo -e "$geo_accession\t${job_ids["$geo_accession,geneJson"]}\tgeneJson" >> "$Logdir/${BIOProjectID}_jobids.tsv"

    # Submit transcript bed JSON job
    job_ids["$geo_accession,txBedJson"]=$(sbatch --job-name="tj_$geo_accession" \
        --export=baseDir="$baseDir",scratch_dir="$scratch_dir",BIOProjectID="$BIOProjectID",Run="$run",geo_accession="$geo_accession",PRJdir="$PRJdir",tissue="$Tissue",strain="$Strain",sex="$Sex",title="$Title",Sample_characteristics="$Sample_characteristics",PMID="$PMID",GEOpath="$GEOpath",unique_name="$unique_name",Logdir="$Logdir",StrainInfo="$StrainInfo" \
        --output="$sample_log_dir/TxBedJSON-%j.out" --error="$sample_log_dir/TxBedJSON-%j.err" tJSON_v4.sh | awk '{print $4}')

    if [ $? -ne 0 ]; then
        echo "Error submitting transcript bed JSON job for $geo_accession" >> "$log_file"
        continue
    fi

    echo -e "$geo_accession\t${job_ids["$geo_accession,txBedJson"]}\ttxBedJson" >> "$Logdir/${BIOProjectID}_jobids.tsv"

done < "$AccList"

echo "All jobs for gene and transcript JSON generation have been submitted." >> "$log_file"
