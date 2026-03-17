#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1mb
#SBATCH --time=00:10:00
#SBATCH --account=akwitek
#SBATCH --partition=normal
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wdemos@mcw.edu

# 18 June 2025: Removed --error and --output to avoid SLURM log clutter.
# Redirect master script output
set -euxo pipefail

# Check command line arguments
if [ $# -ne 2 ]; then
    echo "Usage: $0 <SRA_accession_list.txt> <BioProject_ID>"
    exit 1
fi

AccList="$1"
BIOProjectID="$2"
myDir="/home/wdemos"
PRJdir="$myDir/data/expression/GEO/$BIOProjectID/reads_fastq"
Logdir="$myDir/data/expression/GEO/$BIOProjectID"
scratch_dir="/scratch/g/akwitek/wdemos/$BIOProjectID"

# Create needed directories
mkdir -p "$PRJdir" "$Logdir"

# Central master log file
log_file="${Logdir}/${BIOProjectID}_jsons.out"
exec 1>"$log_file" 2>&1

echo "[$(date)] Master script started"
echo "BIOProject: $BIOProjectID"
echo "Accession list: $AccList"

# Job tracking
declare -A job_ids

while IFS=$'\t' read -r run geo_accession Tissue Strain Sex PMID GEOpath Title Sample_characteristics StrainInfo; do
    [[ "$run" == "Run" ]] && continue

    unique_name="${Tissue}_${Strain}_${Sex}_${geo_accession}"
    sample_output_dir="$PRJdir/$geo_accession"
    mkdir -p "$sample_output_dir"

    echo "[$(date)] Submitting jobs for $geo_accession (Run: $run)"

    # Submit bigWig JSON job
    job_ids["$geo_accession,2bigWigjson"]=$(sbatch --job-name="bgJ_$geo_accession" \
        --export=ALL,BIOProjectID="$BIOProjectID",Run="$run",geo_accession="$geo_accession",PRJdir="$PRJdir",tissue="$Tissue",strain="$Strain",sex="$Sex",title="$Title",Sample_characteristics="$Sample_characteristics",PMID="$PMID",GEOpath="$GEOpath",unique_name="$unique_name",Logdir="$Logdir",scratch_dir="$scratch_dir",StrainInfo="$StrainInfo" \
        BWjsonv2.sh | awk '{print $4}')

    echo "Submitted BW JSON job for $geo_accession (Job ID: ${job_ids["$geo_accession,2bigWigjson"]})"

    # Submit gene JSON job
    job_ids["$geo_accession,geneJson"]=$(sbatch --job-name="gj_$geo_accession" \
        --export=ALL,BIOProjectID="$BIOProjectID",Run="$run",geo_accession="$geo_accession",PRJdir="$PRJdir",tissue="$Tissue",strain="$Strain",sex="$Sex",title="$Title",Sample_characteristics="$Sample_characteristics",PMID="$PMID",GEOpath="$GEOpath",unique_name="$unique_name",Logdir="$Logdir",scratch_dir="$scratch_dir",StrainInfo="$StrainInfo" \
        gJSONv2.sh | awk '{print $4}')

    echo "Submitted Gene JSON job for $geo_accession (Job ID: ${job_ids["$geo_accession,geneJson"]})"

    # Submit transcript JSON job
    job_ids["$geo_accession,txBedJson"]=$(sbatch --job-name="tj_$geo_accession" \
        --export=ALL,BIOProjectID="$BIOProjectID",Run="$run",geo_accession="$geo_accession",PRJdir="$PRJdir",tissue="$Tissue",strain="$Strain",sex="$Sex",title="$Title",Sample_characteristics="$Sample_characteristics",PMID="$PMID",GEOpath="$GEOpath",unique_name="$unique_name",Logdir="$Logdir",scratch_dir="$scratch_dir",StrainInfo="$StrainInfo" \
        tJSONv2.sh | awk '{print $4}')

    echo "Submitted Transcript JSON job for $geo_accession (Job ID: ${job_ids["$geo_accession,txBedJson"]})"

done < "$AccList"

echo "[$(date)] All jobs submitted."
