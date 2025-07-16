#!/bin/bash
#SBATCH --job-name=matrix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=3gb
#SBATCH --time=06:00:00
#SBATCH --account=akwitek
#SBATCH --output=%x-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wdemos@mcw.edu
#SBATCH --error=matrix.err

# Enable debugging by printing each command before execution
set -x

# Check command line arguments
if [ $# -ne 3 ]; then
    echo "Usage: $0 <SRA_accession_list.txt> <BioProject_ID> <Read_Length>"
    exit 1
fi

AccList=$1 
BIOProjectID=$2
length=$3
myDir="/home/wdemos"
PRJdir="$myDir/data/expression/GEO/$BIOProjectID/reads_fastq"
Logdir="$myDir/data/expression/GEO/$BIOProjectID"
scratch_dir="/scratch/g/akwitek/wdemos/$BIOProjectID"
INDEX_DIR="$scratch_dir/RefIndex"
REFdir="/home/wdemos/NCBI_1August2023"


# Create directories
mkdir -p "$PRJdir"
mkdir -p "$Logdir"

# Redirect output to log file
log_file="${Logdir}/${BIOProjectID}_matrix.out"
exec 1>"$log_file" 2>&1

# Function to wait for jobs to complete
wait_for_jobs() {
    for job_id in "${job_ids[@]}"; do
        while squeue -u $USER | grep -q "$job_id"; do
            sleep 5
        done
    done
}

###########################
# Final Steps #
###########################
echo "Generating combined matrices for all samples and final MultiQC report."
job_id=$(sbatch --export=Run="$run",geo_accession="$geo_accession",PRJdir="$PRJdir",Logdir="$Logdir",AccList="$AccList",BIOProjectID="$BIOProjectID" \
              --output="$Logdir/RSEMmatrix-%j.out" \
              --error="$Logdir/RSEMmatrix-%j.err" RSEMmatrix.sh)
if [ $? -ne 0 ]; then
    echo "Error running RSEMmatrix.sh"
    exit 1
fi
# Preserve bam files only upon successful processing updated 10 April 2025
#echo "${Logdir}/reads_fastq/${geo_accession}/"
#cp "${scratch_dir}/${Run}/${geo_accession}_GENOME_SORT.bam" ${Logdir}/reads_fastq/${geo_accession}/
#cp "${scratch_dir}/${Run}/${geo_accession}_GENOME_SORT.bam.bai" ${Logdir}/reads_fastq/${geo_accession}/
#cp "${scratch_dir}/${Run}/${geo_accession}.transcript.sorted.bam" ${Logdir}/reads_fastq/${geo_accession}/
#cp "${scratch_dir}/${Run}/${geo_accession}.transcript.sorted.bam.bai" ${Logdir}/reads_fastq/${geo_accession}/

echo "RGD STAR-RSEM workflow is complete."
