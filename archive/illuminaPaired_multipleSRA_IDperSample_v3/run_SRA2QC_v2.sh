#!/bin/bash
#SBATCH --job-name=SRAtoQC
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=12:00:00
#SBATCH --account=akwitek
#SBATCH --partition=normal
#SBATCH --output=%x-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wdemos@mcw.edu
#SBATCH --error=SRAtoQCerr
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
baseDir="$myDir/data/expression/GEO/${BIOProjectID}"
PRJdir="$baseDir/reads_fastq"
Logdir="$baseDir/log_files"
screenconfig="$myDir/dependencies/FastQ_Screen_Genomes/fastq_screen.conf"
scratchDir="/scratch/g/akwitek/wdemos/$BIOProjectID"

# Create directories
mkdir -p "$baseDir"
mkdir -p "$PRJdir"
mkdir -p "$Logdir"

# Redirect output to a master log file
log_file="${Logdir}/${BIOProjectID}_Step1.out"
exec 1>"$log_file" 2>&1

# Read accessions and submit jobs
job_ids=()  # Array to store job IDs
while IFS=$'\t' read -r Run geo_accession Tissue Strain Sex PMID GEOpath Title Sample_characteristics StrainInfo; do
    # Skip the header row
    if [ "$Run" == "Run" ]; then
        continue
    fi

    # Create a unique subdirectory for each geo_accession under reads_fastq
    geo_log_dir="${PRJdir}/${geo_accession}/log_files/SRA2QC"
    mkdir -p "$geo_log_dir"  # Ensure the directory exists

    # Submit each sample as a separate job and capture the job ID
    job_id=$(sbatch --export=BIOProjectID="$BIOProjectID",Run="$Run",geo_accession="$geo_accession",PRJdir="$PRJdir",screenconfig="$screenconfig" \
        --output="${geo_log_dir}/SRA2QC-%j.out" \
        --error="${geo_log_dir}/SRA2QC-%j.err" \
        SRA2QC.sh | awk '{print $4}')

    # Store the job ID in the array
    job_ids+=("$job_id")
done < "$AccList"

echo "All jobs submitted."

# Wait for all jobs to complete
for job_id in "${job_ids[@]}"; do
    while squeue -u $USER | grep -q "$job_id"; do
        sleep 5
    done
done

# Load MultiQC and run it on the PRJdir
module load python/3.9.1 multiqc/1.18
#python3 -m multiqc "$scratchDir" -o "${Logdir}" || { echo "MultiQC failed"; exit 1; }
python3 -m multiqc "$scratchDir" -o "${baseDir}" -n ${BIOProjectID}_fastq_multiQC_report || { echo "MultiQC failed"; exit 1; }
# Optionally move MultiQC results to PRJdir (if needed)
# mv "$scratchDir/multiqc_report.html" "$PRJdir/"
module unload python/3.9.1 multiqc/1.18
echo "MultiQC completed successfully."

