#!/usr/bin/env bash
#SBATCH --job-name=SRAtoQC_SE
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=12:00:00
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=normal
#SBATCH --output=%x-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=your@email.edu
#SBATCH --error=SRAtoQC_SE.err

###############################################################################
# USER CONFIGURATION
# Update the variables below for your system before running.
###############################################################################
# Your home/base directory
myDir="/path/to/home"
# Directory containing all pipeline scripts
SCRIPT_DIR="$myDir/STAR_RSEM_pipeline/RGD_Illumina_SingleEnd_RNAseq_pipeline"
# Scratch filesystem mount point (for disk checks)
SCRATCH_BASE="/path/to/scratch"
# FastQ-Screen config file
screenconfig="$myDir/dependencies/FastQ_Screen_Genomes/fastq_screen.conf"
###############################################################################

set -x

# Check command line arguments
if [ $# -ne 2 ]; then
    echo "Usage: $0 <SRA_accession_list.txt> <BioProject_ID>"
    exit 1
fi

AccList=$1
BIOProjectID=$2
baseDir="$myDir/data/expression/GEO/${BIOProjectID}"
PRJdir="$baseDir/reads_fastq"
Logdir="$baseDir/log_files"
scratchDir="${SCRATCH_BASE}/${BIOProjectID}"

# Create directories
mkdir -p "$baseDir"
mkdir -p "$PRJdir"
mkdir -p "$Logdir"

# Redirect output to a master log file
log_file="${Logdir}/${BIOProjectID}_Step1.out"
exec 1>"$log_file" 2>&1

echo "===== STEP 1: SRA Download and QC (Single-End) ====="
echo "BioProject: $BIOProjectID"
echo "AccList:    $AccList"
echo "Started:    $(date)"
echo "==========================================="

###############################################################################
# DISK GUARD
###############################################################################
SCRATCH_WARN_PCT="${SCRATCH_WARN_PCT:-85}"
SCRATCH_HOLD_PCT="${SCRATCH_HOLD_PCT:-90}"

scratch_pct_used() {
    local pct
    pct=$(df --output=pcent "$SCRATCH_BASE" 2>/dev/null | tail -1 | tr -d ' %')
    [[ "$pct" =~ ^[0-9]+$ ]] || pct=0
    echo "$pct"
}

check_scratch_or_exit() {
    local context="$1"
    local pct
    pct=$(scratch_pct_used)
    if [[ $pct -ge $SCRATCH_HOLD_PCT ]]; then
        echo "ERROR [$context]: Scratch at ${pct}% (>= ${SCRATCH_HOLD_PCT}% hold threshold). Aborting."
        exit 1
    elif [[ $pct -ge $SCRATCH_WARN_PCT ]]; then
        echo "WARN  [$context]: Scratch at ${pct}% (>= ${SCRATCH_WARN_PCT}% warn threshold). Continuing."
    else
        echo "DISK OK [$context]: Scratch at ${pct}%"
    fi
}

check_scratch_or_exit "pre-flight"

# Read accessions and submit jobs
job_ids=()
while IFS=$'\t' read -r Run geo_accession Tissue Strain Sex PMID GEOpath Title Sample_characteristics StrainInfo; do
    if [ "$Run" == "Run" ]; then
        continue
    fi

    geo_log_dir="${PRJdir}/${geo_accession}/log_files/SRA2QC"
    mkdir -p "$geo_log_dir"

    check_scratch_or_exit "before submitting $geo_accession"

    job_id=$(sbatch \
        --export=BIOProjectID="$BIOProjectID",Run="$Run",geo_accession="$geo_accession",PRJdir="$PRJdir",screenconfig="$screenconfig" \
        --output="${geo_log_dir}/SRA2QC-%j.out" \
        --error="${geo_log_dir}/SRA2QC-%j.err" \
        SRA2QC_SE_v1.sh | awk '{print $4}')

    job_ids+=("$job_id")
    echo "Submitted SRA2QC_SE job for $Run (GSM: $geo_accession) with job ID: $job_id"
done < "$AccList"

echo "All SRA2QC_SE jobs submitted. Total jobs: ${#job_ids[@]}"
echo "==========================================="

# Wait for all jobs to complete
echo "Waiting for all SRA2QC_SE jobs to complete..."
for job_id in "${job_ids[@]}"; do
    while squeue -u $USER | grep -q "$job_id"; do
        sleep 5
    done
done

echo "All SRA2QC_SE jobs have finished. Checking completion status..."

# Check job success
failed_jobs=()
for job_id in "${job_ids[@]}"; do
    sleep 2
    state=$(sacct -j "$job_id" --format=State --noheader 2>/dev/null | head -1 | tr -d ' ')
    if [[ "$state" != "COMPLETED" ]]; then
        echo "ERROR: Job $job_id failed with state: $state"
        failed_jobs+=("$job_id")
    else
        echo "Job $job_id completed successfully"
    fi
done

if [ ${#failed_jobs[@]} -gt 0 ]; then
    echo "==========================================="
    echo "ERROR: ${#failed_jobs[@]} SRA2QC_SE job(s) failed:"
    printf '  - Job ID: %s\n' "${failed_jobs[@]}"
    echo "==========================================="
    echo "Cannot proceed to Step 2. Please check the failed job logs."
    echo "STEP 1 FAILED at $(date)"
    exit 1
fi

echo "All SRA2QC_SE jobs completed successfully"

# Run MultiQC
echo "==========================================="
echo "Running MultiQC to generate QC report..."
module load python/3.9.1 multiqc/1.18
python3 -m multiqc "$scratchDir" -o "${baseDir}" -n ${BIOProjectID}_fastq_multiQC_report || {
    echo "MultiQC failed"
    module unload python/3.9.1 multiqc/1.18
    exit 1
}
module unload python/3.9.1 multiqc/1.18
echo "MultiQC completed successfully."
echo "Report saved to: ${baseDir}/${BIOProjectID}_fastq_multiQC_report.html"

# Create success marker
SUCCESS_DIR="$baseDir/.status"
mkdir -p "$SUCCESS_DIR"
touch "${SUCCESS_DIR}/.step1_complete"

echo "==========================================="
echo "STEP 1 COMPLETED SUCCESSFULLY"
echo "Completed: $(date)"
echo "==========================================="
echo ""
echo "Summary:"
echo "  - Total SRA runs processed: ${#job_ids[@]}"
echo "  - Failed jobs: ${#failed_jobs[@]}"
echo "  - FASTQ files location: $scratchDir/<run>/"
echo "  - MultiQC report: ${baseDir}/${BIOProjectID}_fastq_multiQC_report.html"
echo "  - Success marker created: ${SUCCESS_DIR}/.step1_complete"
echo ""
echo "You may now proceed to Step 2 (run_RNApipeline_SE_diskGuard.sh)"

exit 0
