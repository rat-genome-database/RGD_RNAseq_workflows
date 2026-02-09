#!/bin/bash
#SBATCH --job-name=S2R
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=06:30:00
#SBATCH --account=akwitek
#SBATCH --output=%x-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wdemos@mcw.edu
#SBATCH --error=S2J.err

# Enable strict error handling and debugging
set -euo pipefail
set -x

LIB="/home/wdemos/STAR_RSEM_pipeline/bulkDir_illuminaPaired_oneSRA_IDperSample/lib_v10.sh"

if [[ ! -f "$LIB" ]]; then
    echo "FATAL: lib_v10.sh not found at $LIB" >&2
    exit 1
fi

source "$LIB"

# Configuration: Orchestrated mode disables duplicate notifications
ORCHESTRATED_MODE="${ORCHESTRATED_MODE:-false}"
KEEP_SCRATCH="${KEEP_SCRATCH:-false}"

# Check command line arguments
if [ $# -ne 3 ]; then
    echo "Usage: $0 <SRA_accession_list.txt> <BioProject_ID> <Read_Length>"
    exit 1
fi

AccList=$1
BIOProjectID=$2
length=$3

# Validate inputs
if [[ ! -f "$AccList" ]]; then
    echo "ERROR: Accession list file not found: $AccList"
    exit 1
fi

# Validate BioProject ID for path injection
if [[ "$BIOProjectID" =~ [^a-zA-Z0-9_-] ]]; then
    echo "ERROR: Invalid characters in BioProject ID. Only alphanumeric, underscore, and hyphen allowed."
    exit 1
fi

if [[ ! "$length" =~ ^[0-9]+$ ]]; then
    echo "ERROR: Read length must be a positive integer"
    exit 1
fi

myDir="/home/wdemos"
PRJdir="$myDir/data/expression/GEO/$BIOProjectID/reads_fastq"
Logdir="$myDir/data/expression/GEO/$BIOProjectID/log_files"
scratch_dir="/scratch/g/akwitek/wdemos/$BIOProjectID"
INDEX_DIR="$scratch_dir/RefIndex"
REFdir="/home/wdemos/NCBI_1August2023"
baseDir="$myDir/data/expression/GEO/$BIOProjectID"
SCRATCH_BASE="/scratch/g/akwitek"

# Configuration: Set to "true" to make BWjson failures critical, "false" for warnings only
BWJSON_REQUIRED=true

# Status tracking
SUCCESS_DIR="$baseDir/.status"
mkdir -p "$SUCCESS_DIR"

# Verify STEP1 completed
if [ ! -f "${SUCCESS_DIR}/.step1_complete" ]; then
    echo "ERROR: STEP 1 must complete successfully before STEP 2 can run"
    echo "Missing success marker: ${SUCCESS_DIR}/.step1_complete"
    exit 1
fi

# Create directories
mkdir -p "$PRJdir"
mkdir -p "$Logdir"

# Redirect output to log file
log_file="${Logdir}/${BIOProjectID}_Step2_RNAworkflow.out"
exec 1>"$log_file" 2>&1

echo "===== STEP 2: RNA Processing Pipeline ====="
echo "BioProject: $BIOProjectID"
echo "AccList: $AccList"
echo "Read Length: $length"
echo "Orchestrated Mode: $ORCHESTRATED_MODE"
echo "Started: $(date)"
echo "=========================================="

# Disk space check function
check_disk_space() {
    local mount_point="$1"
    local threshold="$2"

    local usage
    usage=$(df -h "$mount_point" | awk 'NR==2 {print $5}' | sed 's/%//')

    if [[ "$usage" -ge "$threshold" ]]; then
        echo "ERROR: Disk usage at $mount_point is ${usage}% (threshold: ${threshold}%)"
        return 1
    fi

    echo "Disk usage OK: $mount_point at ${usage}%"
    return 0
}

# Check disk space before starting
if ! check_disk_space "$myDir" 95; then
    echo "STEP 2 FAILED for $BIOProjectID - Insufficient disk space in HOME"
    exit 1
fi

if ! check_disk_space "$SCRATCH_BASE" 90; then
    echo "STEP 2 FAILED for $BIOProjectID - Insufficient disk space in SCRATCH"
    exit 1
fi

# Function to wait for jobs to complete
wait_for_jobs() {
    local ids=("$@")

    # If no IDs passed, fall back to global job_ids array
    if [ ${#ids[@]} -eq 0 ]; then
        ids=("${job_ids[@]:-}")
    fi

    for job_id in "${ids[@]:-}"; do
        [ -z "$job_id" ] && continue
        echo "Waiting for job: $job_id"
        # Use only the JOBID column, no header
        while squeue -u "$USER" -h -o "%A" | grep -qx "$job_id"; do
            sleep 5
        done
    done
}

# Function to check if a job completed successfully
check_job_success() {
    local job_id="$1"
    local job_name="$2"

    if [ -z "$job_id" ]; then
        echo "ERROR: No job ID provided for $job_name"
        return 1
    fi

    local max_retries=5
    local retry=0

    while [ $retry -lt $max_retries ]; do
        local state
        state=$(sacct -j "$job_id" --format=State --noheader 2>/dev/null | head -1 | tr -d ' ')

        if [[ -n "$state" ]]; then
            case "$state" in
                COMPLETED)
                    return 0
                    ;;
                FAILED|CANCELLED|TIMEOUT|OUT_OF_MEMORY)
                    echo "ERROR: Job $job_id ($job_name) failed with state: $state"
                    return 1
                    ;;
                *)
                    echo "WARNING: Job $job_id ($job_name) in unexpected state: $state"
                    return 1
                    ;;
            esac
        fi

        sleep 2
        retry=$((retry + 1))
    done

    echo "ERROR: Could not verify completion status for job $job_id ($job_name) after $max_retries attempts"
    return 1
}

# Function to send notifications (skips if in orchestrated mode)
send_notification() {
    local subject="$1"
    local message="$2"

    if [[ "$ORCHESTRATED_MODE" == "true" ]]; then
        echo "NOTIFICATION (suppressed in orchestrated mode): $subject"
        return 0
    fi

    echo -e "$message" | mail -s "$subject" wdemos@mcw.edu 2>/dev/null || true
    echo "NOTIFICATION SENT: $subject"
}

##############################
# STAR Reference Preparation #
##############################
SR_log_dir="$Logdir/STAR"
mkdir -p "$SR_log_dir"
if [ ! -d "$INDEX_DIR" ] || [ ! -f "$INDEX_DIR/SA" ] || [ ! -f "$INDEX_DIR/SAindex" ] || [ ! -f "$INDEX_DIR/Genome" ]; then
    echo "Submitting starRef_v4.sh to generate the STAR reference."
    starRef_job_id=$(sbatch --export=PRJdir="$PRJdir",Logdir="$Logdir",BIOProjectID="$BIOProjectID",length="$length" \
                        --output="$SR_log_dir/starRef-%j.out" \
                        --error="$SR_log_dir/starRef-%j.err" starRef_v4.sh | awk '{print $4}')
    if [ -z "$starRef_job_id" ]; then
        echo "Error: starRef_v4.sh submission failed."
        exit 1
    fi
    echo "STAR reference job ID: $starRef_job_id"

    # Wait and check success
    wait_for_jobs "$starRef_job_id"
    if ! check_job_success "$starRef_job_id" "STAR Reference"; then
        echo "STEP 2 FAILED for $BIOProjectID - STAR reference generation failed"
        exit 1
    fi
else
    echo "STAR reference already exists. Skipping starRef_v4.sh."
    starRef_job_id=""
fi

####################################
# STAR Alignment Jobs for Samples  #
####################################
star_job_ids=()  # Array to store job IDs for STAR alignment
while IFS=$'\t' read -r run geo_accession Tissue Strain Sex PMID GEOpath Title Sample_characteristics StrainInfo || [[ -n "$run" ]]; do
    if [ "$run" == "Run" ]; then continue; fi
    [[ -z "$run" ]] && continue

    unique_name="${Tissue}_${Strain}_${Sex}_${geo_accession}"
    sample_output_dir="$PRJdir/$geo_accession"
    sample_log_dir="$sample_output_dir/log_files/STAR"
    mkdir -p "$sample_log_dir"

    echo "Submitting STAR job for $geo_accession"

    if [ -z "$starRef_job_id" ]; then
        star_job_id=$(sbatch --job-name="${geo_accession}_STAR72" \
                        --export=BIOProjectID="$BIOProjectID",Run="$run",geo_accession="$geo_accession",PRJdir="$PRJdir",length="$length",unique_name="$unique_name",Logdir="$Logdir" \
                        --output="$sample_log_dir/STAR-%j.out" \
                        --error="$sample_log_dir/STAR-%j.err" \
                        star_v5.sh | awk '{print $4}')
    else
        star_job_id=$(sbatch --dependency=afterok:$starRef_job_id \
                        --job-name="${geo_accession}_STAR72" \
                        --export=BIOProjectID="$BIOProjectID",Run="$run",geo_accession="$geo_accession",PRJdir="$PRJdir",length="$length",unique_name="$unique_name",Logdir="$Logdir" \
                        --output="$sample_log_dir/STAR-%j.out" \
                        --error="$sample_log_dir/STAR-%j.err" \
                        star_v5.sh | awk '{print $4}')
    fi

    if [ $? -ne 0 ] || [ -z "$star_job_id" ]; then
        echo "Error submitting STAR job for $geo_accession"
        continue
    fi
    star_job_ids+=("$star_job_id")
done < "$AccList"

# Wait for all STAR alignment jobs
echo "Waiting for STAR alignment jobs to complete..."
wait_for_jobs "${star_job_ids[@]:-}"

# Check that all STAR jobs succeeded
echo "Checking STAR job completion status..."
for star_job_id in "${star_job_ids[@]:-}"; do
    if ! check_job_success "$star_job_id" "STAR alignment"; then
        echo "STEP 2 FAILED for $BIOProjectID - STAR alignment failed"
        exit 1
    fi
done
echo "All STAR alignment jobs completed successfully"

##############################
# STAR Alignment Rate Check  #
##############################
echo "Check sample alignment rates for downstream processing."
STARQC_log_dir="$Logdir/STARQC"
mkdir -p "$STARQC_log_dir"

# Create deduplicated accession list
uniqueAccList="${Logdir}/${BIOProjectID}_Unique_AccList.txt"
{
    head -n 1 "$AccList"
    tail -n +2 "$AccList" | sort -u -t$'\t' -k2,2
} > "$uniqueAccList"
echo "Unique accession list saved to $uniqueAccList"

# Submit STARQC job
starqc_output=$(sbatch \
    --job-name="STARQC" \
    --export=Logdir="$Logdir",AccList="$uniqueAccList",BIOProjectID="$BIOProjectID",PRJdir="$PRJdir" \
    --output="$STARQC_log_dir/STARQC-%j.out" \
    --error="$STARQC_log_dir/STARQC-%j.err" \
    pSTARQC_v1.sh "$uniqueAccList" "$BIOProjectID" 2>&1)

# Extract and validate job ID
starqc_job_id=$(echo "$starqc_output" | grep -oP 'Submitted batch job \K\d+')

if [[ -z "$starqc_job_id" ]]; then
    echo "ERROR: Failed to submit STARQC job"
    echo "Output: $starqc_output"
    exit 1
fi

echo "STARQC job ID: $starqc_job_id"

# Wait for STARQC to complete
echo "Waiting for STARQC job to complete..."
wait_for_jobs "$starqc_job_id"
if ! check_job_success "$starqc_job_id" "STARQC"; then
    echo "STEP 2 FAILED for $BIOProjectID - STARQC failed"
    exit 1
fi
echo "STARQC completed successfully"

#####################################################
# Filter to PASS samples for downstream processing  #
#####################################################
starqc_report="${Logdir}/STARQC/${BIOProjectID}_STAR_Align_sum.txt"
passAccList="${Logdir}/STARQC/${BIOProjectID}_Unique_AccList_PASS.txt"

# Sanity check: STARQC report exists
for i in {1..6}; do
  if [ -s "$starqc_report" ]; then break; fi
  echo "Waiting for STARQC report to appear... ($i/6)"
  sleep 5
done

if [ ! -s "$starqc_report" ]; then
  echo "ERROR: STARQC report not found or empty: $starqc_report"
  echo "Hint: check STARQC logs:"
  echo "  OUT: $STARQC_log_dir/STARQC-$starqc_job_id.out"
  echo "  ERR: $STARQC_log_dir/STARQC-$starqc_job_id.err"
  # Quick peek to help debugging
  tail -n +1 "$STARQC_log_dir/STARQC-$starqc_job_id".{out,err} 2>/dev/null | sed 's/^/STARQC log> /'
  exit 1
fi

# Build PASS-only accession list (keep header; filter by geo_accession in col 2)
awk 'NR==FNR { if (FNR>1 && $5=="PASS") pass[$1]=1; next }
     FNR==1 { print; next }
     ($2 in pass)' OFS='\t' "$starqc_report" "$uniqueAccList" > "$passAccList"

# Quick stats
total_unique=$(($(wc -l < "$uniqueAccList") - 1))
total_pass=$(($(wc -l < "$passAccList") - 1))
total_fail=$(( total_unique - total_pass ))
echo "STARQC gating: PASS=$total_pass  FAIL=$total_fail"
echo "PASS-only accession list saved to: $passAccList"

# If nothing passed, stop the pipeline early
if [ "$total_pass" -le 0 ]; then
    echo "ERROR: No samples passed alignment rate threshold. Exiting."
    exit 1
fi

##############################
# Sample Sex Estimation      #
##############################
echo "Alignment with STAR is complete for all samples. Computing sample sex."

ComputeSex_log_dir="$Logdir/ComputeSex"
mkdir -p "$ComputeSex_log_dir"

# Submit ComputeSex job
job_output=$(sbatch \
    --job-name="computeSex" \
    --export=PRJdir="$PRJdir",Logdir="$Logdir",AccList="$passAccList",BIOProjectID="$BIOProjectID",baseDir="$baseDir" \
    --output="$ComputeSex_log_dir/ComputeSex-%j.out" \
    --error="$ComputeSex_log_dir/ComputeSex-%j.err" \
    ComputeSex_v4.sh "$passAccList" "$BIOProjectID")

compute_sex_job_id=$(echo "$job_output" | awk '{print $4}')
echo "ComputeSex job submission output: $job_output"

if [[ -z "$compute_sex_job_id" || ! "$compute_sex_job_id" =~ ^[0-9]+$ ]]; then
    echo "Error: Invalid ComputeSex job ID: '$compute_sex_job_id'"
    exit 1
fi

echo "ComputeSex_v4.sh submitted with job ID: $compute_sex_job_id"

# Wait for the ComputeSex job to complete
wait_for_jobs "$compute_sex_job_id"
echo "ComputeSex job completed. Continuing to next step..."

##############################
# RSEM Reference Preparation #
##############################

Rref_log_dir="$Logdir/RSEM"
mkdir -p "$Rref_log_dir"

RSEM_REF_READY=1
for f in \
    rsemref.chrlist \
    rsemref.grp \
    rsemref.idx.fa \
    rsemref.n2g.idx.fa \
    rsemref.seq \
    rsemref.ti \
    rsemref.transcripts.fa
do
    [[ -s "$scratch_dir/$f" ]] || RSEM_REF_READY=0
done

if [[ "$RSEM_REF_READY" -eq 1 ]]; then
    echo "RSEM reference already exists skipping RSEMref_v4.sh"
    rsemRef_job_id=""
else
    echo "Submitting RSEMref_v4.sh to generate RSEM reference"

    rsemRef_job_id=$(sbatch \
        --export=BIOProjectID="$BIOProjectID" \
        --output="$Rref_log_dir/RSEMref-%j.out" \
        --error="$Rref_log_dir/RSEMref-%j.err" \
        RSEMref_v4.sh | awk '{print $4}')

    [[ -n "$rsemRef_job_id" ]] || {
        echo "ERROR: Failed to submit RSEMref_v4.sh"
        exit 1
    }

    echo "RSEM reference job ID: $rsemRef_job_id"

    wait_for_jobs "$rsemRef_job_id"

    if ! check_job_success "$rsemRef_job_id" "RSEM Reference"; then
        echo "STEP 2 FAILED RSEM reference generation failed"
        exit 1
    fi
fi

##############################
# RSEM Quantification        #
##############################
rsem_job_ids=()
while IFS=$'\t' read -r run geo_accession Tissue Strain Sex PMID GEOpath Title Sample_characteristics StrainInfo; do
    if [ "$run" == "Run" ]; then continue; fi

    unique_name="${Tissue}_${Strain}_${Sex}_${geo_accession}"
    sample_output_dir="$PRJdir/$geo_accession"
    sample_log_dir="$sample_output_dir/log_files/RSEM"
    mkdir -p "$sample_log_dir"

    echo "Submitting RSEM job for $geo_accession"

    if [ -z "$rsemRef_job_id" ]; then
        rsem_job_id=$(sbatch --job-name="${geo_accession}_RSEM72" \
                        --export=BIOProjectID="$BIOProjectID",Run="$run",geo_accession="$geo_accession",PRJdir="$PRJdir",unique_name="$unique_name",Logdir="$Logdir" \
                        --output="$sample_log_dir/RSEM-%j.out" \
                        --error="$sample_log_dir/RSEM-%j.err" \
                        RSEM_v5.sh | awk '{print $4}')
    else
        rsem_job_id=$(sbatch --dependency=afterok:$rsemRef_job_id \
                        --job-name="${geo_accession}_RSEM72" \
                        --export=BIOProjectID="$BIOProjectID",Run="$run",geo_accession="$geo_accession",PRJdir="$PRJdir",unique_name="$unique_name",Logdir="$Logdir" \
                        --output="$sample_log_dir/RSEM-%j.out" \
                        --error="$sample_log_dir/RSEM-%j.err" \
                        RSEM_v5.sh | awk '{print $4}')
    fi

    if [ $? -ne 0 ]; then
        echo "Error submitting RSEM job for $geo_accession"
        continue
    fi
    rsem_job_ids+=("$rsem_job_id")
done < "$passAccList"

# Wait for all RSEM sample jobs
echo "Waiting for RSEM jobs to complete..."
wait_for_jobs "${rsem_job_ids[@]}"

# Check RSEM job success
echo "Checking RSEM job completion status..."
for rsem_job_id in "${rsem_job_ids[@]}"; do
    if ! check_job_success "$rsem_job_id" "RSEM"; then
        echo "STEP 2 FAILED for $BIOProjectID - RSEM failed"
        exit 1
    fi
done
echo "All RSEM jobs completed successfully"


##############################
# Generate RSEM Matrix       #
##############################
matrix_log_dir="$Logdir/RSEMmatrix"
mkdir -p "$matrix_log_dir"

echo "Submitting RSEMmatrix job" >> "$log_file"

matrix_output=$(sbatch \
    --job-name="RSEMmatrix_$BIOProjectID" \
    --export=BIOProjectID="$BIOProjectID",PRJdir="$PRJdir",Logdir="$Logdir" \
    --output="$matrix_log_dir/RSEMmatrix-%j.out" \
    --error="$matrix_log_dir/RSEMmatrix-%j.err" \
    RSEMmatrix_v4.sh "$passAccList" "$BIOProjectID" 2>&1)

matrix_job_id=$(echo "$matrix_output" | grep -oP 'Submitted batch job \K\d+')

if [[ -z "$matrix_job_id" || ! "$matrix_job_id" =~ ^[0-9]+$ ]]; then
    echo "ERROR: Failed to submit RSEMmatrix job"
    echo "Output: $matrix_output"
    exit 1
fi

echo "RSEMmatrix job ID: $matrix_job_id" >> "$log_file"

echo "Waiting for RSEMmatrix job to complete..." >> "$log_file"
wait_for_jobs "$matrix_job_id"
if ! check_job_success "$matrix_job_id" "RSEM Matrix"; then
    echo "STEP 2 FAILED for $BIOProjectID - RSEM matrix generation failed"
    exit 1
fi
echo "RSEMmatrix job completed. Proceeding with BWjson generation..." >> "$log_file"

###################################
# Generate BigWig JSON files      #
###################################
declare -A bwjson_job_ids
declare -A seen_samples
failed_bwjson_samples=()

while IFS=$'\t' read -r run geo_accession Tissue Strain Sex PMID GEOpath Title Sample_characteristics StrainInfo || [[ -n "$run" ]]; do
    if [ "$run" == "Run" ]; then continue; fi
    [[ -z "$run" ]] && continue

    if [[ -n "${seen_samples[$geo_accession]:-}" ]]; then continue; fi
    seen_samples["$geo_accession"]=1

    unique_name="${Tissue}_${Strain}_${Sex}_${geo_accession}"
    sample_output_dir="$PRJdir/$geo_accession"
    sample_log_dir="$sample_output_dir/log_files/BWjson"

    mkdir -p "$sample_log_dir"

    echo "Submitting BWjson job for $geo_accession" >> "$log_file"

    if [[ -z "$matrix_job_id" || ! "$matrix_job_id" =~ ^[0-9]+$ ]]; then
        echo "Skipping BWjson for $geo_accession due to invalid matrix_job_id" >> "$log_file"
        continue
    fi

    # Export all necessary variables and submit job
    bwjson_job_ids["$geo_accession"]=$(
        BIOProjectID="$BIOProjectID" \
        Run="$run" \
        geo_accession="$geo_accession" \
        baseDir="$baseDir" \
        PRJdir="$PRJdir" \
        tissue="$Tissue" \
        strain="$Strain" \
        sex="$Sex" \
        title="$Title" \
        Sample_characteristics="$Sample_characteristics" \
        PMID="$PMID" \
        GEOpath="$GEOpath" \
        unique_name="$unique_name" \
        Logdir="$Logdir" \
        scratch_dir="$scratch_dir" \
        StrainInfo="$StrainInfo" \
        sbatch --dependency=afterok:$matrix_job_id \
            --job-name="BWjson_$geo_accession" \
            --export=ALL \
            --output="$sample_log_dir/BWjson-%j.out" \
            --error="$sample_log_dir/BWjson-%j.err" \
            BWjson_v7.sh | awk '{print $4}'
    )

    if [[ -z "${bwjson_job_ids["$geo_accession"]:-}" ]]; then
        echo "Error submitting BWjson job for $geo_accession" >> "$log_file"
    else
        echo "Submitted BWjson for $geo_accession with job ID: ${bwjson_job_ids["$geo_accession"]:-}" >> "$log_file"
    fi

done < "$passAccList"

echo "All BWjson generation jobs have been submitted." >> "$log_file"

# Wait for all BWjson jobs and check success
echo "Waiting for BWjson jobs to complete..."
for geo_accession in "${!bwjson_job_ids[@]}"; do
    jid="${bwjson_job_ids[$geo_accession]:-}"
    if [[ -n "$jid" ]]; then
        wait_for_jobs "$jid"
        if ! check_job_success "$jid" "BWjson ($geo_accession)"; then
            echo "WARNING: BWjson job failed for $geo_accession" >> "$log_file"
            failed_bwjson_samples+=("$geo_accession")
        fi
    fi
done

# Handle BWjson failures based on configuration
if [ ${#failed_bwjson_samples[@]} -gt 0 ]; then
    failed_list=$(printf '  - %s\n' "${failed_bwjson_samples[@]:-}")
    
    if [ "$BWJSON_REQUIRED" = true ]; then
        # CRITICAL: Fail the pipeline
        echo "ERROR: ${#failed_bwjson_samples[@]} samples had BWjson failures (CRITICAL)" >> "$log_file"
        send_notification "Pipeline FAILURE: $BIOProjectID - BWjson failed" \
            "BigWig JSON generation FAILED for $BIOProjectID:\n\nFAILED samples (${#failed_bwjson_samples[@]}):\n${failed_list}\n\nBWjson generation is configured as REQUIRED for this pipeline.\n\nCheck logs in:\n$PRJdir/<sample>/log_files/BWjson/\n\nCommon issues:\n- Wrong script path/version (check sbatch command)\n- Missing BigWig files from STAR\n- Disk space issues\n- Memory errors\n\nTo make BWjson optional, set BWJSON_REQUIRED=false at line 58 of the script."
        echo "STEP 2 FAILED for $BIOProjectID - BWjson generation failed"
        echo "Completed: $(date)"
        exit 1
    else
        # NON-CRITICAL: Send warning and continue
        echo "WARNING: ${#failed_bwjson_samples[@]} samples had BWjson failures (non-critical)" >> "$log_file"
        send_notification "Pipeline WARNING: $BIOProjectID - Some BWjson jobs failed" \
            "BigWig JSON generation results for $BIOProjectID:\n\nFAILED samples (${#failed_bwjson_samples[@]}):\n${failed_list}\n\nThese samples will not have visualization tracks in JBrowse.\nCore analysis (RSEM matrices) completed successfully.\n\nCheck logs in:\n$PRJdir/<sample>/log_files/BWjson/"
    fi
else
    echo "All BWjson jobs completed successfully"
fi

###############################################
# Submit JBrowse session builder job         #
###############################################
echo "Preparing to submit JBrowse session builder job..." >> "$log_file"

# Collect all BWjson job IDs into a dependency list
dep_ids=()
for geo_accession in "${!bwjson_job_ids[@]}"; do
    jid="${bwjson_job_ids[$geo_accession]:-}"
    if [[ -n "$jid" ]]; then
        dep_ids+=("$jid")
    fi
done

if ((${#dep_ids[@]} > 0)); then
    # Join job IDs with ':'
    dep_string=$(IFS=:; echo "${dep_ids[*]:-}")

    session_log_dir="$Logdir/JBrowseSession"
    mkdir -p "$session_log_dir"

    echo "Submitting JBrowseSession_v1.sh with dependency afterok:$dep_string" >> "$log_file"

    session_submit=$(sbatch \
        --dependency=afterok:$dep_string \
        --job-name="JBrowseSession_${BIOProjectID}" \
        --export=BIOProjectID="$BIOProjectID",PRJdir="$PRJdir",baseDir="$baseDir",Logdir="$Logdir" \
        --output="$session_log_dir/JBrowseSession-%j.out" \
        --error="$session_log_dir/JBrowseSession-%j.err" \
      JBrowseSession_v1.sh 2>&1)

    session_job_id=$(echo "$session_submit" | grep -oP 'Submitted batch job \K\d+')

    if [[ -z "$session_job_id" || ! "$session_job_id" =~ ^[0-9]+$ ]]; then
        echo "ERROR: Invalid JBrowseSession job ID: '$session_job_id'" >> "$log_file"
        echo "WARNING: JBrowse session will not be created" >> "$log_file"
    else
        echo "JBrowseSession_v1.sh submitted with job ID: $session_job_id" >> "$log_file"
        echo "JBrowse session JSON will be written to: ${baseDir}/${BIOProjectID}_jbrowse_session.json" >> "$log_file"

        # Wait for JBrowse session
        wait_for_jobs "$session_job_id"
        if ! check_job_success "$session_job_id" "JBrowse Session"; then
            echo "WARNING: JBrowse session generation failed (non-critical)" >> "$log_file"
            send_notification "Pipeline WARNING: $BIOProjectID - JBrowse session failed" \
                "JBrowse session generation failed for $BIOProjectID.\n\nThis is non-critical - all analysis data is complete.\n\nCheck logs at:\n$session_log_dir/JBrowseSession-$session_job_id.out\n$session_log_dir/JBrowseSession-$session_job_id.err"
        fi
    fi
else
    echo "Warning: No BWjson job IDs recorded; JBrowse session will not be generated." >> "$log_file"
    send_notification "Pipeline WARNING: $BIOProjectID - No JBrowse session (no BWjson)" \
        "No BWjson jobs succeeded for $BIOProjectID.\n\nJBrowse session cannot be created without BigWig tracks.\nCore analysis (alignment and quantification) completed successfully."
fi

echo "=========================================="
echo "STEP 2 COMPLETED SUCCESSFULLY for $BIOProjectID"
echo "Completed: $(date)"
echo "=========================================="

###############################################
# Automatic Scratch Space Cleanup            #
###############################################
if [[ "$KEEP_SCRATCH" == "false" ]]; then
    echo ""
    echo "Starting automatic scratch space cleanup..."

    # Calculate size before cleanup
    if [ -d "$scratch_dir" ]; then
        scratch_size_before=$(du -sh "$scratch_dir" 2>/dev/null | awk '{print $1}')
        echo "Scratch space before cleanup: $scratch_size_before"

        # Remove Project Scratch directory
        echo "Removing $scratch_dir"
        rm -rf "$scratch_dir"

        # Calculate size after cleanup
        if [ -d "$scratch_dir" ]; then
            scratch_size_after=$(du -sh "$scratch_dir" 2>/dev/null | awk '{print $1}')
            echo "Project scratch space after cleanup: $scratch_size_after"
        else
            echo "Scratch directory removed completely."
        fi

        echo "Scratch cleanup complete"
    else
        echo "Scratch directory not found, no cleanup needed."
    fi
else
    echo ""
    echo "Scratch space cleanup SKIPPED (KEEP_SCRATCH=true)"
    echo "Scratch directory preserved at: $scratch_dir"
fi

echo "=========================================="

# Create overall success marker
touch "${SUCCESS_DIR}/.step2_complete"

exit 0
