#!/bin/bash
#SBATCH --job-name=S2R
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=06:30:00
#SBATCH --account=your-slurm-account
#SBATCH --output=%x-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=your@email.edu
#SBATCH --error=regensj.err

###############################################################################
# USER CONFIGURATION
# Update the variables below for your system before running.
###############################################################################
# Directory containing all pipeline scripts
SCRIPT_DIR="/path/to/your/pipeline/scripts"
# Your home/base directory
myDir="/path/to/your/home"
# Scratch filesystem mount point
SCRATCH_BASE="/path/to/your/scratch/mount"
###############################################################################


#updated 9 Feb 2026
#ensure to submit the passed accession list so regeneration does not fail due to unprocessed samples

# Enable strict error handling and debugging
set -euo pipefail
set -x

LIB="${SCRIPT_DIR}/lib_v10.sh"

if [[ ! -f "$LIB" ]]; then
    echo "FATAL: lib_v10.sh not found at $LIB" >&2
    exit 1
fi

source "$LIB"

# Configuration: Orchestrated mode disables duplicate notifications
ORCHESTRATED_MODE="${ORCHESTRATED_MODE:-false}"
KEEP_SCRATCH="${KEEP_SCRATCH:-false}"

# Check command line arguments
if [ $# -ne 2 ]; then
    echo "Usage: $0 <SRA_accession_list.txt> <BioProject_ID> "
    exit 1
fi

AccList=$1
BIOProjectID=$2

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


PRJdir="$myDir/data/expression/GEO/$BIOProjectID/reads_fastq"
Logdir="$myDir/data/expression/GEO/$BIOProjectID/log_files"
# FIX: Define scratch_dir before using it
INDEX_DIR="$scratch_dir/RefIndex"
baseDir="$myDir/data/expression/GEO/$BIOProjectID"

# Configuration: Set to "true" to make BWjson failures critical, "false" for warnings only
BWJSON_REQUIRED=true

# Status tracking
SUCCESS_DIR="$baseDir/.status"
mkdir -p "$SUCCESS_DIR"

# Create directories
mkdir -p "$PRJdir"
mkdir -p "$Logdir"
mkdir -p "$scratch_dir"

# Redirect output to log file
log_file="${Logdir}/${BIOProjectID}_regenSexConflictandjson.out"
exec 1>"$log_file" 2>&1

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

    echo -e "$message" | mail -s "$subject" your@email.edu 2>/dev/null || true
    echo "NOTIFICATION SENT: $subject"
}


##################################
# Run Sample Sex Conflict Script #
##################################
sh "${SCRIPT_DIR}/ConflictedSampleReport_v4.sh" ${BIOProjectID}
###################################
# Generate BigWig JSON files      #
###################################
declare -A bwjson_job_ids
declare -A seen_samples
failed_bwjson_samples=()

# FIX: Changed $passAccList to $AccList (the actual defined variable)
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

    # FIX: Removed the matrix_job_id check - if you need matrix jobs, define matrix_job_id earlier
    # If BWjson jobs don't depend on matrix jobs, submit them without dependency
    # If they DO depend on matrix jobs, you need to set matrix_job_id before this section

    # OPTION 1: No dependency (if matrix jobs aren't needed)
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
        sbatch --job-name="BWjson_$geo_accession" \
            --export=ALL \
            --output="$sample_log_dir/BWjson-%j.out" \
            --error="$sample_log_dir/BWjson-%j.err" \
            BWjson_v7.sh | awk '{print $4}'
    )

    # OPTION 2: If you need matrix job dependency, uncomment this instead:
    # if [[ -n "${matrix_job_id:-}" && "$matrix_job_id" =~ ^[0-9]+$ ]]; then
    #     bwjson_job_ids["$geo_accession"]=$(
    #         BIOProjectID="$BIOProjectID" \
    #         ... (same exports as above) ...
    #         sbatch --dependency=afterok:$matrix_job_id \
    #             --job-name="BWjson_$geo_accession" \
    #             --export=ALL \
    #             --output="$sample_log_dir/BWjson-%j.out" \
    #             --error="$sample_log_dir/BWjson-%j.err" \
    #             BWjson_v7.sh | awk '{print $4}'
    #     )
    # else
    #     echo "WARNING: Skipping BWjson for $geo_accession - matrix_job_id not set or invalid" >> "$log_file"
    #     continue
    # fi

    if [[ -z "${bwjson_job_ids["$geo_accession"]:-}" ]]; then
        echo "Error submitting BWjson job for $geo_accession" >> "$log_file"
    else
        echo "Submitted BWjson for $geo_accession with job ID: ${bwjson_job_ids["$geo_accession"]:-}" >> "$log_file"
    fi

done < "$AccList"

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


exit 0
