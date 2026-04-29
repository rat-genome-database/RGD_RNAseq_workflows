#!/usr/bin/env bash
#SBATCH --job-name=S2R
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=06:30:00
#SBATCH --account=YOUR_SLURM_ACCOUNT       # <-- replace with your SLURM account
#SBATCH --output=%x-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=YOUR_EMAIL@example.com  # <-- replace with your email
#SBATCH --error=regensj.err

#####################################################################################################################
# sex_json_regen_v2.sh
#
# Utility script: Regenerates the sex conflict report and per-sample BigWig JSON
# (JBrowse2 track definition) files for a BioProject. Use this when you need to
# re-run these steps independently of the full pipeline — for example, after
# manually correcting a sex assignment in the sex result file.
#
# Usage:
#   bash sex_json_regen_v2.sh <SRA_accession_list.txt> <BioProject_ID>
#   sbatch sex_json_regen_v2.sh <SRA_accession_list.txt> <BioProject_ID>
#
# Arguments:
#   <SRA_accession_list.txt>  Tab-delimited accession list used by the main pipeline.
#                             Must include a header row with: Run, geo_accession, Tissue,
#                             Strain, Sex, PMID, GEOpath, Title, Sample_characteristics, StrainInfo
#   <BioProject_ID>           GEO/BioProject accession (e.g. GSE70012). Used to locate
#                             the project directory and label output files.
#
# Important: Pass the accession list that reflects the final/corrected sample set so
# that regeneration does not fail due to unprocessed samples.
#
# Steps performed:
#   1. Run ConflictedSampleReport to regenerate the sex conflict report
#   2. Submit per-sample BWjson jobs (one per unique geo_accession)
#   3. Wait for all BWjson jobs to finish; report failures
#   4. Submit JBrowseSession job (dependent on BWjson completion) to rebuild the session JSON
#
# Dependencies:
#   - lib_v<N>.sh          : shared library sourced at the top (set LIB below)
#   - ConflictedSampleReport_v<N>.sh : sex conflict report generator (set SCRIPT_DIR below)
#   - BWjson_v<N>.sh       : per-sample BigWig JSON generator (set SCRIPT_DIR below)
#   - JBrowseSession_v<N>.sh : JBrowse2 session builder (set SCRIPT_DIR below)
#   - SLURM (squeue, sacct, sbatch)
#   - mail utility (optional, for failure notifications)
#
# Configuration:
#   Edit the CONFIGURATION block below before running.
#####################################################################################################################

set -euo pipefail
set -x

# ---------------------------------------------------------------------------
# CONFIGURATION — edit these variables before running
# ---------------------------------------------------------------------------
# Directory containing pipeline helper scripts and the shared library
SCRIPT_DIR="/path/to/your/pipeline/scripts"   # <-- replace

# Shared library filename (must exist inside SCRIPT_DIR)
LIB_NAME="lib_v10.sh"                          # <-- replace with your current lib version

# Root directory under which all BioProject data directories live
# Final project path will be: DATA_BASE/<BioProject_ID>/
DATA_BASE="/path/to/your/expression/data"      # <-- replace

# Scratch filesystem base for temporary index files
SCRATCH_BASE="/scratch/your/scratch/area"      # <-- replace

# Reference genome directory
REFdir="/path/to/your/reference/genome"        # <-- replace

# Notification email (used by send_notification; ignored in orchestrated mode)
NOTIFY_EMAIL="YOUR_EMAIL@example.com"          # <-- replace

# Script names for sub-steps (update version numbers as needed)
CONFLICT_SCRIPT="ConflictedSampleReport_v4.sh" # <-- replace with your current version
BWJSON_SCRIPT="BWjson_v7.sh"                   # <-- replace with your current version
SESSION_SCRIPT="JBrowseSession_v1.sh"          # <-- replace with your current version

# Set to "true" to make BWjson failures abort the pipeline; "false" to warn and continue
BWJSON_REQUIRED=true
# ---------------------------------------------------------------------------

LIB="${SCRIPT_DIR}/${LIB_NAME}"

if [[ ! -f "$LIB" ]]; then
    echo "FATAL: ${LIB_NAME} not found at $LIB" >&2
    exit 1
fi

source "$LIB"

# Orchestrated mode disables duplicate email notifications when called from a parent script
ORCHESTRATED_MODE="${ORCHESTRATED_MODE:-false}"
KEEP_SCRATCH="${KEEP_SCRATCH:-false}"

# --- Argument parsing ---
if [ $# -ne 2 ]; then
    echo "Usage: $0 <SRA_accession_list.txt> <BioProject_ID>"
    exit 1
fi

AccList=$1
BIOProjectID=$2

# --- Validate inputs ---
if [[ ! -f "$AccList" ]]; then
    echo "ERROR: Accession list file not found: $AccList"
    exit 1
fi

if [[ "$BIOProjectID" =~ [^a-zA-Z0-9_-] ]]; then
    echo "ERROR: Invalid characters in BioProject ID. Only alphanumeric, underscore, and hyphen allowed."
    exit 1
fi

# --- Derive paths from configuration ---
myDir="${DATA_BASE}"
PRJdir="${myDir}/${BIOProjectID}/reads_fastq"
Logdir="${myDir}/${BIOProjectID}/log_files"
scratch_dir="${SCRATCH_BASE}/${BIOProjectID}_$$"
INDEX_DIR="${scratch_dir}/RefIndex"
baseDir="${myDir}/${BIOProjectID}"

# --- Status tracking ---
SUCCESS_DIR="$baseDir/.status"
mkdir -p "$SUCCESS_DIR"

# --- Create directories ---
mkdir -p "$PRJdir"
mkdir -p "$Logdir"
mkdir -p "$scratch_dir"

# --- Redirect output to log ---
log_file="${Logdir}/${BIOProjectID}_regenSexConflictandjson.out"
exec 1>"$log_file" 2>&1

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

# Wait for a list of SLURM job IDs to leave the queue
wait_for_jobs() {
    local ids=("$@")

    if [ ${#ids[@]} -eq 0 ]; then
        ids=("${job_ids[@]:-}")
    fi

    for job_id in "${ids[@]:-}"; do
        [ -z "$job_id" ] && continue
        echo "Waiting for job: $job_id"
        while squeue -u "$USER" -h -o "%A" | grep -qx "$job_id"; do
            sleep 5
        done
    done
}

# Return 0 if the job completed successfully, 1 otherwise
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

# Send an email notification; suppressed when running in orchestrated mode
send_notification() {
    local subject="$1"
    local message="$2"

    if [[ "$ORCHESTRATED_MODE" == "true" ]]; then
        echo "NOTIFICATION (suppressed in orchestrated mode): $subject"
        return 0
    fi

    echo -e "$message" | mail -s "$subject" "$NOTIFY_EMAIL" 2>/dev/null || true
    echo "NOTIFICATION SENT: $subject"
}

# ---------------------------------------------------------------------------
# STEP 1: Regenerate sex conflict report
# ---------------------------------------------------------------------------
sh "${SCRIPT_DIR}/${CONFLICT_SCRIPT}" "${BIOProjectID}"

# ---------------------------------------------------------------------------
# STEP 2: Submit per-sample BigWig JSON generation jobs
# ---------------------------------------------------------------------------
declare -A bwjson_job_ids
declare -A seen_samples
failed_bwjson_samples=()

while IFS=$'\t' read -r run geo_accession Tissue Strain Sex PMID GEOpath Title Sample_characteristics StrainInfo || [[ -n "$run" ]]; do
    if [ "$run" == "Run" ]; then continue; fi
    [[ -z "$run" ]] && continue

    # Process only the first run per geo_accession (unique sample)
    if [[ -n "${seen_samples[$geo_accession]:-}" ]]; then continue; fi
    seen_samples["$geo_accession"]=1

    unique_name="${Tissue}_${Strain}_${Sex}_${geo_accession}"
    sample_output_dir="$PRJdir/$geo_accession"
    sample_log_dir="$sample_output_dir/log_files/BWjson"

    mkdir -p "$sample_log_dir"

    echo "Submitting BWjson job for $geo_accession" >> "$log_file"

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
            "${SCRIPT_DIR}/${BWJSON_SCRIPT}" | awk '{print $4}'
    )

    # If matrix job dependency is needed, replace the sbatch block above with:
    # if [[ -n "${matrix_job_id:-}" && "$matrix_job_id" =~ ^[0-9]+$ ]]; then
    #     bwjson_job_ids["$geo_accession"]=$(
    #         ... (same exports) ...
    #         sbatch --dependency=afterok:$matrix_job_id \
    #             ...
    #     )
    # else
    #     echo "WARNING: Skipping BWjson for $geo_accession - matrix_job_id not set" >> "$log_file"
    #     continue
    # fi

    if [[ -z "${bwjson_job_ids["$geo_accession"]:-}" ]]; then
        echo "Error submitting BWjson job for $geo_accession" >> "$log_file"
    else
        echo "Submitted BWjson for $geo_accession: job ID ${bwjson_job_ids["$geo_accession"]:-}" >> "$log_file"
    fi

done < "$AccList"

echo "All BWjson generation jobs submitted." >> "$log_file"

# ---------------------------------------------------------------------------
# Wait for BWjson jobs and report failures
# ---------------------------------------------------------------------------
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

if [ ${#failed_bwjson_samples[@]} -gt 0 ]; then
    failed_list=$(printf '  - %s\n' "${failed_bwjson_samples[@]:-}")

    if [ "$BWJSON_REQUIRED" = true ]; then
        echo "ERROR: ${#failed_bwjson_samples[@]} samples had BWjson failures (CRITICAL)" >> "$log_file"
        send_notification "Pipeline FAILURE: $BIOProjectID - BWjson failed" \
            "BigWig JSON generation FAILED for $BIOProjectID:\n\nFailed samples (${#failed_bwjson_samples[@]}):\n${failed_list}\n\nBWjson is configured as REQUIRED.\n\nCheck logs in:\n$PRJdir/<sample>/log_files/BWjson/\n\nCommon causes:\n- Wrong script path/version\n- Missing BigWig files from STAR\n- Disk space or memory errors\n\nTo make BWjson optional, set BWJSON_REQUIRED=false in the CONFIGURATION block."
        echo "STEP 2 FAILED for $BIOProjectID - BWjson generation failed"
        echo "Completed: $(date)"
        exit 1
    else
        echo "WARNING: ${#failed_bwjson_samples[@]} samples had BWjson failures (non-critical)" >> "$log_file"
        send_notification "Pipeline WARNING: $BIOProjectID - Some BWjson jobs failed" \
            "BigWig JSON generation results for $BIOProjectID:\n\nFailed samples (${#failed_bwjson_samples[@]}):\n${failed_list}\n\nThese samples will not have visualization tracks in JBrowse.\nCore analysis (RSEM matrices) completed successfully.\n\nCheck logs in:\n$PRJdir/<sample>/log_files/BWjson/"
    fi
else
    echo "All BWjson jobs completed successfully."
fi

# ---------------------------------------------------------------------------
# STEP 3: Submit JBrowse session builder (dependent on all BWjson jobs)
# ---------------------------------------------------------------------------
echo "Preparing to submit JBrowse session builder job..." >> "$log_file"

dep_ids=()
for geo_accession in "${!bwjson_job_ids[@]}"; do
    jid="${bwjson_job_ids[$geo_accession]:-}"
    if [[ -n "$jid" ]]; then
        dep_ids+=("$jid")
    fi
done

if ((${#dep_ids[@]} > 0)); then
    dep_string=$(IFS=:; echo "${dep_ids[*]:-}")

    session_log_dir="$Logdir/JBrowseSession"
    mkdir -p "$session_log_dir"

    echo "Submitting ${SESSION_SCRIPT} with dependency afterok:$dep_string" >> "$log_file"

    session_submit=$(sbatch \
        --dependency=afterok:$dep_string \
        --job-name="JBrowseSession_${BIOProjectID}" \
        --export=BIOProjectID="$BIOProjectID",PRJdir="$PRJdir",baseDir="$baseDir",Logdir="$Logdir" \
        --output="$session_log_dir/JBrowseSession-%j.out" \
        --error="$session_log_dir/JBrowseSession-%j.err" \
        "${SCRIPT_DIR}/${SESSION_SCRIPT}" 2>&1)

    session_job_id=$(echo "$session_submit" | grep -oP 'Submitted batch job \K\d+')

    if [[ -z "$session_job_id" || ! "$session_job_id" =~ ^[0-9]+$ ]]; then
        echo "ERROR: Invalid JBrowseSession job ID: '$session_job_id'" >> "$log_file"
        echo "WARNING: JBrowse session will not be created" >> "$log_file"
    else
        echo "JBrowse session job submitted: $session_job_id" >> "$log_file"
        echo "Output will be written to: ${baseDir}/${BIOProjectID}_jbrowse_session.json" >> "$log_file"

        wait_for_jobs "$session_job_id"
        if ! check_job_success "$session_job_id" "JBrowse Session"; then
            echo "WARNING: JBrowse session generation failed (non-critical)" >> "$log_file"
            send_notification "Pipeline WARNING: $BIOProjectID - JBrowse session failed" \
                "JBrowse session generation failed for $BIOProjectID.\n\nThis is non-critical — all analysis data is complete.\n\nCheck logs at:\n$session_log_dir/JBrowseSession-$session_job_id.out\n$session_log_dir/JBrowseSession-$session_job_id.err"
        fi
    fi
else
    echo "WARNING: No BWjson job IDs recorded; JBrowse session will not be generated." >> "$log_file"
    send_notification "Pipeline WARNING: $BIOProjectID - No JBrowse session (no BWjson jobs)" \
        "No BWjson jobs succeeded for $BIOProjectID.\n\nJBrowse session cannot be created without BigWig tracks.\nCore analysis (alignment and quantification) completed successfully."
fi

echo "=========================================="
echo "STEP 2 COMPLETED SUCCESSFULLY for $BIOProjectID"
echo "Completed: $(date)"
echo "=========================================="

exit 0
