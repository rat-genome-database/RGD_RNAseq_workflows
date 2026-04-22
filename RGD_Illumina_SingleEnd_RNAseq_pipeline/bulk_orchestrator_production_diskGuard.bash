#!/usr/bin/env bash

###############################################################################
# USER CONFIGURATION
# Update the variables below for your system before running.
###############################################################################
# Directory containing all pipeline scripts
SCRIPT_DIR="/path/to/RGD_Illumina_SingleEnd_RNAseq_pipeline"
# Root directory for GEO project output data
BASE_DATA_DIR="/path/to/data/expression/GEO"
# Scratch filesystem root for temporary files
SCRATCH_DIR="/path/to/scratch"
###############################################################################

###############################################################################
# bulk_orchestrator_production_diskGuard.bash
# Hybrid scheduler with robust input parsing
###############################################################################

set -uo pipefail

###############################################################################
# CONFIG
###############################################################################

## 12 Feb 2026 update PATH

LIB="${SCRIPT_DIR}/lib_v10.sh"
STEP1_SCRIPT="${SCRIPT_DIR}/run_SRA2QC_SE_v1.bash"
STEP2_SCRIPT="${SCRIPT_DIR}/run_RNApipeline_SE_diskGuard_v1.bash"


MAX_CONCURRENT_SMALL="${MAX_CONCURRENT_SMALL:-4}"
SMALL_PROJECT_THRESHOLD="${SMALL_PROJECT_THRESHOLD:-20}"
POLL_INTERVAL="${POLL_INTERVAL:-60}"

# Scratch disk guard
SCRATCH_WARN_PCT="${SCRATCH_WARN_PCT:-85}"   # log a warning above this %
SCRATCH_HOLD_PCT="${SCRATCH_HOLD_PCT:-90}"   # pause new submissions above this %
SCRATCH_POLL_WAIT="${SCRATCH_POLL_WAIT:-300}" # seconds to wait before re-checking when held

###############################################################################
# PRE-FLIGHT CHECKS
###############################################################################

for f in "$LIB" "$STEP1_SCRIPT" "$STEP2_SCRIPT"; do
    [[ -f "$f" ]] || { echo "ERROR: missing $f" >&2; exit 1; }
done

source "$LIB"

[[ $# -eq 1 ]] || { echo "Usage: $0 projectlist.txt" >&2; exit 1; }
PROJECT_LIST="$1"
[[ -f "$PROJECT_LIST" ]] || { echo "ERROR: project list not found" >&2; exit 1; }

###############################################################################
# INTERNAL STATE
###############################################################################

declare -A PROJECT_ACCLIST=()
declare -A PROJECT_READLEN=()
declare -A PROJECT_SAMPLE_COUNT=()
declare -A PROJECT_SIZE_CLASS=()
declare -A STEP1_JOBS=()
declare -A STEP2_JOBS=()
declare -A PROJECT_STATUS=()
PROJECTS_TO_PROCESS=()

###############################################################################
# HELPERS
###############################################################################

count_samples_in_acclist() {
    local acclist_path="$1"

    [[ -f "$acclist_path" ]] || { echo "0"; return 1; }
    [[ -r "$acclist_path" ]] || { echo "0"; return 1; }

    local count
    count=$(awk -F'\t' '
        NR > 1 && !/^#/ && $2 != "" && $2 !~ /^geo_accession$/ {
            gsm[$2] = 1
        }
        END {
            print length(gsm) + 0
        }
    ' "$acclist_path" 2>/dev/null)

    [[ "$count" =~ ^[0-9]+$ ]] || count=0
    echo "$count"
}

count_active_small_projects() {
    local count=0
    for project in "${!PROJECT_STATUS[@]}"; do
        [[ "${PROJECT_SIZE_CLASS[$project]}" == "small" ]] || continue
        case "${PROJECT_STATUS[$project]}" in
            step1_running|step1_done|step2_running) count=$((count + 1)) ;;
        esac
    done
    echo "$count"
}

count_active_projects() {
    local count=0
    for project in "${!PROJECT_STATUS[@]}"; do
        case "${PROJECT_STATUS[$project]}" in
            step1_running|step1_done|step2_running) count=$((count + 1)) ;;
        esac
    done
    echo "$count"
}

has_active_large_project() {
    for project in "${!PROJECT_STATUS[@]}"; do
        [[ "${PROJECT_SIZE_CLASS[$project]}" == "large" ]] || continue
        case "${PROJECT_STATUS[$project]}" in
            step1_running|step1_done|step2_running) echo "true"; return ;;
        esac
    done
    echo "false"
}

check_job_status() {
    local jid="$1"
    squeue -j "$jid" -h &>/dev/null && { echo "RUNNING"; return; }
    local state=$(sacct -j "$jid" --format=State --noheader 2>/dev/null | head -1 | tr -d ' ')
    case "$state" in
        COMPLETED) echo "COMPLETED" ;;
        FAILED|CANCELLED|TIMEOUT|OUT_OF_MEMORY) echo "FAILED" ;;
        *) echo "UNKNOWN" ;;
    esac
}

###############################################################################
# DISK GUARD HELPERS
###############################################################################

# Returns the integer % used for the filesystem containing SCRATCH_DIR.
# Falls back to 0 on any error so the pipeline is never blocked by a df failure.
scratch_pct_used() {
    local pct
    pct=$(df --output=pcent "$SCRATCH_DIR" 2>/dev/null | tail -1 | tr -d ' %')
    [[ "$pct" =~ ^[0-9]+$ ]] || pct=0
    echo "$pct"
}

# Logs current usage and returns 0 (safe) or 1 (hold).
# Also emits a warning-level log when between WARN and HOLD thresholds.
scratch_is_safe() {
    local pct
    pct=$(scratch_pct_used)

    if [[ $pct -ge $SCRATCH_HOLD_PCT ]]; then
        log "DISK HOLD: scratch ${pct}% used (>= ${SCRATCH_HOLD_PCT}% threshold) — pausing new submissions"
        return 1
    elif [[ $pct -ge $SCRATCH_WARN_PCT ]]; then
        log "DISK WARN: scratch ${pct}% used (>= ${SCRATCH_WARN_PCT}% warn threshold)"
        return 0
    else
        log "DISK OK:   scratch ${pct}% used"
        return 0
    fi
}

submit_step1() {
    local acclist="$1" project="$2" readlen="$3"
    local sample_count="${PROJECT_SAMPLE_COUNT[$project]}"
    local size_class="${PROJECT_SIZE_CLASS[$project]}"

    log "Submitting STEP 1 for $project ($sample_count samples, $size_class)"

    local jid
    jid=$(sbatch --parsable --job-name="STEP1_${project}" "$STEP1_SCRIPT" "$acclist" "$project")

    STEP1_JOBS["$project"]="$jid"
    PROJECT_STATUS["$project"]="step1_running"

    log "STEP 1 submitted for $project (job: $jid)"
}

submit_step2() {
    local acclist="$1" project="$2" readlen="$3"
    local sample_count="${PROJECT_SAMPLE_COUNT[$project]}"
    local size_class="${PROJECT_SIZE_CLASS[$project]}"

    log "Submitting STEP 2 for $project ($sample_count samples, $size_class)"

    local jid
    jid=$(sbatch --parsable --job-name="STEP2_${project}" \
        --export=ORCHESTRATED_MODE=true,KEEP_SCRATCH=false \
        "$STEP2_SCRIPT" "$acclist" "$project" "$readlen")

    STEP2_JOBS["$project"]="$jid"
    PROJECT_STATUS["$project"]="step2_running"

    log "STEP 2 submitted for $project (job: $jid)"
}

###############################################################################
# LOAD PROJECT LIST
###############################################################################

log "Starting bulk_orchestrator_production_diskGuard.bash(hybrid mode)"
log "Small project threshold: ≤${SMALL_PROJECT_THRESHOLD} samples"
log "Max concurrent small projects: $MAX_CONCURRENT_SMALL"
log "Large projects run in ISOLATION"
log "Poll interval: ${POLL_INTERVAL}s"
log "Scratch dir: $SCRATCH_DIR"
log "Scratch warn threshold: ${SCRATCH_WARN_PCT}%"
log "Scratch hold threshold: ${SCRATCH_HOLD_PCT}%"
log ""
log "Loading and analyzing projects..."
log "=========================================="

total_small=0
total_large=0

# Read project list using awk to properly handle tabs/spaces
while IFS= read -r line; do
    # Skip empty lines and comments
    [[ -z "$line" ]] && continue
    [[ "$line" =~ ^[[:space:]]*# ]] && continue

    # Parse using awk to handle any whitespace
    read -r ACCLIST PROJECT READLEN <<< "$(echo "$line" | awk '{print $1, $2, $3}')"

    [[ -z "$ACCLIST" ]] && continue
    [[ -z "$PROJECT" ]] && continue
    [[ -z "$READLEN" ]] && READLEN="150"  # default

    log "Processing: $PROJECT"
    log "  AccList: $ACCLIST"

    # Store metadata
    PROJECT_ACCLIST["$PROJECT"]="$ACCLIST"
    PROJECT_READLEN["$PROJECT"]="$READLEN"
    PROJECTS_TO_PROCESS+=("$PROJECT")

    # Count samples
    sample_count=$(count_samples_in_acclist "$ACCLIST")

    if [[ $sample_count -eq 0 ]]; then
        log "  ERROR: Found 0 samples - SKIPPING"
        PROJECT_STATUS["$PROJECT"]="failed"
        PROJECT_SAMPLE_COUNT["$PROJECT"]=0
        PROJECT_SIZE_CLASS["$PROJECT"]="unknown"
        continue
    fi

    PROJECT_SAMPLE_COUNT["$PROJECT"]=$sample_count
    log "  Samples: $sample_count"

    # Classify by size
    if [[ $sample_count -le $SMALL_PROJECT_THRESHOLD ]]; then
        PROJECT_SIZE_CLASS["$PROJECT"]="small"
        total_small=$((total_small + 1))
        log "  Size: SMALL"
    else
        PROJECT_SIZE_CLASS["$PROJECT"]="large"
        total_large=$((total_large + 1))
        log "  Size: LARGE"
    fi

    # Check completion status
    STATUS_DIR="${BASE_DATA_DIR}/${PROJECT}/.status"

    log "  Creating status directory: $STATUS_DIR"
    if ! mkdir -p "$STATUS_DIR" 2>/dev/null; then
        log "  ERROR: Failed to create status directory"
        log "  Continuing anyway..."
    fi

    if [[ -f "${STATUS_DIR}/.step2_complete" ]]; then
        PROJECT_STATUS["$PROJECT"]="complete"
        log "  Status: Already complete (skipping)"
    elif [[ -f "${STATUS_DIR}/.step1_complete" ]]; then
        PROJECT_STATUS["$PROJECT"]="step1_done"
        log "  Status: Step 1 done, Step 2 pending"
    else
        PROJECT_STATUS["$PROJECT"]="pending"
        log "  Status: Pending (needs both steps)"
    fi
    log ""

done < "$PROJECT_LIST"

log "=========================================="
log "Loaded ${#PROJECTS_TO_PROCESS[@]} projects: $total_small small, $total_large large"
log "=========================================="
log ""

###############################################################################
# MAIN PROCESSING LOOP
###############################################################################

all_complete=false

while [[ "$all_complete" == "false" ]]; do

    total_active=$(count_active_projects)
    active_small=$(count_active_small_projects)
    large_active=$(has_active_large_project)

    log "Status: $total_active active ($active_small small) | Large active: $large_active"

    #########################################
    # DISK GUARD: check scratch before starting anything new
    #########################################
    if ! scratch_is_safe; then
        log "Scratch at or above ${SCRATCH_HOLD_PCT}% — skipping new submissions this cycle."
        log "Active jobs will continue; waiting ${SCRATCH_POLL_WAIT}s for space to free up..."
        sleep "$SCRATCH_POLL_WAIT"
        # Still need to fall through to phases 2-5 to collect completions,
        # but Phase 1 is gated by disk_ok below.
    fi
    disk_ok=$(scratch_is_safe && echo "true" || echo "false")

    #########################################
    # PHASE 1: Start new projects
    #########################################
    if [[ "$disk_ok" == "false" ]]; then
        log "PHASE 1 skipped — scratch disk above hold threshold"
    else
    for project in "${PROJECTS_TO_PROCESS[@]}"; do
        [[ "${PROJECT_STATUS[$project]}" != "pending" ]] && continue

        # Recalculate active counts on each iteration
        total_active=$(count_active_projects)
        active_small=$(count_active_small_projects)
        large_active=$(has_active_large_project)

        size_class="${PROJECT_SIZE_CLASS[$project]}"
        sample_count="${PROJECT_SAMPLE_COUNT[$project]}"
        acclist="${PROJECT_ACCLIST[$project]}"
        readlen="${PROJECT_READLEN[$project]}"

        # RULE 1: If large project active, wait
        if [[ "$large_active" == "true" ]]; then
            log "Large project active - waiting"
            break
        fi

        # RULE 2: Starting large project - need clear system
        if [[ "$size_class" == "large" ]]; then
            if [[ $total_active -gt 0 ]]; then
                log "Cannot start large $project - waiting for active projects"
                continue
            else
                log "Starting large $project ($sample_count samples) in ISOLATION"
                submit_step1 "$acclist" "$project" "$readlen"
                break
            fi
        fi

        # RULE 3: Starting small project - check capacity
        if [[ "$size_class" == "small" ]]; then
            if [[ $active_small -lt $MAX_CONCURRENT_SMALL ]]; then
                log "Starting small $project ($sample_count samples)"
                submit_step1 "$acclist" "$project" "$readlen"
            else
                log "At capacity: $active_small/$MAX_CONCURRENT_SMALL small projects"
            fi
        fi
    done
    fi  # end disk_ok gate

    #########################################
    # PHASE 2: Check Step 1 completions
    #########################################
    for project in "${!STEP1_JOBS[@]}"; do
        jid="${STEP1_JOBS[$project]}"
        [[ -z "$jid" ]] && continue

        job_status=$(check_job_status "$jid")

        if [[ "$job_status" == "COMPLETED" ]]; then
            STATUS_DIR="${BASE_DATA_DIR}/${project}/.status"
            if [[ -f "${STATUS_DIR}/.step1_complete" ]]; then
                log "STEP 1 completed for $project"
                PROJECT_STATUS["$project"]="step1_done"
                STEP1_JOBS["$project"]=""
            else
                log "STEP 1 marker missing for $project"
                PROJECT_STATUS["$project"]="failed"
                STEP1_JOBS["$project"]=""
            fi
        elif [[ "$job_status" == "FAILED" ]]; then
            log "STEP 1 failed for $project"
            PROJECT_STATUS["$project"]="failed"
            STEP1_JOBS["$project"]=""
        fi
    done

    #########################################
    # PHASE 3: Submit Step 2
    #########################################
    for project in "${PROJECTS_TO_PROCESS[@]}"; do
        if [[ "${PROJECT_STATUS[$project]}" == "step1_done" ]]; then
            acclist="${PROJECT_ACCLIST[$project]}"
            readlen="${PROJECT_READLEN[$project]}"
            submit_step2 "$acclist" "$project" "$readlen"
        fi
    done

    #########################################
    # PHASE 4: Check Step 2 completions
    #########################################
    for project in "${!STEP2_JOBS[@]}"; do
        jid="${STEP2_JOBS[$project]}"
        [[ -z "$jid" ]] && continue

        job_status=$(check_job_status "$jid")

        if [[ "$job_status" == "COMPLETED" ]]; then
            STATUS_DIR="${BASE_DATA_DIR}/${project}/.status"
            if [[ -f "${STATUS_DIR}/.step2_complete" ]]; then
                sample_count="${PROJECT_SAMPLE_COUNT[$project]}"
                log "STEP 2 completed for $project ($sample_count samples)"
                log "PROJECT $project COMPLETE"
                PROJECT_STATUS["$project"]="complete"
                STEP2_JOBS["$project"]=""
            else
                log "STEP 2 marker missing for $project"
                PROJECT_STATUS["$project"]="failed"
                STEP2_JOBS["$project"]=""
            fi
        elif [[ "$job_status" == "FAILED" ]]; then
            log "STEP 2 failed for $project"
            PROJECT_STATUS["$project"]="failed"
            STEP2_JOBS["$project"]=""
        fi
    done

    #########################################
    # PHASE 5: Check completion
    #########################################
    all_complete=true
    for project in "${PROJECTS_TO_PROCESS[@]}"; do
        status="${PROJECT_STATUS[$project]}"
        if [[ "$status" != "complete" && "$status" != "failed" ]]; then
            all_complete=false
            break
        fi
    done

    if [[ "$all_complete" == "false" ]]; then
        log "Sleeping ${POLL_INTERVAL}s..."
        log ""
        sleep "$POLL_INTERVAL"
    fi
done

###############################################################################
# FINAL REPORT
###############################################################################

log "=========================================="
log "Final Results:"
log "=========================================="

failed_projects=()
successful_projects=()
successful_small=0
successful_large=0

for project in "${PROJECTS_TO_PROCESS[@]}"; do
    status="${PROJECT_STATUS[$project]}"
    sample_count="${PROJECT_SAMPLE_COUNT[$project]}"
    size_class="${PROJECT_SIZE_CLASS[$project]}"

    if [[ "$status" == "complete" ]]; then
        log "$project: SUCCESS ($sample_count samples, $size_class)"
        successful_projects+=("$project")
        if [[ "$size_class" == "small" ]]; then
            successful_small=$((successful_small + 1))
        else
            successful_large=$((successful_large + 1))
        fi
    elif [[ "$status" == "failed" ]]; then
        log "$project: FAILED ($sample_count samples, $size_class)"
        failed_projects+=("$project")
    else
        log "? $project: $status"
    fi
done

log "=========================================="
log "Summary: ${#successful_projects[@]} succeeded ($successful_small small, $successful_large large)"
log "         ${#failed_projects[@]} failed"
log "=========================================="

[[ ${#failed_projects[@]} -gt 0 ]] && exit 1 || exit 0
