#!/usr/bin/env bash
# Shared helpers for Single end GRCr8 bulk RNA-seq pipelines (v10)

set -eo pipefail

log() {
    echo "[$(date '+%F %T')] $*" >&2
}

require_nonempty() {
    local f="$1"
    [[ -s "$f" ]] || {
        log "ERROR: expected non-empty file: $f"
        return 1
    }
}

submit_job() {
    local jobfile="$1"
    shift
    local jid
    jid=$(sbatch "$jobfile" "$@" | awk '{print $4}')
    echo "$jid"
}

wait_for_job_final() {
    local jid="$1"
    local state
    state=$(sacct -j "$jid" --format=State --noheader | tail -n1 | tr -d ' ')
    [[ "$state" == "COMPLETED" ]] || {
        log "Job $jid failed with state $state"
        return 1
    }
}

