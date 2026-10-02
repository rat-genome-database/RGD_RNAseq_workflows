#!/usr/bin/env bash
# Re-call library prep for every sample of an already-processed project with libprep_call.sh
# (same script STAR_bigwig3.sh uses), from files kept in the project directory. No realignment.
# Then run ConflictedSampleReport_v8.sh <BIOProjectID> to refresh the report and study summary.
#
# Written 30 Sep 2026 by WMD
# Usage: libprep_recall.sh <BIOProjectID> [path/to/libprep_call.sh]

###############################################################################
# USER CONFIGURATION
# Update the variables below for your system before running.
###############################################################################
# Your home/base directory (projects live in $myDir/data/expression/GEO/<BIOProjectID>)
myDir="/path/to/your/home"
###############################################################################

if [ $# -lt 1 ]; then
    echo "Usage: $0 <BIOProjectID> [libprep_call.sh]"
    exit 1
fi
BIOProjectID=$1
LIBPREP_CALL=${2:-$(dirname "$(readlink -f "$0")")/libprep_call.sh}
reads_dir=${myDir}/data/expression/GEO/${BIOProjectID}/reads_fastq

if [ ! -f "$LIBPREP_CALL" ]; then
    echo "ERROR: libprep_call.sh not found: $LIBPREP_CALL"
    exit 1
fi
LIBPREP_CALL=$(readlink -f "$LIBPREP_CALL")   # a bare name like libprep_call.sh is not on PATH
if [ ! -d "$reads_dir" ]; then
    echo "ERROR: no reads_fastq directory: $reads_dir"
    exit 1
fi

n=0; failed=()
for sdir in "$reads_dir"/*/; do
    sample=$(basename "$sdir")
    [ -d "$sdir/log_files/STAR" ] || continue
    n=$((n + 1))
    bash "$LIBPREP_CALL" "$sample" "$sdir" || failed+=("$sample")
done

echo "$BIOProjectID: $n samples, $(( n - ${#failed[@]} )) called, ${#failed[@]} failed"
if [ ${#failed[@]} -gt 0 ]; then
    echo "Failed: ${failed[*]}"
    exit 1
fi
