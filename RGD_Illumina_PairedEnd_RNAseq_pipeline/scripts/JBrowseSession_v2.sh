#!/usr/bin/env bash
#SBATCH --job-name=JBrowseSession
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --account=your-slurm-account

###############################################################################
# USER CONFIGURATION
# Update the variables below for your system before running.
###############################################################################
# Directory containing all pipeline scripts
SCRIPT_DIR="/path/to/your/pipeline/scripts"
###############################################################################


set -euo pipefail

# v2 (3 June 2026): build the session from the STARQC PASS accession list so that
# QC-failed samples are never included. The PASS AccList is the same single
# source of truth the orchestrator uses to drive RSEM / ComputeSex / BWjson:
#   ${Logdir}/STARQC/${BIOProjectID}_Unique_AccList_PASS.txt
# It is passed to make_jbrowse_session_for_bioproject.py as the 4th argument.
# Requires BIOProjectID, PRJdir, baseDir and Logdir exported by the calling script.

# Load Python module for this job
module load python/3.12.10

# PASS-only accession list produced by the STARQC filtering step
passAccList="${Logdir}/STARQC/${BIOProjectID}_Unique_AccList_PASS.txt"

echo "Starting JBrowse session build for BIOProjectID=${BIOProjectID}"
echo "PRJdir=${PRJdir}"
echo "baseDir=${baseDir}"
echo "passAccList=${passAccList}"

if [ ! -s "${passAccList}" ]; then
  echo "ERROR: PASS accession list not found or empty: ${passAccList}" >&2
  echo "Cannot build JBrowse session without the STARQC PASS list." >&2
  exit 1
fi

python "${SCRIPT_DIR}/make_jbrowse_session_for_bioproject.py" \
  "${BIOProjectID}" \
  "${PRJdir}" \
  "${baseDir}" \
  "${passAccList}"

echo "JBrowse session build completed."
