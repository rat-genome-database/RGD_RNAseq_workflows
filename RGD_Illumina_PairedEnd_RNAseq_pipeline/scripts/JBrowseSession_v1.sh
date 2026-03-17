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

# Script directory — update if you move the pipeline scripts

# Load Python module for this job
module load python/3.12.10

echo "Starting JBrowse session build for BIOProjectID=${BIOProjectID}"
echo "PRJdir=${PRJdir}"
echo "baseDir=${baseDir}"

#python "${baseDir}/make_jbrowse_session_for_bioproject.py" \
python "${SCRIPT_DIR}/make_jbrowse_session_for_bioproject.py" \
  "${BIOProjectID}" \
  "${PRJdir}" \
  "${baseDir}"

echo "JBrowse session build completed."
