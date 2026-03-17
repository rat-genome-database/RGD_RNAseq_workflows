#!/usr/bin/env bash
#SBATCH --job-name=JBrowseSession
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --account=akwitek

set -euo pipefail

# Load Python module for this job
module load python/3.12.10

echo "Starting JBrowse session build for BIOProjectID=${BIOProjectID}"
echo "PRJdir=${PRJdir}"
echo "baseDir=${baseDir}"

#python "${baseDir}/make_jbrowse_session_for_bioproject.py" \
python "/home/wdemos/STAR_RSEM_pipeline/GRC8_bulk_illuminaPaired_oneSRA/make_jbrowse_session_for_bioproject.py" \
  "${BIOProjectID}" \
  "${PRJdir}" \
  "${baseDir}"

echo "JBrowse session build completed."
