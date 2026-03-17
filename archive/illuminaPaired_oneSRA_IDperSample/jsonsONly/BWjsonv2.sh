#!/bin/sh
#SBATCH --job-name=bwjson
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5000kb
#SBATCH --time=00:01:00
#SBATCH --account=akwitek
#SBATCH --partition=normal
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wdemos@mcw.edu

# Enable strict mode and command tracing
set -euxo pipefail

# Retrieve environment variables
Run="${Run}"
geo_accession="${geo_accession}"
BIOProjectID="${BIOProjectID}"
Tissue="${tissue}"
Strain="${strain}"
StrainInfo="${StrainInfo}"
Sex="${sex}"
Title="${title}"
Sample_characteristics="${Sample_characteristics}"
PMID="${PMID}"
GEOpath="${GEOpath}"
unique_name="${unique_name}"
PRJdir="${PRJdir}"
scratch_dir="${scratch_dir}"
Logdir="${Logdir}"

# Set log file and redirect stdout/stderr
log_file="${PRJdir}/${geo_accession}/${geo_accession}_bwJSON.out"
mkdir -p "${PRJdir}/${geo_accession}"
exec 1>"$log_file" 2>&1

# Timestamp start
echo "[$(date)] Starting BW JSON creation for $geo_accession"

# Move into scratch directory
cd "$scratch_dir"
echo "Entered directory: $scratch_dir"

# Path to sex result file
sex_results_file="${Logdir}/${BIOProjectID}_sex_result.txt"

# Lookup computed sex
CalcSex=$(awk -v sample="$geo_accession" -F'\t' '$1 == sample {print $3}' "$sex_results_file" || true)
CalcSex="${CalcSex:-Unknown}"

# Path to bigWig file
bw_file="${PRJdir}/${geo_accession}/RNAseq_${unique_name}.bigwig"
echo "Checking for BigWig file: $bw_file"

if [ -f "$bw_file" ]; then
    echo "BigWig found. Generating JSON..."

    json_file="RNAseq_${unique_name}.json"
    cat <<EOF > "${PRJdir}/${geo_accession}/$json_file"
{
  "type": "FeatureTrack",
  "trackId": "RNAseq_${unique_name}",
  "name": "RNAseq_${unique_name}",
  "assemblyNames": ["mRatBN7.2"],
  "metadata": {
    "Project Title": "${Title}",
    "GEO Project Accession ID": "${BIOProjectID}",
    "Tissue": "${Tissue}",
    "Strain": "${Strain}, ${StrainInfo}",
    "Sex": "${Sex}",
    "Computed Sex": "${CalcSex}",
    "GEO Sample Accession ID": "${geo_accession}",
    "Sample Characteristic": "${Sample_characteristics}",
    "PMID": "${PMID}",
    "Project Link": "${GEOpath}",
    "Data Processing": "HPC RGD workflow",
    "Read alignment": "STAR version=2.7.10b_alpha_220111",
    "Genome version": "mRatBN7.2",
    "Expression Quantification": "RSEM v1.3.1"
  },
  "adapter": {
    "type": "BigWigAdapter",
    "bigWigLocation": {
      "locationType": "UriLocation",
      "uri": "json/RNAseq_${unique_name}.bigwig"
    }
  },
  "displays": [
    {
      "type": "LinearWiggleDisplay",
      "displayId": "RNAseq_${unique_name}-LinearWiggleDisplay"
    }
  ]
}
EOF

    if [ -f "${PRJdir}/${geo_accession}/$json_file" ]; then
        echo "[$(date)] JSON created: $json_file"
    else
        echo "[$(date)] Error: Failed to create JSON file."
    fi
else
    echo "[$(date)] Error: BigWig not found. Skipping JSON creation."
fi

# Timestamp end
echo "[$(date)] Finished BW JSON script for $geo_accession"
