#!/bin/sh
#SBATCH --job-name=genejson
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5000kb
#SBATCH --time=00:01:00
#SBATCH --account=akwitek
#SBATCH --partition=normal
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wdemos@mcw.edu
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

# Strict mode for debugging
set -euxo pipefail

# ---------------------------------------------------------
# Load environment variables passed by sbatch --export
# ---------------------------------------------------------
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

# ---------------------------------------------------------
# Logging
# ---------------------------------------------------------
log_file="${PRJdir}/${geo_accession}/${geo_accession}_geneBedJSON.out"
mkdir -p "${PRJdir}/${geo_accession}"
exec 1>"$log_file" 2>&1

echo "[$(date)] gJSON script started for $geo_accession"

# ---------------------------------------------------------
# Change to working directory
# ---------------------------------------------------------
cd "$scratch_dir"
echo "Entered scratch directory: $scratch_dir"

# ---------------------------------------------------------
# Get computed sex
# ---------------------------------------------------------
sex_results_file="${Logdir}/${BIOProjectID}_sex_result.txt"
CalcSex=$(awk -v sample="$geo_accession" -F'\t' '$1 == sample {print $3}' "$sex_results_file" || true)
CalcSex="${CalcSex:-Unknown}"

# ---------------------------------------------------------
# Generate JSON file
# ---------------------------------------------------------
json_file="RNAseq_${unique_name}.geneTPMbed.json"
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
    "Read alignment": "STAR v2.7.10b_alpha_220111",
    "Genome version": "mRatBN7.2",
    "Expression Quantification": "RSEM v1.3.1",
    "Below Cutoff": ">0 TPM <0.5 TPM (gray)",
    "Low": ">= 0.5 TPM <10 TPM (aqua)",
    "Medium": ">=10 TPM <1000 TPM (sky blue)",
    "High": ">=1000 TPM (dark blue)"
  },
  "adapter": {
    "type": "BedTabixAdapter",
    "bedGzLocation": {
      "locationType": "UriLocation",
      "uri": "/jbrowse2/RNAseq_${unique_name}_geneTPMfinalOP.bed.gz"
    },
    "index": {
      "location": {
        "locationType": "UriLocation",
        "uri": "/jbrowse2/RNAseq_${unique_name}_geneTPMfinalOP.bed.gz.tbi"
      }
    }
  },
  "displays": [
    {
      "type": "LinearBasicDisplay",
      "displayId": "geneTPM-LinearBasicDisplay"
    },
    {
      "type": "LinearArcDisplay",
      "displayId": "geneTPM-LinearArcDisplay"
    }
  ]
}
EOF

# ---------------------------------------------------------
# Verify that JSON file was created
# ---------------------------------------------------------
if [ -f "${PRJdir}/${geo_accession}/$json_file" ]; then
    echo "[$(date)] JSON created successfully: $json_file"
else
    echo "[$(date)] ERROR: Failed to create JSON file."
    exit 1
fi

echo "[$(date)] gJSON script completed for $geo_accession"
