#!/bin/sh
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1mb
#SBATCH --time=00:10:00
#SBATCH --account=akwitek
#SBATCH --partition=normal
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wdemos@mcw.edu

# Debugging mode
set -euxo pipefail

# ---------------------------------------------------------
# Retrieve variables passed via sbatch --export
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
# Setup logging
# ---------------------------------------------------------
log_file="${PRJdir}/${geo_accession}/${geo_accession}_TxBedJSON.out"
mkdir -p "${PRJdir}/${geo_accession}"
exec 1>"$log_file" 2>&1

echo "[$(date)] tJSON script started for $geo_accession"

# ---------------------------------------------------------
# Change to scratch working directory
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
# Write transcript-level JSON
# ---------------------------------------------------------
TXbed_file="${PRJdir}/${geo_accession}/RNAseq_${unique_name}_TXTPMfinalOP.bed.gz"
echo "Generating transcript JSON for: $TXbed_file"

json_file="RNAseq_${unique_name}.TXTPMbed.json"
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
    "Expression Quantification": "RSEM v1.3.1",
    "Below Cutoff Threshold and bar color": ">0 TPM <0.5 TPM , gray",
    "Low Threshold and bar color": ">=0.5 TPM  <10 TPM , aqua",
    "Medium Threshold and bar color": ">= 10 TPM  <1000 TPM , sky blue",
    "High Threshold and bar color": ">=1000 TPM, dark blue"
  },
  "adapter": {
    "type": "BedTabixAdapter",
    "bedGzLocation": {
      "locationType": "UriLocation",
      "uri": "/jbrowse2/RNAseq_${unique_name}_TXTPMfinalOP.bed.gz"
    },
    "index": {
      "location": {
        "locationType": "UriLocation",
        "uri": "/jbrowse2/RNAseq_${unique_name}TXTPMfinalOP.bed.gz.tbi"
      }
    }
  },
  "displays": [
    {
      "type": "LinearBasicDisplay",
      "displayId": "txTPM-LinearBasicDisplay",
      "renderer": {
        "type": "SvgFeatureRenderer",
        "color1": "jexl:get(feature, 'field8')=='128,128,128' ? '#808080' : '#000000'\n&&get(feature, 'field8')=='173,216,230' ? '#99FFFF' : '#000000'\n&&get(feature, 'field8')=='0,0,205' ? '#3399FF' : '#000000'\n&&get(feature, 'field8')=='0,0,139' ? '#0000CC' : '#000000'"
      }
    },
    {
      "type": "LinearArcDisplay",
      "displayId": "TXTPM-LinearArcDisplay"
    }
  ]
}
EOF

# ---------------------------------------------------------
# Verify JSON creation
# ---------------------------------------------------------
if [ -f "${PRJdir}/${geo_accession}/$json_file" ]; then
    echo "[$(date)] Transcript JSON created: $json_file"
else
    echo "[$(date)] ERROR: Transcript JSON not created"
    exit 1
fi

echo "[$(date)] tJSON script completed for $geo_accession"
