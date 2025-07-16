#!/bin/sh
#SBATCH --job-name=Genejson
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5000kb
#SBATCH --time=00:01:00
#SBATCH --account=akwitek
#SBATCH --partition=normal
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wdemos@mcw.edu
set -x
##written 13 September 2024 by WMD
##CHANGE LOG
## 7 Nov 2024
## updated json format and added catch to  check json was generated
## updated to run on RCC cluster
## 18 June removed form header: #SBATCH --output=%x-%j.out, #SBATCH --error=Genejson.err
# retrieve environment variables
Run=${Run}
geo_accession=${geo_accession}
BIOProjectID="$BIOProjectID" # This should be the GEO Project Name
Tissue="$tissue"
Strain="$strain"
StrainInfo="$StrainInfo"
Sex="$sex"
Title="$title"
Sample_characteristics="$Sample_characteristics"
PMID="$PMID"
GEOpath="$GEOpath"
unique_name=${unique_name}
PRJdir=${PRJdir} #$myDir/data/expression/GEO/$BIOProjectID/reads_fastq
scratch_dir=${scratch_dir} #/scratch/g/akwitek/wdemos/$BIOProjectID
Logdir=${Logdir} #$myDir/data/expression/GEO/$BIOProjectID
#Logdir="/home/wdemos/data/expression/GEO/$BIOProjectID"

# Final logging
echo "Logdir is set to: $Logdir"
echo "Scratch directory: $scratch_dir"
echo "gene bed json generations 'Run' variable is ${Run}"

# Move into the directory
cd "$scratch_dir" || { echo "Error: Failed to move into directory $scratch_dir"; exit 1; }
echo "Moved into directory: $scratch_dir"

# Redirect output to the log file in the project directory after it is created
log_file="${PRJdir}/${geo_accession}/${geo_accession}_geneBedJSON.out"
exec 1>"$log_file" 2>&1

# Capture start time
start_time=$(date +%s)

# Path to the sex results file
sex_results_file="${Logdir}/${BIOProjectID}_sex_result.txt"

# Create a function to get CalcSex for a given sample
get_calc_sex() {
    local sample_id=$1
    awk -v sample="$sample_id" -F'\t' '$1 == sample {print $3}' "$sex_results_file"
}

# Define sample ID
sample="$geo_accession"
echo "Processing sample for project: $Run"

# Get CalcSex for the sample
CalcSex=$(get_calc_sex "$sample")
if [ -z "$CalcSex" ]; then
    CalcSex="Unknown"
fi

genebed_file="${PRJdir}/${sample}/RNAseq_${unique_name}_geneTPMfinalOP.bed.gz"
echo "Generating JSON for $genebed_file."

    # Create JSON content with variables substituted
    json_content=$(cat <<EOF
{
  "type": "FeatureTrack",
  "trackId": "RNAseq_${unique_name}",
  "name": "RNAseq_${unique_name}",
  "assemblyNames": [
    "mRatBN7.2"
  ],
  "metadata": {
    "Project Title": "${Title}",
    "GEO Project Accession ID": "${BIOProjectID}",
    "Tissue": "${Tissue}",
    "Strain": "${Strain}, ${StrainInfo}",
    "Sex": "${Sex}",
    "Computed Sex": "${CalcSex}",
    "GEO Sample Accession ID": "${sample}",
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
      "uri":  "/jbrowse2/RNAseq_${unique_name}_geneTPMfinalOP.bed.gz"
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
      "displayId": "SetRGBColor_filt${numID}",
      "renderer": {
        "type": "SvgFeatureRenderer",
        "color1": "jexl:get(feature, 'field8')=='128,128,128' ? '#808080' : '#000000'\n&&get(feature, 'field8')=='173,216,230' ? '#99FFFF' : '#000000'\n&&get(feature, 'field8')=='0,0,205' ? '#3399FF' : '#000000'\n&&get(feature, 'field8')=='0,0,139' ? '#0000CC' : '#000000'"
      }
    },
    {
      "type": "LinearArcDisplay",
      "displayId": "geneTPM-LinearArcDisplay"
    }
  ]
}
EOF
)
    # Write the JSON content to a file
    json_file="RNAseq_${unique_name}.geneTPMbed.json"
    echo "$json_content" > "${PRJdir}/${geo_accession}/$json_file"

    # Verify if the JSON file was created
    if [ -f "$json_file" ]; then
        echo "JSON file created successfully: $json_file"
    else
        echo "Error: Failed to create JSON file: $json_file"
    fi
