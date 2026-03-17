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

## Debugging: Enable command tracing
set -x

# Retrieve environment variables
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

echo "Logdir is set to: $Logdir"
echo "Scratch directory: $scratch_dir"
echo "bigwig Run variable is ${Run}"

# Move into the directory
cd "$scratch_dir" || { echo "Error: Failed to move into directory $scratch_dir"; exit 1; }
echo "Moved into directory: $scratch_dir"

# Redirect output to the log file in the project directory after it is created
#log_file="${Logdir}/${geo_accession}/${geo_accession}_bwJSON.out"
#mkdir -p "${Logdir}/${geo_accession}"  # Ensure the directory exists
#exec 1>"$log_file" 2>&1
# Redirect output to the log file in the project directory after it is created
log_file="${PRJdir}/${geo_accession}/${geo_accession}_bwJSON.out"
mkdir -p "${PRJdir}/${geo_accession}"  # Ensure the directory exists
exec 1>"$log_file" 2>&1

# Capture start time
start_time=$(date +%s)

# Path to the sex results file
#sex_results_file="${Logdir}/${BIOProjectID}_sex_result.txt"
sex_results_file="/home/wdemos/data/expression/GEO/${BIOProjectID}/${BIOProjectID}_sex_result.txt"

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

# Path to the BigWig file
#bw_file="${Logdir}/RNAseq_${unique_name}.bigwig"
bw_file="${PRJdir}/${geo_accession}/RNAseq_${unique_name}.bigwig"
echo "Checking if BigWig file exists at: $bw_file"

# Check if the BigWig file exists
if [ -f "$bw_file" ]; then
    echo "$bw_file exists. Proceeding to create JSON."

    # Create JSON content with variables substituted
    json_content=$(cat <<EOF
{
  "type": "FeatureTrack",
  "trackId": "RNAseq_${unique_name}",
  "name": "RNAseq_${unique_name}",
  "category": [
        "RNA-Seq",
        "${Tissue}",
        "${Strain}"
      ],
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
)
    
    # Write the JSON content to a file
    json_file="RNAseq_${unique_name}.json"
    echo "$json_content" > "${PRJdir}/${geo_accession}/$json_file"

    # Verify if the JSON file was created
    if [ -f "${PRJdir}/${geo_accession}/$json_file" ]; then
        echo "JSON file created successfully: $json_file"
    else
        echo "Error: Failed to create JSON file: $json_file"
    fi

else
    echo "Error: $bw_file does not exist. Skipping JSON creation for this sample."
fi

# Capture end time and calculate wall clock time
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
echo "Processing completed in $elapsed_time seconds."

