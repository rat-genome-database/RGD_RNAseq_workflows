#!/bin/bash
#SBATCH --job-name=RSEMmatrix
#SBATCH --ntasks=1
#SBATCH --mem=1gb
#SBATCH --time=00:05:00
#SBATCH --account=akwitek
#SBATCH --partition=normal
#SBATCH --output=%x-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wdemos@mcw.edu
#SBATCH --error=%x-%j.err

# Written by Wendy Demos
# Last modified on 1 Nov 2024 by WMD
# Enable debugging by printing each command before execution
set -x
# Load necessary modules
module load rsem/1.3.3 samtools/1.20 

# Retrieve environment variables
AccList=${AccList}
PRJdir=${PRJdir}
WrkDir="/home/wdemos/data/expression/GEO/${BIOProjectID}"
BIOProjectID=${BIOProjectID}
echo "Project directory is set to $PRJdir"
scratch_dir="/scratch/g/akwitek/wdemos/$BIOProjectID"
echo "Scratch working directory is set to $scratch_dir"
baseDir=${baseDir} #"${myDir}/data/expression/GEO/${BIOProjectID}"
Logdir="${baseDir}/log_files"
echo "Logdir is set to: $Logdir"
mkdir -p "$scratch_dir"

# Ensure project directory exists
if [ ! -d "$PRJdir" ]; then
    echo "Creating project directory: $PRJdir"
    mkdir -p "$WrkDir" || { echo "Error: Failed to create directory $WrkDir"; exit 1; }
else
    echo "Project directory already exists: $PRJdir"
fi

# Redirect output to the log file in the project directory
log_file="${Logdir}/${BIOProjectID}_RSEMmatrix.out"
exec 1>"$log_file" 2>&1

# Capture start time
start_time=$(date +%s)

# Move into the scratch directory
echo "Moving into the working directory ${scratch_dir}"
cd "$scratch_dir" || { echo "Error: Failed to move into directory $scratch_dir"; exit 1; }

# Arrays to store results file paths
# Arrays to store results file paths and track unique geo_accessions
transcripts_files=()
genes_files=()
declare -A seen_geo_accession

# Read the AccList file to process each unique geo_accession
echo "Reading accession list: $AccList"
while IFS=$'\t' read -r Run geo_accession Tissue Strain Sex PMID GEOpath Title Sample_characteristics StrainInfo; do
    # Skip if this geo_accession has already been processed
    if [[ -n "${seen_geo_accession[$geo_accession]}" ]]; then
        continue
    fi

    # Mark this geo_accession as processed
    seen_geo_accession[$geo_accession]=1

    # Locate isoform and gene result files
    transcripts_file="${scratch_dir}/${geo_accession}/${geo_accession}.transcripts.results"
    genes_file="${scratch_dir}/${geo_accession}/${geo_accession}.genes.results"
    
    if [ -e "$transcripts_file" ]; then
        transcripts_files+=("$transcripts_file")
    else
        echo "Warning: transcripts results file not found for Run $Run"
    fi

    if [ -e "$genes_file" ]; then
        genes_files+=("$genes_file")
    else
        echo "Warning: Gene results file not found for Run $Run"
    fi
done < "$AccList"


# Check if any transcripts or genes files were found
if [ ${#transcripts_files[@]} -eq 0 ]; then
    echo "Error: No transcripts results files found for any samples."
    exit 1
else
    echo "Combining transcripts results files: ${transcripts_files[@]}"

    /home/wdemos/dependencies/rsem-generate-data-matrix "${transcripts_files[@]}" > "${WrkDir}/${BIOProjectID}.transcripts.TPM.matrix"
fi

if [ ${#genes_files[@]} -eq 0 ]; then
    echo "Error: No gene results files found for any samples."
    exit 1
else
    echo "Combining gene results files: ${genes_files[@]}"
    /home/wdemos/dependencies/rsem-generate-data-matrix "${genes_files[@]}" > "${WrkDir}/${BIOProjectID}.genes.TPM.matrix"
fi



#######################################
# Gather Data QC metrics with MultiQC #
#######################################

echo "Gathering QC metrics with MultiQC"
module load python/3.9.1 multiqc/1.18
python3 -m multiqc "$scratch_dir" -o "${WrkDir}" -n "${BIOProjectID}_final_multiQC_repor"t || { echo "MultiQC failed"; exit 1; }

##################################
# Run Sample Sex Conflict Script #
##################################

sh /home/wdemos/STAR_RSEM_pipeline/illuminaPaired_multipleSRA_IDperSample/ConflictedSampleReport.sh ${BIOProjectID}

#####################################################
# copy Acc List and fastq list to project directory #
#####################################################


cp "$AccList" "${WrkDir}"
cp -R ${scratch_dir}/logs "${WrkDir}"
cp ${scratch_dir}/run_combination_log.txt "${WrkDir}"


###########################################
# Capture end time and print elapsed time #
###########################################

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
echo "Total elapsed time: $(date -d@$elapsed_time -u +%H:%M:%S)"

