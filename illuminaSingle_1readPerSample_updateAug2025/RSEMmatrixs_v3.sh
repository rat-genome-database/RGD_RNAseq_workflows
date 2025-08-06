#!/bin/bash
#SBATCH --job-name=RmtxS
#SBATCH --ntasks=1
#SBATCH --mem=1gb
#SBATCH --time=00:05:00
#SBATCH --account=akwitek
#SBATCH --partition=normal
#SBATCH --output=%x-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wdemos@mcw.edu
#SBATCH --error=RmtxS.err

# Written by Wendy Demos
# Last modified on 1 Nov 2024 by WMD
# Enable debugging by printing each command before execution
set -x
# Load necessary modules
module load rsem/1.3.3 samtools/1.20 

# Retrieve environment variables
AccList=${AccList}
PRJdir=${PRJdir}
BIOProjectID=${BIOProjectID}
echo "Project directory is set to $PRJdir"
scratch_dir="/scratch/g/akwitek/wdemos/$BIOProjectID"
echo "Scratch working directory is set to $scratch_dir"
baseDir=${baseDir}
echo "Logdir is set to: $Logdir"
mkdir -p "$scratch_dir"

# Ensure project directory exists
if [ ! -d "$PRJdir" ]; then
    echo "Creating project directory: $PRJdir"
    mkdir -p "$PRJdir" || { echo "Error: Failed to create directory $PRJdir"; exit 1; }
else
    echo "Project directory already exists: $PRJdir"
fi

# Redirect output to the log file in the project directory
log_file="${Logdir}/${BIOProjectID}_RmtxS.out"
exec 1>"$log_file" 2>&1

# Capture start time
start_time=$(date +%s)

# Move into the scratch directory
echo "Moving into the working directory ${scratch_dir}"
cd "$scratch_dir" || { echo "Error: Failed to move into directory $scratch_dir"; exit 1; }

# Arrays to store results file paths
transcripts_files=()
genes_files=()

# Read the AccList file to process each sample
echo "Reading accession list: $AccList"
while IFS=$'\t' read -r Run geo_accession Tissue Strain Sex PMID GEOpath Title Sample_characteristics StrainInfo; do
    # Locate isoform and gene result files for each sample
    transcripts_file="${scratch_dir}/${Run}/${geo_accession}.transcripts.results"
    genes_file="${scratch_dir}/${Run}/${geo_accession}.genes.results"
    
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
    #rsem-generate-data-matrix "${transcripts_files[@]}" > "${scratch_dir}/${BIOProjectID}.transcripts.TPM.matrix"
    /home/wdemos/dependencies/rsem-generate-data-matrix "${transcripts_files[@]}" > "${scratch_dir}/${BIOProjectID}.transcripts.TPM.matrix"
fi

if [ ${#genes_files[@]} -eq 0 ]; then
    echo "Error: No gene results files found for any samples."
    exit 1
else
    echo "Combining gene results files: ${genes_files[@]}"
#    rsem-generate-data-matrix "${genes_files[@]}" > "${scratch_dir}/${BIOProjectID}.genes.TPM.matrix"
    /home/wdemos/dependencies/rsem-generate-data-matrix "${genes_files[@]}" > "${scratch_dir}/${BIOProjectID}.genes.TPM.matrix"

fi
#Copy output files to static GEO directory
echo "copying transcript matric to static GEO project directory:$baseDir "
cp "${scratch_dir}/${BIOProjectID}.transcripts.TPM.matrix" ${baseDir}
cp "${scratch_dir}/${BIOProjectID}.genes.TPM.matrix" ${baseDir}

#######################################
# Gather Data QC metrics with MultiQC #
#######################################

echo "Gathering QC metrics with MultiQC"
module load python/3.9.1 multiqc/1.18
#python3 -m multiqc "$scratch_dir" -o "${Logdir}" || { echo "MultiQC failed"; exit 1; }
python3 -m multiqc "$scratch_dir" -o "${baseDir}" -n ${BIOProjectID}_final_multiQC_report || { echo "MultiQC failed"; exit 1; }

##################################
# Run Sample Sex Conflict Script #
##################################
sh /home/wdemos/STAR_RSEM_pipeline/illuminaSingle_1readPerSample/ConflictedSampleReport_v3.sh ${BIOProjectID}

############################################
# copy Accession List to project directory #
############################################

cp /home/wdemos/STAR_RSEM_pipeline/AccessionLists/${BIOProjectID}_AccList.txt $baseDir

###########################################
# Capture end time and print elapsed time #
###########################################

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
echo "Total elapsed time: $(date -d@$elapsed_time -u +%H:%M:%S)"
module unload rsem/1.3.3 samtools/1.20 
