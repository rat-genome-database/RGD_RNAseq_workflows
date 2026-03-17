#!/bin/bash
#SBATCH --job-name=S2Rs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=3gb
#SBATCH --time=06:00:00
#SBATCH --account=akwitek
#SBATCH --output=%x-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wdemos@mcw.edu
#SBATCH --error=S2Rs.err

# Enable debugging by printing each command before execution
set -x

# Check command line arguments
if [ $# -ne 3 ]; then
    echo "Usage: $0 <SRA_accession_list.txt> <BioProject_ID> <Read_Length>"
    exit 1
fi

AccList=$1 
BIOProjectID=$2
length=$3
myDir="/home/wdemos"
PRJdir="$myDir/data/expression/GEO/$BIOProjectID/reads_fastq"
Logdir="$myDir/data/expression/GEO/$BIOProjectID"
scratch_dir="/scratch/g/akwitek/wdemos/$BIOProjectID"
INDEX_DIR="$scratch_dir/RefIndex"
REFdir="/home/wdemos/NCBI_1August2023"


# Create directories
mkdir -p "$PRJdir"
mkdir -p "$Logdir"

# Redirect output to log file
log_file="${Logdir}/${BIOProjectID}_Step2.out"
exec 1>"$log_file" 2>&1

# Function to wait for jobs to complete
wait_for_jobs() {
    for job_id in "${job_ids[@]}"; do
        while squeue -u $USER | grep -q "$job_id"; do
            sleep 5
        done
    done
}

##############################
# STAR Reference Preparation #
##############################
if [ ! -d "$INDEX_DIR" ] || [ ! -f "$INDEX_DIR/SA" ] || [ ! -f "$INDEX_DIR/SAindex" ] || [ ! -f "$INDEX_DIR/Genome" ]; then
    echo "Submitting starRef.sh to generate the STAR reference."
    starRef_job_id=$(sbatch --export=PRJdir="$PRJdir",Logdir="$Logdir",BIOProjectID="$BIOProjectID",length="$length" \
                        --output="$Logdir/starRef-%j.out" \
                        --error="$Logdir/starRef-%j.err" starRef.sh | awk '{print $4}')
    if [ -z "$starRef_job_id" ]; then
        echo "Error: starRef.sh submission failed."
        exit 1
    fi
    echo "STAR reference job ID: $starRef_job_id"
else
    echo "STAR reference already exists. Skipping starRef.sh."
    starRef_job_id=""
fi

####################################
# STAR Alignment Jobs for Samples #
####################################
job_ids=()  # Array to store job IDs for STAR alignment
while IFS=$'\t' read -r run geo_accession Tissue Strain Sex PMID GEOpath Title Sample_characteristics StrainInfo; do
    # Skip the header row
    if [ "$run" == "Run" ]; then
        continue
    fi

    unique_name="${Tissue}_${Strain}_${Sex}_${geo_accession}"
    sample_output_dir="$PRJdir/$geo_accession"
    mkdir -p "$sample_output_dir"  # Ensure directory exists

    # Debugging: Print the paths and variables
    echo "Submitting STAR job for $geo_accession"

    if [ -z "$starRef_job_id" ]; then
        job_id=$(sbatch --job-name="${geo_accession}_STAR72" \
                        --export=BIOProjectID="$BIOProjectID",Run="$run",geo_accession="$geo_accession",PRJdir="$PRJdir",length="$length",unique_name="$unique_name",Logdir="$Logdir" \
                        --output="$sample_output_dir/STAR-%j.out" \
                        --error="$sample_output_dir/STAR-%j.err" \
                        stars.sh | awk '{print $4}')
    else
        job_id=$(sbatch --dependency=afterok:$starRef_job_id \
                        --job-name="${geo_accession}_STAR72" \
                        --export=BIOProjectID="$BIOProjectID",Run="$run",geo_accession="$geo_accession",PRJdir="$PRJdir",length="$length",unique_name="$unique_name",Logdir="$Logdir" \
                        --output="$sample_output_dir/STAR-%j.out" \
                        --error="$sample_output_dir/STAR-%j.err" \
                        stars.sh | awk '{print $4}')
    fi

    # Check if sbatch failed
    if [ $? -ne 0 ]; then
        echo "Error submitting STAR job for $geo_accession"
        continue
    fi
    job_ids+=("$job_id")
done < "$AccList"

wait_for_jobs

##############################
# Sample Sex Estimation #
##############################
echo "Alignment with STAR is complete for all samples. Computing sample sex."
job_id=$(sbatch --export=PRJdir="$PRJdir",Logdir="$Logdir",AccList="$AccList",BIOProjectID="$BIOProjectID" \
              --output="$Logdir/ComputeSex-%j.out" \
              --error="$Logdir/ComputeSex-%j.err" ComputeSex.sh "$AccList" "$BIOProjectID" | awk '{print $4}')
if [ $? -ne 0 ]; then
    echo "Error running ComputeSex.sh"
    exit 1
fi

#wait_for_jobs

##############################
# RSEM Reference Preparation #
##############################
if [ -f "$scratch_dir/rsemref.chrlist" ] && [ -f "$scratch_dir/rsemref.grp" ] && [ -f "$scratch_dir/rsemref.idx.fa" ] && [ -f "$scratch_dir/rsemref.n2g.idx.fa" ] && [ -f "$scratch_dir/rsemref.seq" ] && [ -f "$scratch_dir/rsemref.ti" ] && [ -f "$scratch_dir/rsemref.transcripts.fa" ]; then
    echo "RSEM reference already exists. Skipping RSEMref.sh."
else
    echo "Submitting RSEMref.sh to generate the RSEM reference."
    rsemRef_job_id=$(sbatch --export=BIOProjectID="$BIOProjectID" \
                           --output="$Logdir/RSEMref-%j.out" \
                           --error="$Logdir/RSEMref-%j.err" RSEMref.sh | awk '{print $4}')
    if [ -z "$rsemRef_job_id" ]; then
        echo "Error: RSEMref.sh submission failed."
        exit 1
    fi
    echo "RSEM reference job ID: $rsemRef_job_id"
fi

###########################
# RSEM Jobs for Samples #
###########################
job_ids=()  # Array to store job IDs for RSEM
while IFS=$'\t' read -r run geo_accession Tissue Strain Sex PMID GEOpath Title Sample_characteristics StrainInfo; do
    # Skip the header row
    if [ "$run" == "Run" ]; then
        continue
    fi

    sample_output_dir="$PRJdir/$geo_accession"
    mkdir -p "$sample_output_dir"  # Ensure directory exists

    # Debugging: Print the paths and variables
    echo "Submitting RSEM job for $geo_accession"

    if [ -z "$rsemRef_job_id" ]; then
        job_id=$(sbatch --job-name="${geo_accession}_RSEM72" \
                        --export=BIOProjectID="$BIOProjectID",Run="$run",geo_accession="$geo_accession",PRJdir="$PRJdir",Logdir="$Logdir" \
                        --output="$sample_output_dir/RSEM-%j.out" \
                       --error="$sample_output_dir/RSEM-%j.err" \
                        RSEMs.sh | awk '{print $4}')
    else
        job_id=$(sbatch --dependency=afterok:$rsemRef_job_id \
                        --job-name="${geo_accession}_RSEM72" \
                        --export=BIOProjectID="$BIOProjectID",Run="$run",geo_accession="$geo_accession",PRJdir="$PRJdir",Logdir="$Logdir" \
                        --output="$sample_output_dir/RSEM-%j.out" \
                        --error="$sample_output_dir/RSEM-%j.err" \
                        RSEMs.sh | awk '{print $4}')
    fi

    # Check if sbatch failed
    if [ $? -ne 0 ]; then
        echo "Error submitting RSEM job for $geo_accession"
        continue
    fi
    job_ids+=("$job_id")
done < "$AccList"

wait_for_jobs

#####################################
# Final Steps - Matrix and QC files #
#####################################
echo "Generating combined matrices for all samples and final MultiQC report."

# Submit the RSEMmatrix.sh job and capture output
job_id_output=$(sbatch --export=Run="$run",geo_accession="$geo_accession",PRJdir="$PRJdir",Logdir="$Logdir",AccList="$AccList",BIOProjectID="$BIOProjectID" \
              --output="$Logdir/RSEMmatrix-%j.out" \
              --error="$Logdir/RSEMmatrix-%j.err" RSEMmatrixs.sh 2>&1)

# Extract job ID
matrix_job_id=$(echo "$job_id_output" | awk '{print $4}')

# Log the submission output for debugging
echo "RSEMmatrixs.sh submission output: $job_id_output" >> "$log_file"

# Validate that a proper job ID was obtained
if [[ -z "$matrix_job_id" || ! "$matrix_job_id" =~ ^[0-9]+$ ]]; then
    echo "Error: Invalid matrix_job_id: '$matrix_job_id'" >> "$log_file"
    exit 1
fi

echo "RSEMmatrixs.sh submitted with job ID: $matrix_job_id" >> "$log_file"

# Wait for jobs before proceeding
echo "Waiting for RSEMmatrix job to complete..." >> "$log_file"
wait_for_jobs
echo "RSEMmatrix job completed. Proceeding with JSON and BED generation..." >> "$log_file"

###################################
# Generate BED and JSON files
###################################

declare -A json_bed_job_ids  # Track job IDs for JSON and BED files
declare -A seen_samples  # Track processed geo_accession values

while IFS=$'\t' read -r run geo_accession Tissue Strain Sex PMID GEOpath Title Sample_characteristics StrainInfo; do
    # Skip the header row
    if [ "$run" == "Run" ]; then
        continue
    fi

    # Process each geo_accession only once
    if [[ -n "${seen_samples[$geo_accession]}" ]]; then
        continue  # Skip if already processed
    fi
    seen_samples["$geo_accession"]=1  # Mark as processed

    # Create a unique filename using Tissue, Sex, and Strain
    unique_name="${Tissue}_${Strain}_${Sex}_${geo_accession}"

    # Define output paths for this sample
    sample_output_dir="$PRJdir/$geo_accession"
    mkdir -p "$sample_output_dir"

    echo "Submitting JSON and BED jobs for $geo_accession" >> "$log_file"

    # Submit jobs with dependency on RSEMmatrix.sh completion
    for job_type in "2bigWigjson" "geneJson" "txBedJson" "GeneTPMbed" "txTPMbed"; do
        script_name=""
        case "$job_type" in
            "2bigWigjson") script_name="BWjson.sh" ;;
            "geneJson") script_name="gJSON.sh" ;;
            "txBedJson") script_name="tJSON.sh" ;;
            "GeneTPMbed") script_name="GeneTPMbed.sh" ;;
            "txTPMbed") script_name="txTPMbed.sh" ;;
        esac

        # Ensure matrix_job_id is valid before proceeding
        if [[ -z "$matrix_job_id" || ! "$matrix_job_id" =~ ^[0-9]+$ ]]; then
            echo "Skipping $job_type for $geo_accession due to invalid matrix_job_id" >> "$log_file"
            continue
        fi

        # Submit dependent job
        json_bed_job_ids["$geo_accession,$job_type"]=$(sbatch --dependency=afterok:$matrix_job_id --job-name="${job_type}_$geo_accession" \
            --export=BIOProjectID="$BIOProjectID",Run="$run",geo_accession="$geo_accession",PRJdir="$PRJdir",tissue="$Tissue",strain="$Strain",sex="$Sex",title="$Title",Sample_characteristics="$Sample_characteristics",PMID="$PMID",GEOpath="$GEOpath",unique_name="$unique_name",Logdir="$Logdir",scratch_dir="$scratch_dir",StrainInfo="$StrainInfo" \
            --output="$sample_output_dir/${job_type}-%j.out" --error="$sample_output_dir/${job_type}-%j.err" "$script_name" | awk '{print $4}')

        # Check job submission success
        if [[ -z "${json_bed_job_ids["$geo_accession,$job_type"]}" ]]; then
            echo "Error submitting $job_type job for $geo_accession" >> "$log_file"
        else
            echo "Submitted $job_type for $geo_accession with job ID: ${json_bed_job_ids["$geo_accession,$job_type"]}" >> "$log_file"
        fi
    done

done < "$AccList"

echo "All JSON and BED generation jobs have been submitted." >> "$log_file"
