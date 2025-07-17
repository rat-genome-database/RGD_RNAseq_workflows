#!/bin/sh

# Check command line arguments
if [ $# -ne 3 ]; then
    echo "Usage: $0 <SRA_accession_list.txt> <BioProject_ID> <Read_Length>"
    exit 1
fi

AccList=$1  # Accession list file
BIOProjectID=$2  # Project ID for output directory
length=$3
scratch_dir="/scratch/g/akwitek/wdemos/${BIOProjectID}"
log_file="${scratch_dir}/run_combination_log.txt"
INDEX_DIR="${scratch_dir}/RefIndex"
Logdir="${scratch_dir}/logs"
PRJdir="/home/wdemos/data/expression/GEO/${BIOProjectID}/reads_fastq"


# Ensure necessary directories exist
mkdir -p "$scratch_dir" "$Logdir"

# Clear the log file before starting
echo "Log of combined runs for each sample" > "$log_file"

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

###########################
# STAR Alignment Jobs for Each Sample
###########################

temp_file=$(mktemp)
trap 'rm -f "$temp_file"' EXIT

cut -f2 "$AccList" | tail -n +2 | sort | uniq | grep -v '^$' > "$temp_file"
mapfile -t geo_accessions < "$temp_file"

if [ ! -x 3STAR_multiSRR.sh ]; then
    echo "ERROR: 3STAR_multiSRR.sh not found or not executable. Exiting."
    exit 1
fi

star_job_ids=()

for geo_accession in "${geo_accessions[@]}"; do
    RUNS=$(awk -v sample="$geo_accession" -F'\t' '$2 == sample {print $1}' "$AccList")
    read Tissue Strain Sex <<< $(awk -v sample="$geo_accession" -F'\t' '$2 == sample {print $3, $4, $5; exit}' "$AccList")
    unique_name="${Tissue}_${Strain}_${Sex}_${geo_accession}"

    READ1_FILES=$(for run in $RUNS; do ls "${scratch_dir}/${run}"/*_1.fastq.gz 2>/dev/null; done | paste -sd ",")
    READ2_FILES=$(for run in $RUNS; do ls "${scratch_dir}/${run}"/*_2.fastq.gz 2>/dev/null; done | paste -sd ",")

    echo "Processing Sample: $geo_accession"
    echo "  Runs: $RUNS"
    echo "  READ1: $READ1_FILES"
    echo "  READ2: $READ2_FILES" >> "$log_file"

    if [[ -n "$READ1_FILES" && -n "$READ2_FILES" ]]; then
        echo "  Submitting STAR alignment job..."

        if [ -n "$starRef_job_id" ]; then
            sbatch_output=$(sbatch --dependency=afterok:${starRef_job_id} \
                            3STAR_multiSRR.sh "$geo_accession" "$READ1_FILES" "$READ2_FILES" "$BIOProjectID" "$unique_name")
        else
            sbatch_output=$(sbatch \
                            3STAR_multiSRR.sh "$geo_accession" "$READ1_FILES" "$READ2_FILES" "$BIOProjectID" "$unique_name")
        fi

        echo "  sbatch output: $sbatch_output"
        job_id=$(echo "$sbatch_output" | awk '{print $4}')

        if [[ "$job_id" =~ ^[0-9]+$ ]]; then
            echo "  STAR job ID: $job_id"
            star_job_ids+=("$job_id")
        else
            echo "  Warning: Failed to extract valid job ID for $geo_accession"
        fi
    else
        echo "  Warning: No FASTQ files found for $geo_accession" >> "$log_file"
    fi
done

# Final check on job IDs
valid_star_job_ids=()
for id in "${star_job_ids[@]}"; do
    if [[ "$id" =~ ^[0-9]+$ ]]; then
        valid_star_job_ids+=("$id")
    fi
done

if [ ${#valid_star_job_ids[@]} -eq 0 ]; then
    echo "ERROR: No valid STAR jobs submitted. Exiting."
    exit 1
fi

star_job_ids_str=$(IFS=:; echo "${valid_star_job_ids[*]}")
echo "STAR job dependency string: $star_job_ids_str"

##############################
# Sample Sex Estimation
##############################

if [ ! -x ComputeSex.sh ]; then
    echo "ERROR: ComputeSex.sh not found or not executable. Exiting."
    exit 1
fi

echo "Submitting ComputeSex.sh job with dependency on STAR alignment jobs..."

sbatch_output=$(sbatch --dependency=afterok:${star_job_ids_str} \
    --export=PRJdir="$PRJdir",Logdir="$Logdir",AccList="$AccList",BIOProjectID="$BIOProjectID" \
    --output="$Logdir/ComputeSex-%j.out" \
    --error="$Logdir/ComputeSex-%j.err" ComputeSex.sh "$BIOProjectID")

echo "ComputeSex sbatch output: $sbatch_output"
sex_job_id=$(echo "$sbatch_output" | awk '{print $4}')

if [[ ! "$sex_job_id" =~ ^[0-9]+$ ]]; then
    echo "ERROR: ComputeSex.sh submission failed. SLURM output: $sbatch_output"
    exit 1
fi

echo "ComputeSex job ID: $sex_job_id"

##############################
# RSEM Reference Preparation
##############################

if [ -f "$scratch_dir/rsemref.chrlist" ] && [ -f "$scratch_dir/rsemref.grp" ] && \
   [ -f "$scratch_dir/rsemref.idx.fa" ] && [ -f "$scratch_dir/rsemref.n2g.idx.fa" ] && \
   [ -f "$scratch_dir/rsemref.seq" ] && [ -f "$scratch_dir/rsemref.ti" ] && \
   [ -f "$scratch_dir/rsemref.transcripts.fa" ]; then
    echo "RSEM reference already exists. Skipping RSEMref.sh."
else
    echo "RSEM reference files not found. Submitting RSEMref.sh..."

    if [ -z "$star_job_ids_str" ]; then
        echo "ERROR: Dependency string for RSEMref is empty. Exiting."
        exit 1
    fi

    if [ ! -x RSEMref.sh ]; then
        echo "ERROR: RSEMref.sh not found or not executable. Exiting."
        exit 1
    fi

    mkdir -p "$Logdir"

    sbatch_output=$(sbatch --dependency=afterok:${star_job_ids_str} \
        --export=BIOProjectID="$BIOProjectID" \
        --output="$Logdir/RSEMref-%j.out" \
        --error="$Logdir/RSEMref-%j.err" RSEMref.sh)

    echo "RSEMref sbatch output: $sbatch_output"
    rsemRef_job_id=$(echo "$sbatch_output" | awk '{print $4}')

    if [[ ! "$rsemRef_job_id" =~ ^[0-9]+$ ]]; then
        echo "ERROR: RSEMref.sh submission failed. SLURM output: $sbatch_output"
        exit 1
    fi

    echo "RSEMref job ID: $rsemRef_job_id"
fi

###########################
# RSEM Jobs for Samples #
###########################
job_ids=()  # Array to store job IDs for RSEM
declare -A processed_samples  # Associative array to track processed geo_accessions

while IFS=$'\t' read -r run geo_accession Tissue Strain Sex PMID GEOpath Title Sample_characteristics StrainInfo; do
    # Skip the header row
    if [ "$run" == "Run" ]; then
        continue
    fi

    # Skip if this geo_accession has already been processed
    if [[ -n "${processed_samples[$geo_accession]}" ]]; then
        continue
    fi

    sample_output_dir="$PRJdir/$geo_accession"
    mkdir -p "$sample_output_dir"  # Ensure directory exists

    echo "Submitting RSEM job for $geo_accession"

    if [ -z "$rsemRef_job_id" ]; then
        job_id=$(sbatch --job-name="${geo_accession}_RSEM72" \
                        --export=BIOProjectID="$BIOProjectID",geo_accession="$geo_accession",PRJdir="$PRJdir",Logdir="$Logdir" \
                        --output="$sample_output_dir/RSEM-%j.out" \
                        --error="$sample_output_dir/RSEM-%j.err" \
                        RSEM.sh | awk '{print $4}')
    else
        job_id=$(sbatch --dependency=afterok:$rsemRef_job_id \
                        --job-name="${geo_accession}_RSEM72" \
                        --export=BIOProjectID="$BIOProjectID",geo_accession="$geo_accession",PRJdir="$PRJdir",Logdir="$Logdir" \
                        --output="$sample_output_dir/RSEM-%j.out" \
                        --error="$sample_output_dir/RSEM-%j.err" \
                        RSEM.sh | awk '{print $4}')
    fi

    # Check if sbatch failed
    if [ -z "$job_id" ]; then
        echo "Error submitting RSEM job for $geo_accession"
        continue
    fi

    job_ids+=("$job_id")
    processed_samples["$geo_accession"]=1  # Mark sample as processed

done < "$AccList"

wait_for_jobs

#######################
# Matrix and QC files #
#######################
echo "Generating combined matrices for all samples and final MultiQC report."
matrix_job_id=$(sbatch --export=Run="$run",geo_accession="$geo_accession",PRJdir="$PRJdir",Logdir="$Logdir",AccList="$AccList",BIOProjectID="$BIOProjectID" \
              --output="$Logdir/RSEMmatrix-%j.out" \
              --error="$Logdir/RSEMmatrix-%j.err" 2RSEMmatrix.sh | awk '{print $4}')

if [ -z "$matrix_job_id" ]; then
    echo "Error running 2RSEMmatrix.sh"
    exit 1
fi

wait_for_jobs

###################################
# Generate bed and JSON files #
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
        continue  # Skip if this geo_accession was already processed
    fi
    seen_samples["$geo_accession"]=1  # Mark as processed

    # Create a unique filename using Tissue, Sex, and Strain
    unique_name="${Tissue}_${Strain}_${Sex}_${geo_accession}"

    # Define output paths for this sample
    sample_output_dir="$PRJdir/$geo_accession"
    mkdir -p "$sample_output_dir"

    echo "Submitting JSON and BED jobs for $geo_accession" >> "$log_file"

    # Submit jobs with dependency on 2RSEMmatrix.sh completion
    for job_type in "2bigWigjson" "geneJson" "txBedJson" "GeneTPMbed" "txTPMbed"; do
        script_name=""
        case "$job_type" in
            "2bigWigjson") script_name="BWjson.sh" ;;
            "geneJson") script_name="gJSON.sh" ;;
            "txBedJson") script_name="tJSON.sh" ;;
            "GeneTPMbed") script_name="GeneTPMbed.sh" ;;
            "txTPMbed") script_name="txTPMbed.sh" ;;
        esac

        json_bed_job_ids["$geo_accession,$job_type"]=$(sbatch --dependency=afterok:$matrix_job_id --job-name="${job_type}_$geo_accession" \
            --export=BIOProjectID="$BIOProjectID",Run="$run",geo_accession="$geo_accession",PRJdir="$PRJdir",tissue="$Tissue",strain="$Strain",sex="$Sex",title="$Title",Sample_characteristics="$Sample_characteristics",PMID="$PMID",GEOpath="$GEOpath",unique_name="$unique_name",Logdir="$Logdir",scratch_dir="$scratch_dir",StrainInfo="$StrainInfo" \
            --output="$sample_output_dir/${job_type}-%j.out" --error="$sample_output_dir/${job_type}-%j.err" "$script_name" | awk '{print $4}')

        if [ -z "${json_bed_job_ids["$geo_accession,$job_type"]}" ]; then
            echo "Error submitting $job_type job for $geo_accession" >> "$log_file"
        else
            echo "Submitted $job_type for $geo_accession with job ID: ${json_bed_job_ids["$geo_accession,$job_type"]}" >> "$log_file"
        fi
    done

done < "$AccList"

echo "All JSON and BED generation jobs have been submitted." >> "$log_file"
