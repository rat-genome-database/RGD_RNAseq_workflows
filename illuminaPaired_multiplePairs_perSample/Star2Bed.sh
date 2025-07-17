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

######################
# STAR Alignment Jobs for Each Sample #
######################
cut -f2 "$AccList" | tail -n +2 | sort | uniq | while read geo_accession; do
    # Get all runs associated with this geo_accession
    RUNS=$(awk -v sample="$geo_accession" -F'\t' '$2 == sample {print $1}' "$AccList")

    # Extract metadata (Tissue, Strain, Sex) for the sample
    read Tissue Strain Sex <<< $(awk -v sample="$geo_accession" -F'\t' '$2 == sample {print $3, $4, $5; exit}' "$AccList")

    # Construct the unique name
    unique_name="${Tissue}_${Strain}_${Sex}_${geo_accession}"

    # Collect all READ1 and READ2 FASTQ files
    READ1_FILES=$(for run in $RUNS; do ls ${scratch_dir}/${run}/*_1.fastq.gz 2>/dev/null; done | paste -sd ",")
    READ2_FILES=$(for run in $RUNS; do ls ${scratch_dir}/${run}/*_2.fastq.gz 2>/dev/null; done | paste -sd ",")

    # Log the combined runs
    echo "Sample: ${geo_accession}, Runs: ${RUNS}" >> "$log_file"
    echo "these are $READ1_FILES" >> "$log_file"
    echo "these are $READ2_FILES" >> "$log_file"
    # Ensure we have reads before submitting
    if [[ -n "$READ1_FILES" && -n "$READ2_FILES" ]]; then
        echo "Submitting STAR alignment for Sample: ${geo_accession}"
        star_job_id=$(sbatch --dependency=afterok:${starRef_job_id} 2STAR_multiSRR.sh "${geo_accession}" "${READ1_FILES}" "${READ2_FILES}" "${BIOProjectID}" "${unique_name}" | awk '{print $4}')
        star_job_ids+=($star_job_id)  # Store the job IDs for STAR alignments
    else
        echo "Warning: No FASTQ files found for ${geo_accession}" >> "$log_file"
    fi
done
wait_for_jobs

##############################
# Sample Sex Estimation #
##############################
echo "Alignment with STAR is complete for all samples. Computing sample sex."

# Submit the ComputeSex job only after all STAR alignments are complete
job_id=$(sbatch --dependency=afterok:${star_job_ids[*]} \
    --export=PRJdir="$PRJdir",Logdir="$Logdir",AccList="$AccList",BIOProjectID="$BIOProjectID" \
    --output="$Logdir/ComputeSex-%j.out" \
    --error="$Logdir/ComputeSex-%j.err" ComputeSex.sh "$BIOProjectID" | awk '{print $4}')

wait_for_jobs

##############################
# RSEM Reference Preparation #
##############################
if [ -f "$scratch_dir/rsemref.chrlist" ] && [ -f "$scratch_dir/rsemref.grp" ] && [ -f "$scratch_dir/rsemref.idx.fa" ] && [ -f "$scratch_dir/rsemref.n2g.idx.fa" ] && [ -f "$scratch_dir/rsemref.seq" ] && [ -f "$scratch_dir/rsemref.ti" ] && [ -f "$scratch_dir/rsemref.transcripts.fa" ]; then
    echo "RSEM reference already exists. Skipping RSEMref.sh."
else
    echo "Submitting RSEMref.sh to generate the RSEM reference."
    rsemRef_job_id=$(sbatch --dependency=afterok:${star_job_ids[*]} \
                           --export=BIOProjectID="$BIOProjectID" \
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
    if [ $? -ne 0 ]; then
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
              --error="$Logdir/RSEMmatrix-%j.err" 2RSEMmatrix.sh)
if [ $? -ne 0 ]; then
    echo "Error running 2RSEMmatrix.sh"
    exit 1
fi
wait_for_jobs

###################################
# Generate bed and JSON files #
###################################

# Read accessions and prepare jobs for each unique sample
declare -A job_ids
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

    # Debugging: Print the value of `Run` before submitting jobs
    echo "Debug: Submitting jobs for geo_accession $geo_accession" >> "$log_file"

    # Submit bigWig JSON job with dependency on 2RSEMmatrix.sh completion
    job_ids["$geo_accession,2bigWigjson"]=$(sbatch --dependency=afterok:$matrix_job_id --job-name="bgJ_$geo_accession" \
        --export=BIOProjectID="$BIOProjectID",Run="$run",geo_accession="$geo_accession",PRJdir="$PRJdir",tissue="$Tissue",strain="$Strain",sex="$Sex",title="$Title",Sample_characteristics="$Sample_characteristics",PMID="$PMID",GEOpath="$GEOpath",unique_name="$unique_name",Logdir="$Logdir",scratch_dir="$scratch_dir",StrainInfo="$StrainInfo" \
        --output="$sample_output_dir/bwjson-%j.out" --error="$sample_output_dir/bwjson-%j.err" BWjson.sh | awk '{print $4}')

    if [ $? -ne 0 ]; then
        echo "Error submitting bigWig JSON job for $geo_accession" >> "$log_file"
        continue
    fi

    # Submit gene TPM bed JSON job with dependency on 2RSEMmatrix.sh completion
    job_ids["$geo_accession,geneJson"]=$(sbatch --dependency=afterok:$matrix_job_id --job-name="gj_$geo_accession" \
        --export=BIOProjectID="$BIOProjectID",Run="$run",geo_accession="$geo_accession",PRJdir="$PRJdir",tissue="$Tissue",strain="$Strain",sex="$Sex",title="$Title",Sample_characteristics="$Sample_characteristics",PMID="$PMID",GEOpath="$GEOpath",unique_name="$unique_name",Logdir="$Logdir",scratch_dir="$scratch_dir",StrainInfo="$StrainInfo" \
        --output="$sample_output_dir/GbedJSON-%j.out" --error="$sample_output_dir/GbedJSON-%j.err" gJSON.sh | awk '{print $4}')

    if [ $? -ne 0 ]; then
        echo "Error submitting gene TPM bed JSON job for $geo_accession" >> "$log_file"
        continue
    fi

    # Submit transcript TPM bed JSON job with dependency on 2RSEMmatrix.sh completion
    job_ids["$geo_accession,txBedJson"]=$(sbatch --dependency=afterok:$matrix_job_id --job-name="tj_$geo_accession" \
        --export=BIOProjectID="$BIOProjectID",Run="$run",geo_accession="$geo_accession",PRJdir="$PRJdir",tissue="$Tissue",strain="$Strain",sex="$Sex",title="$Title",Sample_characteristics="$Sample_characteristics",PMID="$PMID",GEOpath="$GEOpath",unique_name="$unique_name",Logdir="$Logdir",scratch_dir="$scratch_dir",StrainInfo="$StrainInfo" \
        --output="$sample_output_dir/tbedJSON-%j.out" --error="$sample_output_dir/GbedJSON-%j.err" tJSON.sh | awk '{print $4}')
    if [ $? -ne 0 ]; then
        echo "Error submitting transcript bed JSON job for $geo_accession" >> "$log_file"
        continue
    fi

    # Submit gene TPM bed job with dependency on 2RSEMmatrix.sh completion
    job_ids["$geo_accession,GeneTPMbed"]=$(sbatch --dependency=afterok:$matrix_job_id --job-name="gb_$geo_accession" \
        --export=scratch_dir="$scratch_dir",BIOProjectID="$BIOProjectID",geo_accession="$geo_accession",PRJdir="$PRJdir",tissue="$Tissue",strain="$Strain",sex="$Sex",unique_name="$unique_name",Logdir="$Logdir" \
        --output="$log_file" --error="$log_file" GeneTPMbed.sh | awk '{print $4}')

    if [ $? -ne 0 ]; then
        echo "Error submitting gene TPM bed job for $geo_accession" >> "$log_file"
        continue
    fi

    # Submit transcript TPM bed job with dependency on 2RSEMmatrix.sh completion
    job_ids["$geo_accession,txTPMbed"]=$(sbatch --dependency=afterok:$matrix_job_id --job-name="txb_$geo_accession" \
        --export=scratch_dir="$scratch_dir",BIOProjectID="$BIOProjectID",geo_accession="$geo_accession",PRJdir="$PRJdir",tissue="$Tissue",strain="$Strain",sex="$Sex",unique_name="$unique_name",Logdir="$Logdir" \
        --output="$log_file" --error="$log_file" txTPMbed.sh | awk '{print $4}')
    echo "Jobs submitted for geo_accession $geo_accession: BigWig JSON, Gene JSON, and Transcript JSON." >> "$log_file"

done < "$AccList"

echo "All JSON generation jobs have been submitted." >> "$log_file"
