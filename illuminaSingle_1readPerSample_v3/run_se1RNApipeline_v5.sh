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
Logdir="$myDir/data/expression/GEO/$BIOProjectID/log_files"
scratch_dir="/scratch/g/akwitek/wdemos/$BIOProjectID"
INDEX_DIR="$scratch_dir/RefIndex"
REFdir="/home/wdemos/NCBI_1August2023"
baseDir="$myDir/data/expression/GEO/$BIOProjectID"
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

STAR_log_dir="$Logdir/STAR"
mkdir -p "$STAR_log_dir"

if [ ! -d "$INDEX_DIR" ] || [ ! -f "$INDEX_DIR/SA" ] || [ ! -f "$INDEX_DIR/SAindex" ] || [ ! -f "$INDEX_DIR/Genome" ]; then
    echo "Submitting starRef_v3.sh to generate the STAR reference."
    starRef_job_id=$(sbatch --export=PRJdir="$PRJdir",Logdir="$Logdir",BIOProjectID="$BIOProjectID",length="$length" \
                        --output="$STAR_log_dir/starRef-%j.out" \
                        --error="$STAR_log_dir/starRef-%j.err" starRef_v3.sh | awk '{print $4}')
    if [ -z "$starRef_job_id" ]; then
        echo "Error: starRef_v3.sh submission failed."
        exit 1
    fi
    echo "STAR reference job ID: $starRef_job_id"
else
    echo "STAR reference already exists. Skipping starRef_v3.sh."
    starRef_job_id=""
fi

####################################
# STAR Alignment Jobs for Samples #
####################################

job_ids=()  # Array to store job IDs for STAR alignment
while IFS=$'\t' read -r run geo_accession Tissue Strain Sex PMID GEOpath Title Sample_characteristics StrainInfo; do
    if [ "$run" == "Run" ]; then continue; fi

    unique_name="${Tissue}_${Strain}_${Sex}_${geo_accession}"
    sample_output_dir="$PRJdir/$geo_accession"
    sample_log_dir="$sample_output_dir/log_files/STAR"
    mkdir -p "$sample_log_dir"

    echo "Submitting STAR job for $geo_accession"

    if [ -z "$starRef_job_id" ]; then
        job_id=$(sbatch --job-name="${geo_accession}_STAR72" \
                        --export=BIOProjectID="$BIOProjectID",Run="$run",geo_accession="$geo_accession",PRJdir="$PRJdir",length="$length",unique_name="$unique_name",Logdir="$Logdir" \
                        --output="$sample_log_dir/STAR-%j.out" \
                        --error="$sample_log_dir/STAR-%j.err" \
                        stars_v3.sh | awk '{print $4}')
    else
        job_id=$(sbatch --dependency=afterok:$starRef_job_id \
                        --job-name="${geo_accession}_STAR72" \
                        --export=BIOProjectID="$BIOProjectID",Run="$run",geo_accession="$geo_accession",PRJdir="$PRJdir",length="$length",unique_name="$unique_name",Logdir="$Logdir" \
                        --output="$sample_log_dir/STAR-%j.out" \
                        --error="$sample_log_dir/STAR-%j.err" \
                        stars_v3.sh | awk '{print $4}')
    fi

    if [ $? -ne 0 ]; then echo "Error submitting STAR job for $geo_accession"; continue; fi
    job_ids+=("$job_id")
done < "$AccList"

wait_for_jobs
##############################
# STAR Alignment Rate Check  #
##############################
echo "Check sample alignment rates for downstream processing."
STARQC_log_dir="$Logdir/STARQC"
mkdir -p "$STARQC_log_dir"

# Submit STARQC job
starqc_output=$(sbatch \
    --job-name="STARQC" \
    --export=Logdir="$Logdir",AccList="$AccList",BIOProjectID="$BIOProjectID",PRJdir="$PRJdir" \
    --output="$STARQC_log_dir/STARQC-%j.out" \
    --error="$STARQC_log_dir/STARQC-%j.err" \
    sSTARQC_v1.sh "$AccList" "$BIOProjectID" )

# Extract and validate job ID
starqc_job_id=$(echo "$starqc_output" | awk '{print $4}')
echo "sSTARQC_v1.sh job submission output: $starqc_output"

if [[ -z "$starqc_job_id" || ! "$starqc_job_id" =~ ^[0-9]+$ ]]; then
    echo "Error: Invalid STARQC job ID: '$starqc_job_id'"
    exit 1
fi

# Wait for STARQC to finish before consuming its output
wait_for_jobs "$starqc_job_id"
echo "sSTARQC_v1.sh completed. Output expected at: ${STARQC_log_dir}/${BIOProjectID}_STAR_Align_sum.txt"

#####################################################
# Filter to PASS samples for downstream processing  #
#####################################################
starqc_report="${Logdir}/STARQC/${BIOProjectID}_STAR_Align_sum.txt"
passAccList="${Logdir}/STARQC/${BIOProjectID}_Unique_AccList_PASS.txt"

# Sanity check: STARQC report exists
for i in {1..6}; do
  if [ -s "$starqc_report" ]; then break; fi
  echo "Waiting for STARQC report to appear... ($i/6)"
  sleep 5
done

if [ ! -s "$starqc_report" ]; then
  echo "ERROR: STARQC report not found or empty: $starqc_report"
  echo "Hint: check STARQC logs:"
  echo "  OUT: $STARQC_log_dir/STARQC-$starqc_job_id.out"
  echo "  ERR: $STARQC_log_dir/STARQC-$starqc_job_id.err"
  # Quick peek to help debugging (optional)
  tail -n +1 "$STARQC_log_dir/STARQC-$starqc_job_id".{out,err} 2>/dev/null | sed 's/^/STARQC log> /'
  exit 1
fi
# Build PASS-only accession list (keep header; filter by geo_accession in col 2)
awk 'NR==FNR { if (FNR>1 && $5=="PASS") pass[$1]=1; next }
     FNR==1 { print; next }
     ($2 in pass)' OFS='\t' "$starqc_report" "$AccList" > "$passAccList"

# Quick stats
total_unique=$(($(wc -l < "$AccList") - 1))
total_pass=$(($(wc -l < "$passAccList") - 1))
total_fail=$(( total_unique - total_pass ))

echo "STARQC gating: PASS=$total_pass  FAIL=$total_fail"
echo "PASS-only accession list saved to: $passAccList"

# If nothing passed, stop the pipeline early
if [ "$total_pass" -le 0 ]; then
    echo "ERROR: No samples passed alignment rate threshold. Exiting."
    exit 1
fi


##############################
# Sample Sex Estimation #
##############################
echo "Alignment with STAR is complete for all samples. Computing sample sex for samples with greater than 50% STAR alignment rate."

ComputeSex_log_dir="$Logdir/ComputeSex"
mkdir -p "$ComputeSex_log_dir"

# Sanity check
if [ ! -s "$passAccList" ]; then
    echo "ERROR: PASS-only accession list not found or empty: $passAccList"
    exit 1
fi

# Ensure the Logdir exists
mkdir -p "$Logdir"

# Submit ComputeSex job
job_output=$(sbatch \
    --job-name="computeSex" \
    --export=PRJdir="$PRJdir",Logdir="$Logdir",AccList="$passAccList",BIOProjectID="$BIOProjectID",baseDir="$baseDir" \
    --output="$ComputeSex_log_dir/ComputeSex-%j.out" \
    --error="$ComputeSex_log_dir/ComputeSex-%j.err" \
    ComputeSex_v3.sh "$AccList" "$BIOProjectID")

# Extract and validate job ID
compute_sex_job_id=$(echo "$job_output" | awk '{print $4}')
echo "ComputeSex job submission output: $job_output"

if [[ -z "$compute_sex_job_id" || ! "$compute_sex_job_id" =~ ^[0-9]+$ ]]; then
    echo "Error: Invalid ComputeSex job ID: '$compute_sex_job_id'"
    exit 1
fi

echo "ComputeSex_v3.sh submitted with job ID: $compute_sex_job_id"

# Wait for the ComputeSex job to complete
wait_for_jobs
echo "ComputeSex job completed. Continuing to next step..."

##############################
# RSEM Reference Preparation #
##############################
RSEM_log_dir="$Logdir/RSEM"
mkdir -p "$RSEM_log_dir"


if [ -f "$scratch_dir/rsemref.chrlist" ] && [ -f "$scratch_dir/rsemref.grp" ] && [ -f "$scratch_dir/rsemref.idx.fa" ] && [ -f "$scratch_dir/rsemref.n2g.idx.fa" ] && [ -f "$scratch_dir/rsemref.seq" ] && [ -f "$scratch_dir/rsemref.ti" ] && [ -f "$scratch_dir/rsemref.transcripts.fa" ]; then
    echo "RSEM reference already exists. Skipping RSEMref_v3.sh."
else
    echo "Submitting RSEMref_v3.sh to generate the RSEM reference."
    rsemRef_job_id=$(sbatch --export=BIOProjectID="$BIOProjectID" \
                           --output="$RSEM_log_dir/RSEMref-%j.out" \
                           --error="$RSEM_log_dir/RSEMref-%j.err" RSEMref_v3.sh | awk '{print $4}')
    if [ -z "$rsemRef_job_id" ]; then echo "Error: RSEMref_v3.sh submission failed."; exit 1; fi
    echo "RSEM reference job ID: $rsemRef_job_id"
fi

###########################
# RSEM Jobs for Samples #
###########################
job_ids=()
while IFS=$'\t' read -r run geo_accession Tissue Strain Sex PMID GEOpath Title Sample_characteristics StrainInfo; do
    if [ "$run" == "Run" ]; then continue; fi
    sample_output_dir="$PRJdir/$geo_accession"
    sample_log_dir="$sample_output_dir/log_files/RSEM"
    mkdir -p "$sample_log_dir"

    echo "Submitting RSEM job for $geo_accession"

    if [ -z "$rsemRef_job_id" ]; then
        job_id=$(sbatch --job-name="${geo_accession}_RSEM72" \
                        --export=BIOProjectID="$BIOProjectID",Run="$run",geo_accession="$geo_accession",PRJdir="$PRJdir",Logdir="$Logdir" \
                        --output="$sample_log_dir/RSEM-%j.out" \
                        --error="$sample_log_dir/RSEM-%j.err" \
                        RSEMs_v3.sh | awk '{print $4}')
    else
        job_id=$(sbatch --dependency=afterok:$rsemRef_job_id \
                        --job-name="${geo_accession}_RSEM72" \
                        --export=BIOProjectID="$BIOProjectID",Run="$run",geo_accession="$geo_accession",PRJdir="$PRJdir",Logdir="$Logdir" \
                        --output="$sample_log_dir/RSEM-%j.out" \
                        --error="$sample_log_dir/RSEM-%j.err" \
                        RSEMs_v3.sh | awk '{print $4}')
    fi

    # Check if sbatch failed
    if [ $? -ne 0 ]; then
        echo "Error submitting RSEM job for $geo_accession"
        continue
    fi
    job_ids+=("$job_id")
done < "$passAccList"

wait_for_jobs

#####################################
# Final Steps - Matrix and QC files #
#####################################
echo "Generating combined matrices for all samples and final MultiQC report."

matrix_log_dir="$Logdir/RSEM"
mkdir -p "$matrix_log_dir"

job_id_output=$(sbatch --export=Run="$run",geo_accession="$geo_accession",PRJdir="$PRJdir",Logdir="$Logdir",AccList="$passAccList",BIOProjectID="$BIOProjectID",baseDir="$baseDir" \
              --output="$matrix_log_dir/RSEMmatrix-%j.out" \
              --error="$matrix_log_dir/RSEMmatrix-%j.err" RSEMmatrixs_v3.sh 2>&1)

matrix_job_id=$(echo "$job_id_output" | awk '{print $4}')
echo "RSEMmatrixs_v3.sh submission output: $job_id_output" >> "$log_file"

if [[ -z "$matrix_job_id" || ! "$matrix_job_id" =~ ^[0-9]+$ ]]; then
    echo "Error: Invalid matrix_job_id: '$matrix_job_id'" >> "$log_file"
    exit 1
fi

echo "RSEMmatrixs_v3.sh submitted with job ID: $matrix_job_id" >> "$log_file"

echo "Waiting for RSEMmatrix job to complete..." >> "$log_file"
wait_for_jobs
echo "RSEMmatrix job completed. Proceeding with JSON and BED generation..." >> "$log_file"

###################################
# Generate BED and JSON files
###################################

declare -A json_bed_job_ids
declare -A seen_samples

while IFS=$'\t' read -r run geo_accession Tissue Strain Sex PMID GEOpath Title Sample_characteristics StrainInfo; do
    if [ "$run" == "Run" ]; then continue; fi

    if [[ -n "${seen_samples[$geo_accession]}" ]]; then continue; fi
    seen_samples["$geo_accession"]=1

    unique_name="${Tissue}_${Strain}_${Sex}_${geo_accession}"
    sample_output_dir="$PRJdir/$geo_accession"
    sample_log_dir="$sample_output_dir/log_files/bedsNjsons"

    mkdir -p "$sample_log_dir"

    echo "Submitting JSON and BED jobs for $geo_accession" >> "$log_file"

    for job_type in "2bigWigjson" "geneJson" "txBedJson" "GeneTPMbed" "txTPMbed"; do
        case "$job_type" in
            "2bigWigjson") script_name="BWjson_v3.sh" ;;
            "geneJson") script_name="gJSON_v3.sh" ;;
            "txBedJson") script_name="tJSON_v3.sh" ;;
            "GeneTPMbed") script_name="GeneTPMbed_v3.sh" ;;
            "txTPMbed") script_name="txTPMbed_v3.sh" ;;
        esac

        if [[ -z "$matrix_job_id" || ! "$matrix_job_id" =~ ^[0-9]+$ ]]; then
            echo "Skipping $job_type for $geo_accession due to invalid matrix_job_id" >> "$log_file"
            continue
        fi

        json_bed_job_ids["$geo_accession,$job_type"]=$(sbatch --dependency=afterok:$matrix_job_id \
            --job-name="${job_type}_$geo_accession" \
            --export=BIOProjectID="$BIOProjectID",Run="$run",geo_accession="$geo_accession",baseDir="$baseDir",PRJdir="$PRJdir",tissue="$Tissue",strain="$Strain",sex="$Sex",title="$Title",Sample_characteristics="$Sample_characteristics",PMID="$PMID",GEOpath="$GEOpath",unique_name="$unique_name",Logdir="$Logdir",scratch_dir="$scratch_dir",StrainInfo="$StrainInfo" \
            --output="$sample_log_dir/${job_type}-%j.out" \
            --error="$sample_log_dir/${job_type}-%j.err" \
            "$script_name" | awk '{print $4}')

        if [[ -z "${json_bed_job_ids["$geo_accession,$job_type"]}" ]]; then
            echo "Error submitting $job_type job for $geo_accession" >> "$log_file"
        else
            echo "Submitted $job_type for $geo_accession with job ID: ${json_bed_job_ids["$geo_accession,$job_type"]}" >> "$log_file"
        fi
    done

done < "$passAccList"

echo "All JSON and BED generation jobs have been submitted." >> "$log_file"
