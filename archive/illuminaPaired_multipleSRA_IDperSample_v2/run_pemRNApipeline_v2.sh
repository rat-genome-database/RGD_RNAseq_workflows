#!/bin/bash
#SBATCH --job-name=mRNA
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=06:00:00
#SBATCH --account=akwitek
#SBATCH --output=%x-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wdemos@mcw.edu
#SBATCH --error=%x-%j.err

# Enable debugging by printing each command before execution
set -x
### 12 Aug 2025 this script must be called with "bash" and not "sh" due to awk functions. If it fails, check first to see how the script was called.
# Check command line arguments
if [ $# -ne 3 ]; then
    echo "Usage: $0 <SRA_accession_list.txt> <BioProject_ID> <Read_Length>"
    exit 1
fi

# Set inputs and directories
AccList=$1  # Accession list file
BIOProjectID=$2  # Project ID for output directory
length=$3
myDir="/home/wdemos"
scratch_dir="/scratch/g/akwitek/wdemos/${BIOProjectID}"
baseDir="${myDir}/data/expression/GEO/${BIOProjectID}"
Logdir="${baseDir}/log_files"
PRJdir="${baseDir}/reads_fastq"  # This is where individual sample output lands
INDEX_DIR="${scratch_dir}/RefIndex"

# Ensure necessary directories exist
mkdir -p "$scratch_dir" "$baseDir" "$PRJdir" "$Logdir"

# Function to wait for jobs to complete
wait_for_jobs() {
    for job_id in "$@"; do
        while squeue -u $USER | grep -q "$job_id"; do
            sleep 5
        done
    done
}

##############################
# STAR Reference Preparation #
##############################

StarRef_log_dir="$Logdir/STAR"
mkdir -p "$StarRef_log_dir"

# Submit STAR genome generation if not already done
if [ ! -d "$INDEX_DIR" ] || [ ! -f "$INDEX_DIR/SA" ] || [ ! -f "$INDEX_DIR/SAindex" ] || [ ! -f "$INDEX_DIR/Genome" ]; then
    echo "Submitting starRef_v2.sh to generate the STAR reference."
    starRef_job_id=$(sbatch --export=PRJdir="$PRJdir",Logdir="$Logdir",BIOProjectID="$BIOProjectID",length="$length" \
                        --output="$StarRef_log_dir/starRef-%j.out" \
                        --error="$StarRef_log_dir/starRef-%j.err" starRef_v2.sh | awk '{print $4}')
    if [ -z "$starRef_job_id" ]; then
        echo "Error: starRef_v2.sh submission failed."
        exit 1
    fi
    echo "STAR reference job ID: $starRef_job_id"
else
    echo "STAR reference already exists. Skipping starRef_v2.sh."
    starRef_job_id=""
fi

###########################
# STAR Alignment Jobs for Each Sample
###########################

temp_file=$(mktemp)
trap 'rm -f "$temp_file"' EXIT

cut -f2 "$AccList" | tail -n +2 | sort | uniq | grep -v '^$' > "$temp_file"
mapfile -t geo_accessions < "$temp_file"

if [ ! -x STAR_mv2.sh ]; then
    echo "ERROR: STAR_mv2.sh not found or not executable. Please check the script path and permissions."
    exit 1
fi

star_job_ids=()

for geo_accession in "${geo_accessions[@]}"; do
    RUNS=$(awk -v sample="$geo_accession" -F'\t' '$2 == sample {print $1}' "$AccList")
    read -r Tissue Strain Sex <<< $(awk -v sample="$geo_accession" -F'\t' '$2 == sample {print $3, $4, $5; exit}' "$AccList")
    unique_name="${Tissue}_${Strain}_${Sex}_${geo_accession}"

    READ1_FILES=$(for run in $RUNS; do ls "${scratch_dir}/${run}"/*_1.fastq.gz 2>/dev/null; done | paste -sd ",")
    READ2_FILES=$(for run in $RUNS; do ls "${scratch_dir}/${run}"/*_2.fastq.gz 2>/dev/null; done | paste -sd ",")

    echo "$(date): Processing Sample: $geo_accession"
    echo "  Runs: $RUNS"
    echo "  READ1: $READ1_FILES"
    echo "  READ2: $READ2_FILES"

    if [[ -n "$READ1_FILES" && -n "$READ2_FILES" ]]; then
        echo "  Submitting STAR alignment job..."

        sample_log_dir="$PRJdir/$geo_accession/log_files/STAR"
        mkdir -p "$sample_log_dir"

        if [ -n "$starRef_job_id" ]; then
            sbatch_output=$(sbatch --job-name="STAR_$geo_accession" \
                --export=baseDir="$baseDir",PRJdir="$PRJdir",Logdir="$Logdir" \
                --output="$sample_log_dir/STAR-%j.out" \
                --error="$sample_log_dir/STAR-%j.err" \
                --dependency=afterok:${starRef_job_id} \
                STAR_mv2.sh "$geo_accession" "$READ1_FILES" "$READ2_FILES" "$BIOProjectID" "$unique_name")
        else
            sbatch_output=$(sbatch --job-name="STAR_$geo_accession" \
                --export=baseDir="$baseDir",PRJdir="$PRJdir",Logdir="$Logdir" \
                --output="$sample_log_dir/STAR-%j.out" \
                --error="$sample_log_dir/STAR-%j.err" \
                STAR_mv2.sh "$geo_accession" "$READ1_FILES" "$READ2_FILES" "$BIOProjectID" "$unique_name")
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
        echo "  Warning: No FASTQ files found for $geo_accession. Skipping STAR alignment."
        echo "  Runs attempted: $RUNS"
        for run in $RUNS; do
            echo "    Checked: ${scratch_dir}/${run}/*_1.fastq.gz"
            echo "             ${scratch_dir}/${run}/*_2.fastq.gz"
        done
    fi

done

# Wait for all STAR jobs to complete
if [ ${#star_job_ids[@]} -eq 0 ]; then
    echo "ERROR: No valid STAR jobs submitted. Exiting."
    exit 1
fi

echo "Waiting for all STAR jobs to complete before proceeding..."
wait_for_jobs "${star_job_ids[@]}"
echo "All STAR jobs completed."

##############################
# Create Deduplicated Accession List
##############################

uniqueAccList="${Logdir}/${BIOProjectID}_Unique_AccList.txt"
{
    head -n 1 "$AccList"
    tail -n +2 "$AccList" | sort -u -t$'\t' -k2,2
} > "$uniqueAccList"
echo "Unique accession list saved to $uniqueAccList"

##############################
# Sample Sex Estimation      #
##############################

echo "Computing sample sex."

# Submit ComputeSex job
job_output=$(sbatch \
    --job-name="computeSex" \
    --export=PRJdir="$PRJdir",Logdir="$Logdir",AccList="$uniqueAccList",BIOProjectID="$BIOProjectID",baseDir="$baseDir" \
    --output="$Logdir/ComputeSex-%j.out" \
    --error="$Logdir/ComputeSex-%j.err" \
    ComputeSex_v2.sh "$uniqueAccList" "$BIOProjectID")

# Extract and validate job ID
compute_sex_job_id=$(echo "$job_output" | awk '{print $4}')
echo "ComputeSex_v5 job submission output: $job_output"

if [[ -z "$compute_sex_job_id" || ! "$compute_sex_job_id" =~ ^[0-9]+$ ]]; then
    echo "Error: Invalid ComputeSex job ID: '$compute_sex_job_id'"
    exit 1
fi

echo "ComputeSex_v2.sh submitted with job ID: $compute_sex_job_id"

# Wait for the ComputeSex job to complete
wait_for_jobs "$compute_sex_job_id"
echo "ComputeSex_v5 job completed. Continuing to next step..."
##############################
# RSEM Reference Preparation #
##############################
if [ -f "$scratch_dir/rsemref.chrlist" ] && [ -f "$scratch_dir/rsemref.grp" ] && [ -f "$scratch_dir/rsemref.idx.fa" ] && [ -f "$scratch_dir/rsemref.n2g.idx.fa" ] && [ -f "$scratch_dir/rsemref.seq" ] && [ -f "$scratch_dir/rsemref.ti" ] && [ -f "$scratch_dir/rsemref.transcripts.fa" ]; then
    echo "RSEM reference already exists. Skipping RSEMref.sh."
else
    echo "Submitting RSEMref.sh to generate the RSEM reference."
    rsemRef_job_id=$(sbatch --export=BIOProjectID="$BIOProjectID" \
                           --output="$Logdir/RSEMref-%j.out" \
                           --error="$Logdir/RSEMref-%j.err" RSEMref_v2.sh | awk '{print $4}')
    if [ -z "$rsemRef_job_id" ]; then echo "Error: RSEMref_v2.sh submission failed."; exit 1; fi
    echo "RSEM reference job ID: $rsemRef_job_id"
fi

###########################
# RSEM Jobs for Samples #
###########################
# Get unique geo_accessions from AccList (skip header)
mapfile -t unique_geo_accessions < <(awk -F'\t' 'NR>1 && $2 != "" {print $2}' "$AccList" | sort -u)

job_ids=()

for geo_accession in "${unique_geo_accessions[@]}"; do
    # Get the first Run associated with this geo_accession
    run=$(awk -F'\t' -v geo="$geo_accession" '$2 == geo {print $1; exit}' "$AccList")

    sample_output_dir="$PRJdir/$geo_accession"
    sample_log_dir="$sample_output_dir/log_files/RSEM"
    mkdir -p "$sample_log_dir"

    echo "Submitting RSEM job for sample $geo_accession (run: $run)"

    if [ -z "$rsemRef_job_id" ]; then
        job_output=$(sbatch --job-name="${geo_accession}_RSEM" \
                           --export=BIOProjectID="$BIOProjectID",Run="$run",geo_accession="$geo_accession",PRJdir="$PRJdir",Logdir="$Logdir" \
                           --output="$sample_log_dir/RSEM-%j.out" \
                           --error="$sample_log_dir/RSEM-%j.err" \
                           RSEM_v2.sh)
    else
        job_output=$(sbatch --dependency=afterok:$rsemRef_job_id \
                           --job-name="${geo_accession}_RSEM" \
                           --export=BIOProjectID="$BIOProjectID",Run="$run",geo_accession="$geo_accession",PRJdir="$PRJdir",Logdir="$Logdir" \
                           --output="$sample_log_dir/RSEM-%j.out" \
                           --error="$sample_log_dir/RSEM-%j.err" \
                           RSEM_v2.sh)
    fi

    job_id=$(echo "$job_output" | awk '{print $4}')
    if [[ "$job_id" =~ ^[0-9]+$ ]]; then
        echo "  Submitted RSEM job with job ID: $job_id"
        job_ids+=("$job_id")
    else
        echo "  Failed to get job ID from sbatch output: $job_output"
    fi
done
#######################
# Matrix and QC files #
#######################

echo "Generating combined matrices for all samples and final MultiQC report."

if [ ${#job_ids[@]} -eq 0 ]; then
    echo "Error: No RSEM job IDs collected. Cannot submit matrix generation."
    exit 1
fi

dep_string=$(IFS=:; echo "${job_ids[*]}")

matrix_job_id=$(sbatch \
    --dependency=afterok:$dep_string \
    --export=baseDir="$baseDir",PRJdir="$PRJdir",Logdir="$Logdir",AccList="$AccList",BIOProjectID="$BIOProjectID" \
    --output="$Logdir/RSEMmatrix-%j.out" \
    --error="$Logdir/RSEMmatrix-%j.err" \
    RSEMmatrix_v2.sh | awk '{print $4}')

if [ -z "$matrix_job_id" ]; then
    echo "Error running RSEMmatrix_v1.sh"
    exit 1
fi

###################################
# Generate bed and JSON files #
###################################

# Ensure log file is defined
log_file="${Logdir}/json_bed_submission.log"
touch "$log_file"

# Ensure matrix job ID is valid before proceeding
if [ -z "$matrix_job_id" ]; then
    echo "Matrix job ID not found. Cannot proceed to BED/JSON submission." >> "$log_file"
    exit 1
fi

wait_for_jobs "$matrix_job_id"

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
            "2bigWigjson") script_name="BWjson_v2.sh" ;;
            "geneJson")    script_name="gJSON_v2.sh" ;;
            "txBedJson")   script_name="tJSON_v2.sh" ;;
            "GeneTPMbed")  script_name="GeneTPMbed_v2.sh" ;;
            "txTPMbed")    script_name="txTPMbed_v2.sh" ;;
        esac

        if [ ! -x "$script_name" ]; then
            echo "Script $script_name not found or not executable. Skipping $job_type for $geo_accession" >> "$log_file"
            continue
        fi

        job_id=$(sbatch \
            --dependency=afterok:$matrix_job_id \
            --job-name="${job_type}_$geo_accession" \
            --export=baseDir="$baseDir",BIOProjectID="$BIOProjectID",Run="$run",geo_accession="$geo_accession",PRJdir="$PRJdir",tissue="$Tissue",strain="$Strain",sex="$Sex",title="${Title}",Sample_characteristics="${Sample_characteristics}",PMID="$PMID",GEOpath="$GEOpath",unique_name="$unique_name",Logdir="$Logdir",scratch_dir="$scratch_dir",StrainInfo="${StrainInfo}" \
            --output="$sample_log_dir/${job_type}-%j.out" \
            --error="$sample_log_dir/${job_type}-%j.err" \
            "$script_name" | awk '{print $4}')

        json_bed_job_ids["$geo_accession,$job_type"]=$job_id

        if [ -z "$job_id" ]; then
            echo "Error submitting $job_type job for $geo_accession" >> "$log_file"
        else
            echo "Submitted $job_type for $geo_accession with job ID: $job_id" >> "$log_file"
        fi
    done

done < "$AccList"

echo "All JSON and BED generation jobs have been submitted." >> "$log_file"
