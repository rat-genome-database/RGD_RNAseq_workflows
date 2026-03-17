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

if [ ! -x ./STAR_v1.sh ]; then
  echo "ERROR: STAR_v1.sh not found or not executable. Please check the script path and permissions."
  exit 1
fi

wait_for_jobs() {
  local ids=("$@")
  local joined
  joined=$(IFS=,; echo "${ids[*]}")
  while squeue -h -j "$joined" >/dev/null 2>&1 && [ -n "$(squeue -h -j "$joined")" ]; do
    sleep 10
  done
}

star_job_ids=()

for geo_accession in "${geo_accessions[@]}"; do
  RUNS=$(awk -v sample="$geo_accession" -F'\t' '$2 == sample {print $1}' "$AccList")
  read -r Tissue Strain Sex <<< "$(awk -v sample="$geo_accession" -F'\t' '$2 == sample {print $3, $4, $5; exit}' "$AccList")"
  unique_name="${Tissue}_${Strain}_${Sex}_${geo_accession}"

  READ1_FILES=$(
    for run in $RUNS; do
      ls -1 "${scratch_dir}/${run}"/*.fastq.gz 2>/dev/null
    done | sort -V | paste -sd ","
  )
  READ1_FILES="${READ1_FILES%,}"

  echo "$(date): Processing Sample: $geo_accession"
  echo "  Runs: $RUNS"
  echo "  READ1: $READ1_FILES"

  if [[ -n "$READ1_FILES" ]]; then
    echo "  Submitting STAR alignment job..."
    sample_log_dir="$PRJdir/$geo_accession/log_files/STAR"
    mkdir -p "$sample_log_dir"

    sbatch_output=$(sbatch --parsable --job-name="STAR_$geo_accession" \
      --export=ALL,baseDir="$baseDir",PRJdir="$PRJdir",Logdir="$Logdir" \
      --output="$sample_log_dir/STAR-%j.out" \
      --error="$sample_log_dir/STAR-%j.err" \
      ${starRef_job_id:+--dependency=afterok:${starRef_job_id}} \
      STAR_v1.sh "$geo_accession" "$READ1_FILES" "$BIOProjectID" "$unique_name")

    echo "  sbatch output: $sbatch_output"
    job_id="$sbatch_output"

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
      echo "    Checked: ${scratch_dir}/${run}/*.fastq.gz"
    done
  fi
done

if [ ${#star_job_ids[@]} -eq 0 ]; then
  echo "ERROR: No valid STAR jobs submitted. Exiting."
  exit 1
fi

echo "Waiting for all STAR jobs to complete before proceeding..."
wait_for_jobs "${star_job_ids[@]}"
echo "All STAR jobs completed."



######################################
# Create Deduplicated Accession List #
######################################
STARQC_log_dir="$Logdir/STARQC"
mkdir -p "$STARQC_log_dir"
uniqueAccList="${STARQC_log_dir}/${BIOProjectID}_Unique_AccList.txt"
{
    head -n 1 "$AccList"
    tail -n +2 "$AccList" | sort -u -t$'\t' -k2,2
} > "$uniqueAccList"
echo "Unique accession list saved to $uniqueAccList"

##############################
# STAR Alignment Rate Check  #
##############################
echo "Check sample alignment rates for downstream processing."
STARQC_log_dir="$Logdir/STARQC"
mkdir -p "$STARQC_log_dir"

# Submit STARQC job
starqc_output=$(sbatch \
    --job-name="STARQC" \
    --export=Logdir="$Logdir",AccList="$uniqueAccList",BIOProjectID="$BIOProjectID",PRJdir="$PRJdir" \
    --output="$STARQC_log_dir/STARQC-%j.out" \
    --error="$STARQC_log_dir/STARQC-%j.err" \
    STARQC_v1.sh "$uniqueAccList" "$BIOProjectID" )

# Extract and validate job ID
starqc_job_id=$(echo "$starqc_output" | awk '{print $4}')
echo "STARQC_v1.sh job submission output: $starqc_output"

if [[ -z "$starqc_job_id" || ! "$starqc_job_id" =~ ^[0-9]+$ ]]; then
    echo "Error: Invalid STARQC job ID: '$starqc_job_id'"
    exit 1
fi

# Wait for STARQC to finish before consuming its output
wait_for_jobs "$starqc_job_id"
#echo "STARQC_v1.sh completed. Output expected at: ${Logdir}/${BIOProjectID}_STAR_Align_sum.txt"
echo "STARQC_v1.sh completed. Output expected at: ${STARQC_log_dir}/${BIOProjectID}_STAR_Align_sum.txt"

###########################################
# Filter to PASS samples for downstream   #
###########################################
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
     ($2 in pass)' OFS='\t' "$starqc_report" "$uniqueAccList" > "$passAccList"

# Quick stats
total_unique=$(($(wc -l < "$uniqueAccList") - 1))
total_pass=$(($(wc -l < "$passAccList") - 1))
total_fail=$(( total_unique - total_pass ))

echo "STARQC gating: unique=$total_unique  PASS=$total_pass  FAIL=$total_fail"
echo "PASS-only accession list saved to: $passAccList"

# If nothing passed, stop the pipeline early
if [ "$total_pass" -le 0 ]; then
    echo "ERROR: No samples passed alignment rate threshold. Exiting."
    exit 1
fi

##############################
#Sample Sex Estimation        #
##############################
echo "Computing sample sex."

echo "Check sample alignment rates for downstream processing."
ComputeSex_log_dir="$Logdir/ComputeSex"
mkdir -p "$ComputeSex_log_dir"

# Submit ComputeSex job
job_output=$(sbatch \
    --job-name="computeSex" \
    --export=PRJdir="$PRJdir",Logdir="$Logdir",AccList="$passAccList",BIOProjectID="$BIOProjectID",baseDir="$baseDir" \
    --output="$ComputeSex_log_dir/ComputeSex-%j.out" \
    --error="$ComputeSex_log_dir/ComputeSex-%j.err" \
    ComputeSex_v2.sh "$AccList" "$BIOProjectID")
#    ComputeSex_v2.sh "$uniqueAccList" "$BIOProjectID")

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
    rref_log_dir="$Logdir/RSEM"
    mkdir -p "$rref_log_dir"

if [ -f "$scratch_dir/rsemref.chrlist" ] && [ -f "$scratch_dir/rsemref.grp" ] && [ -f "$scratch_dir/rsemref.idx.fa" ] && [ -f "$scratch_dir/rsemref.n2g.idx.fa" ] && [ -f "$scratch_dir/rsemref.seq" ] && [ -f "$scratch_dir/rsemref.ti" ] && [ -f "$scratch_dir/rsemref.transcripts.fa" ]; then
    echo "RSEM reference already exists. Skipping RSEMref.sh."
else
    echo "Submitting RSEMref.sh to generate the RSEM reference."
    rsemRef_job_id=$(sbatch --export=BIOProjectID="$BIOProjectID" \
                           --output="$rref_log_dir/RSEMref-%j.out" \
                           --error="$rref_log_dir/RSEMref-%j.err" RSEMref_v2.sh | awk '{print $4}')
    if [ -z "$rsemRef_job_id" ]; then echo "Error: RSEMref_v2.sh submission failed."; exit 1; fi
    echo "RSEM reference job ID: $rsemRef_job_id"
fi

###########################
# RSEM Jobs for Samples   #
###########################
# We only submit samples that PASSED STARQC.
# $passAccList has the same columns as AccList (header + geo_accession in col 2).

# Sanity check
if [ ! -s "$passAccList" ]; then
    echo "ERROR: PASS-only accession list not found or empty: $passAccList"
    exit 1
fi

# Build unique GEO sample IDs from PASS list (skip header)
mapfile -t unique_geo_accessions < <(awk -F'\t' 'NR>1 && $2!="" {print $2}' "$passAccList" | sort -u)

job_ids=()

for geo_accession in "${unique_geo_accessions[@]}"; do
    # RSEM script now uses per-sample dirs like /scratch/.../$geo_accession):
     run="$geo_accession"

    sample_output_dir="$PRJdir/$geo_accession"
    sample_log_dir="$sample_output_dir/log_files/RSEM"
    mkdir -p "$sample_log_dir"

    echo "Submitting RSEM job for sample $geo_accession (run: ${run:-N/A})"

    if [ -z "$rsemRef_job_id" ]; then
        job_output=$(sbatch --job-name="${geo_accession}_RSEM" \
                           --export=BIOProjectID="$BIOProjectID",Run="$run",geo_accession="$geo_accession",PRJdir="$PRJdir",Logdir="$Logdir",baseDir="$baseDir" \
                           --output="$sample_log_dir/RSEM-%j.out" \
                           --error="$sample_log_dir/RSEM-%j.err" \
                           RSEM_v1.sh)
    else
        job_output=$(sbatch --dependency=afterok:$rsemRef_job_id \
                           --job-name="${geo_accession}_RSEM" \
                           --export=BIOProjectID="$BIOProjectID",Run="$run",geo_accession="$geo_accession",PRJdir="$PRJdir",Logdir="$Logdir",baseDir="$baseDir" \
                           --output="$sample_log_dir/RSEM-%j.out" \
                           --error="$sample_log_dir/RSEM-%j.err" \
                           RSEM_v1.sh)
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

echo "Generating combined matrices (PASS-only samples) and final MultiQC report."

echo "Check sample alignment rates for downstream processing."
Matrix_log_dir="$Logdir/RSEM"
mkdir -p "$Matrix_log_dir_log_dir"

# We should have RSEM job IDs only for PASS samples
if [ ${#job_ids[@]} -eq 0 ]; then
    echo "Error: No PASS samples (no RSEM job IDs collected). Cannot submit matrix generation."
    exit 1
fi

# Dependency string from PASS-only RSEM jobs
dep_string=$(IFS=:; echo "${job_ids[*]}")

# Use PASS list for downstream matrix script
matrix_job_id=$(sbatch \
    --dependency=afterok:$dep_string \
    --export=baseDir="$baseDir",PRJdir="$PRJdir",Logdir="$Logdir",AccList="$passAccList",BIOProjectID="$BIOProjectID" \
    --output="$Matrix_log_dir/RSEMmatrix-%j.out" \
    --error="$Matrix_log_dir/RSEMmatrix-%j.err" \
    RSEMmatrix_v2.sh | awk '{print $4}')

if [ -z "$matrix_job_id" ]; then
    echo "Error: RSEMmatrix_v2.sh submission failed."
    exit 1
fi

echo "RSEMmatrix_v2.sh submitted with job ID: $matrix_job_id (PASS-only)."


###########################################
# Generate BED and JSON files (PASS only) #
###########################################
echo "Generating BED/JSON files for PASS-only samples."

# Must have PASS RSEM jobs + a finished matrix job
if [ ${#job_ids[@]} -eq 0 ]; then
    echo "Error: No PASS samples (no RSEM job IDs collected). Cannot submit BED/JSON generation."
    exit 1
fi
if [[ -z "$matrix_job_id" || ! "$matrix_job_id" =~ ^[0-9]+$ ]]; then
    echo "Matrix job ID not found or invalid. Cannot proceed to BED/JSON submission."
    exit 1
fi
wait_for_jobs "$matrix_job_id"

bNj_log_dir="$Logdir/bedsNjsons"
 mkdir -p "$bNj_log_dir"

log_file="${bNj_log_dir}/json_bed_submission.log"
touch "$log_file"

declare -A seen_samples

while IFS=$'\t' read -r run geo_accession Tissue Strain Sex PMID GEOpath Title Sample_characteristics StrainInfo; do
    # Skip header & dupes
    [[ "$run" == "Run" ]] && continue
    [[ -n "${seen_samples[$geo_accession]}" ]] && continue
    seen_samples["$geo_accession"]=1

    unique_name="${Tissue}_${Strain}_${Sex}_${geo_accession}"
    sample_dir="$PRJdir/$geo_accession"
    sample_log_dir="$sample_dir/log_files/bedsNjsons"
    mkdir -p "$sample_log_dir"

    # Expected RSEM outputs
#    rsem_genes="$sample_dir/${geo_accession}.genes.results"
rsem_genes="$scratch_dir/${geo_accession}/${geo_accession}.genes.results"
#    rsem_tx="$sample_dir/${geo_accession}.transcripts.results"
    rsem_tx="$scratch_dir/${geo_accession}/${geo_accession}.transcripts.results"
    echo "Preparing BED/JSON submissions for $geo_accession" >> "$log_file"

    # Submit BWjson unconditionally (no BigWig check per request)
    sb_out=$(sbatch \
        --dependency=afterok:$matrix_job_id \
        --job-name="BWjson_${geo_accession}" \
        --export=baseDir="$baseDir",BIOProjectID="$BIOProjectID",Run="$run",geo_accession="$geo_accession",PRJdir="$PRJdir",tissue="$Tissue",strain="$Strain",sex="$Sex",title="$Title",Sample_characteristics="$Sample_characteristics",PMID="$PMID",GEOpath="$GEOpath",unique_name="$unique_name",Logdir="$Logdir",scratch_dir="$scratch_dir",StrainInfo="$StrainInfo" \
        --output="$sample_log_dir/BWjson-%j.out" \
        --error="$sample_log_dir/BWjson-%j.err" \
        BWjson_v2a.sh)
    bwjson_jid=$(echo "$sb_out" | awk '{print $4}')
    [[ "$bwjson_jid" =~ ^[0-9]+$ ]] || bwjson_jid=""
    [ -z "$bwjson_jid" ] && echo "  ERROR: BWjson submit failed: $sb_out" >> "$log_file" || echo "  BWjson JID: $bwjson_jid" >> "$log_file"

    # Gene BED -> gene JSON (guard on genes.results)
    if [ -s "$rsem_genes" ]; then
        sb_out=$(sbatch \
            --dependency=afterok:$matrix_job_id \
            --job-name="GeneTPMbed_${geo_accession}" \
            --export=baseDir="$baseDir",BIOProjectID="$BIOProjectID",Run="$run",geo_accession="$geo_accession",PRJdir="$PRJdir",tissue="$Tissue",strain="$Strain",sex="$Sex",title="$Title",Sample_characteristics="$Sample_characteristics",PMID="$PMID",GEOpath="$GEOpath",unique_name="$unique_name",Logdir="$Logdir",scratch_dir="$scratch_dir",StrainInfo="$StrainInfo" \
            --output="$sample_log_dir/GeneTPMbed-%j.out" \
            --error="$sample_log_dir/GeneTPMbed-%j.err" \
            GeneTPMbed_v2.sh)
        gene_bed_jid=$(echo "$sb_out" | awk '{print $4}')
        [[ "$gene_bed_jid" =~ ^[0-9]+$ ]] || gene_bed_jid=""
        if [ -z "$gene_bed_jid" ]; then
            echo "  ERROR: GeneTPMbed submit failed: $sb_out" >> "$log_file"
        else
            echo "  GeneTPMbed JID: $gene_bed_jid" >> "$log_file"
            # Chain geneJson to GeneTPMbed
            sb_out=$(sbatch \
                --dependency=afterok:$matrix_job_id:$gene_bed_jid \
                --job-name="geneJson_${geo_accession}" \
                --export=baseDir="$baseDir",BIOProjectID="$BIOProjectID",Run="$run",geo_accession="$geo_accession",PRJdir="$PRJdir",tissue="$Tissue",strain="$Strain",sex="$Sex",title="$Title",Sample_characteristics="$Sample_characteristics",PMID="$PMID",GEOpath="$GEOpath",unique_name="$unique_name",Logdir="$Logdir",scratch_dir="$scratch_dir",StrainInfo="$StrainInfo" \
                --output="$sample_log_dir/geneJson-%j.out" \
                --error="$sample_log_dir/geneJson-%j.err" \
                gJSON_v2.sh)
            gj_jid=$(echo "$sb_out" | awk '{print $4}')
            [[ "$gj_jid" =~ ^[0-9]+$ ]] || gj_jid=""
            [ -z "$gj_jid" ] && echo "  ERROR: geneJson submit failed: $sb_out" >> "$log_file" || echo "  geneJson JID: $gj_jid" >> "$log_file"
        fi
    else
        echo "  SKIP GeneTPMbed/geneJson: missing $rsem_genes" >> "$log_file"
    fi

    # Transcript BED -> txBedJson (guard on transcripts.results)
    if [ -s "$rsem_tx" ]; then
        sb_out=$(sbatch \
            --dependency=afterok:$matrix_job_id \
            --job-name="txTPMbed_${geo_accession}" \
            --export=baseDir="$baseDir",BIOProjectID="$BIOProjectID",Run="$run",geo_accession="$geo_accession",PRJdir="$PRJdir",tissue="$Tissue",strain="$Strain",sex="$Sex",title="$Title",Sample_characteristics="$Sample_characteristics",PMID="$PMID",GEOpath="$GEOpath",unique_name="$unique_name",Logdir="$Logdir",scratch_dir="$scratch_dir",StrainInfo="$StrainInfo" \
            --output="$sample_log_dir/txTPMbed-%j.out" \
            --error="$sample_log_dir/txTPMbed-%j.err" \
            txTPMbed_v2.sh)
        tx_bed_jid=$(echo "$sb_out" | awk '{print $4}')
        [[ "$tx_bed_jid" =~ ^[0-9]+$ ]] || tx_bed_jid=""
        if [ -z "$tx_bed_jid" ]; then
            echo "  ERROR: txTPMbed submit failed: $sb_out" >> "$log_file"
        else
            echo "  txTPMbed JID: $tx_bed_jid" >> "$log_file"
            # Chain txBedJson to txTPMbed
            sb_out=$(sbatch \
                --dependency=afterok:$matrix_job_id:$tx_bed_jid \
                --job-name="txBedJson_${geo_accession}" \
                --export=baseDir="$baseDir",BIOProjectID="$BIOProjectID",Run="$run",geo_accession="$geo_accession",PRJdir="$PRJdir",tissue="$Tissue",strain="$Strain",sex="$Sex",title="$Title",Sample_characteristics="$Sample_characteristics",PMID="$PMID",GEOpath="$GEOpath",unique_name="$unique_name",Logdir="$Logdir",scratch_dir="$scratch_dir",StrainInfo="$StrainInfo" \
                --output="$sample_log_dir/txBedJson-%j.out" \
                --error="$sample_log_dir/txBedJson-%j.err" \
                tJSON_v2.sh)
            tj_jid=$(echo "$sb_out" | awk '{print $4}')
            [[ "$tj_jid" =~ ^[0-9]+$ ]] || tj_jid=""
            [ -z "$tj_jid" ] && echo "  ERROR: txBedJson submit failed: $sb_out" >> "$log_file" || echo "  txBedJson JID: $tj_jid" >> "$log_file"
        fi
    else
        echo "  SKIP txTPMbed/txBedJson: missing $rsem_tx" >> "$log_file"
    fi

done < "$passAccList"

echo "All PASS-only BED/JSON submissions attempted; see $log_file for details."
