#!/bin/sh
#SBATCH --job-name=TxBed
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50000kb
#SBATCH --time=01:00:00
#SBATCH --account=akwitek
#SBATCH --partition=normal
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wdemos@mcw.edu
# #SBATCH --output=%x-%j.out 
# #SBATCH --error=TxBed.err

module load htslib samtools/1.20

# retrieve environment variables
Run=${Run}
geo_accession=${geo_accession}
BIOProjectID="$BIOProjectID" # This should be the GEO Project Name
Tissue="$tissue"
Strain="$strain"
Sex="$sex"
unique_name=${unique_name}
PRJdir=${PRJdir} # $myDir/data/expression/GEO/$BIOProjectID/reads_fastq
scratch_dir=${scratch_dir} # /scratch/g/akwitek/wdemos/$BIOProjectID
Logdir=${Logdir} # $myDir/data/expression/GEO/$BIOProjectID/log_files
# FinalOPdir=${Logdir}

echo "Logdir is set to: $Logdir"
echo "Run is ${Run}"
echo "geo_accession is ${geo_accession}"
echo "BIOProjectID is $BIOProjectID"
echo "Tissue is ${Tissue}"
echo "strain is ${Strain}"
echo "Sex is ${Sex}"
echo "unique_name is ${unique_name}"
echo "scratch_dir is ${scratch_dir}"

export LC_ALL=C

cd "$scratch_dir" || { echo "Error: Failed to move into directory $scratch_dir"; exit 1; }
echo "Moved into directory: $scratch_dir"

# Ensure the log directory exists
mkdir -p "${PRJdir}/${geo_accession}/log_files/bedsNjsons"

# Redirect output to the log file in the project directory after it is created
log_file="${PRJdir}/${geo_accession}/log_files/bedsNjsons/${geo_accession}_TxBedTPM.out"
exec 1>"$log_file" 2>&1

# Capture start time
start_time=$(date +%s)

TBed=/home/wdemos/NCBI_1August2023/mod_transcripts_sorted.bed

# Assigns RGB color codes based on TPM values
get_rgb() {
    local value=$1
    if [ "$(echo "$value <= 0.5" | bc -l)" -eq 1 ]; then
        echo "128,128,128"
    elif [ "$(echo "$value > 0.5" | bc -l)" -eq 1 ] && [ "$(echo "$value <= 10" | bc -l)" -eq 1 ]; then
        echo "173,216,230"
    elif [ "$(echo "$value > 10" | bc -l)" -eq 1 ] && [ "$(echo "$value <= 1000" | bc -l)" -eq 1 ]; then
        echo "0,0,205"
    else
        echo "0,0,139"
    fi
}

# Check if the transcripts results input file exists
tx_file="${scratch_dir}/${geo_accession}/${geo_accession}.transcripts.results"
if [ -f "$tx_file" ]; then
    echo "transcripts file exists. Proceeding to create JSON."

    output_file="${scratch_dir}/${geo_accession}/${geo_accession}_TX_TPMs.txt"

    # Extract transcript IDs and TPMs from the transcripts results file
    echo "Parsing transcript name and TPMs from $tx_file"
    awk 'BEGIN {OFS="\t"} FNR > 1 {print $1, $6}' "$tx_file" > "$output_file"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to parse $tx_file" >> "$log_file"
        exit 1
    fi

    if [ -f "$output_file" ] && [ -s "$output_file" ]; then
        echo "File exists and is not empty: $output_file"
    else
        echo "Error: $output_file does not exist or is empty" >> "$log_file"
        exit 1
    fi

    # Sort the TPM output file by transcript ID
    sorted_file="${scratch_dir}/${geo_accession}/${geo_accession}_sorted_TX_TPMs.txt"
    echo "Sorting $output_file"
    sort -k1,1 "$output_file" > "$sorted_file" 2>"${PRJdir}/${geo_accession}/log_files/bedsNjsons/${geo_accession}_TXsort_error.log" 
    if [ $? -ne 0 ]; then
        echo "Error: Failed to sort $output_file" >> "$log_file"
        cat "${PRJdir}/${geo_accession}/log_files/bedsNjsons/${geo_accession}_TXsort_error.log"  >> "$log_file"
        exit 1
    fi

    # Debugging: Log the first few lines of sorted file
    echo "Debug: First few lines of sorted file ($sorted_file):" >> "$log_file"
    head -n 10 "$sorted_file" >> "$log_file"

    # Merge with the reference BED file by transcript name
    if [ -f "$TBed" ]; then
        merged_file="${scratch_dir}/${geo_accession}/${geo_accession}_merged_TX.bed"
        echo "Merging $sorted_file with $TBed"
        awk 'BEGIN {OFS="\t"} NR==FNR {a[$1]=$2; next} $4 in a {print $1, $2, $3, $4, a[$4]}' "$sorted_file" "$TBed" > "$merged_file"
        if [ $? -ne 0 ]; then
            echo "Error: Failed to merge $sorted_file with $TBed" >> "$log_file"
            exit 1
        fi
    else
        echo "Error: $TBed does not exist" >> "$log_file"
        exit 1
    fi

    # Generate the raw output file with all TPM values
    rawOP="${scratch_dir}/${geo_accession}/${geo_accession}_TX_rawTPM.bed"
    echo "Generating raw output file: $rawOP"
    while IFS=$'\t' read -r chr start end name score; do
        rgb=$(get_rgb "$score")
        echo -e "$chr\t$start\t$end\t$name\t$score\t.\t$start\t$end\t$rgb"
    done < "$merged_file" > "$rawOP"

    raw_count=$(wc -l < "$rawOP")
    echo "Number of rows in the raw final output file ($rawOP): $raw_count" >> "$log_file"
    if [ $raw_count -eq 0 ]; then
        echo "Error: $rawOP has no data." >> "$log_file"
        exit 1
    fi

    # Filter out rows starting with "NC_" or "NW_" and TPM equal to 0
    chr_filtered="${scratch_dir}/${geo_accession}/${geo_accession}_NWfiltered.bed"
    chr0_filtered="${scratch_dir}/${geo_accession}/${geo_accession}_TX_TPMfiltchr0.bed"

    echo "Filtering $rawOP to remove NW_ (unlocalized or unplaced scaffold)"
    awk '!/^NW_/ && /^chr/' "$rawOP" > "$chr_filtered"

    echo "Filtering $chr_filtered to remove rows where TPM is 0.00"
    awk '$5 != "0.00"' "$chr_filtered" > "$chr0_filtered"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to filter $chr_filtered" >> "$log_file"
        exit 1
    fi

    filtered_count=$(wc -l < "$chr0_filtered")
    echo "Number of rows in the filtered final output file ($chr0_filtered): $filtered_count" >> "$log_file"
    if [ $filtered_count -eq 0 ]; then
        echo "Error: $chr0_filtered has no data after filtering." >> "$log_file"
        exit 1
    fi

    # Debugging: Log the first few lines of the filtered file
    echo "Debug: First few lines of chr0_filtered before sorting:" >> "$log_file"
    head -n 10 "$chr0_filtered" >> "$log_file"

    # Sorting to create the final zipped bed file and tabix index
#    CHRsorted_file="${Logdir}/${geo_accession}/RNAseq_${unique_name}_TXTPMfinalOP.bed" ##try sending to scratch and only final to home
    CHRsorted_file="${PRJdir}/${geo_accession}/RNAseq_${unique_name}_TXTPMfinalOP.bed"
    echo "Sorting $chr0_filtered by chromosome position" >> "$log_file"
    sort -k1,1 -k2,2n -k3,3n "$chr0_filtered" > "$CHRsorted_file" 2>>"${PRJdir}/${geo_accession}/log_files/bedsNjsons/${geo_accession}_txCHRsort_error.log"
    # Check if sorting succeeded
    if [ $? -ne 0 ]; then
        echo "Error: Sort failed for $chr0_filtered" >> "$log_file"
        cat "${PRJdir}/${geo_accession}/log_files/bedsNjsons/${geo_accession}_txCHRsort_error.log" >> "$log_file"
        exit 1
    fi

    # Ensure that the sorted file exists and is non-empty before compressing
    if [ ! -s "$CHRsorted_file" ]; then
        echo "Error: Sorted file $CHRsorted_file is empty or does not exist." >> "$log_file"
        exit 1
    fi

    # Compress the sorted BED file using bgzip
    echo "compress $CHRsorted_file"
    bgzip "$CHRsorted_file"
    if [ $? -ne 0 ]; then
        echo "Error: bgzip failed to compress $CHRsorted_file" >> "$log_file"
        exit 1
    fi

    # Index the compressed BED file with tabix
    tabix -p bed "${CHRsorted_file}.gz"
    if [ $? -ne 0 ]; then
        echo "Error: tabix indexing failed for ${CHRsorted_file}.gz" >> "$log_file"
        exit 1
    fi

    # Final confirmation of successful completion
    echo "Transcript TPM BED file processing is complete for $geo_accession" >> "$log_file"
    end_time=$(date +%s)
    runtime=$((end_time - start_time))
    echo "Processing time: $runtime seconds" >> "$log_file"

else
    echo "Error: Transcript results file $tx_file does not exist" >> "$log_file"
    exit 1
fi
