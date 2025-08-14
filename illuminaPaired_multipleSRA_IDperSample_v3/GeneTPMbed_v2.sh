#!/bin/sh
#SBATCH --job-name=geneBed
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50000kb
#SBATCH --time=01:00:00
#SBATCH --account=akwitek
#SBATCH --partition=normal
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wdemos@mcw.edu
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

module load htslib samtools/1.20

# === Environment variables ===
Run=${Run}
geo_accession=${geo_accession}
BIOProjectID=${BIOProjectID}
Tissue=${tissue}
Strain=${strain}
Sex=${sex}
unique_name=${unique_name}
PRJdir=${PRJdir}
scratch_dir=${scratch_dir}
Logdir=${Logdir}

echo "Logdir is set to: $Logdir"
echo "Run is ${Run}"
echo "geo_accession is ${geo_accession}"
echo "BIOProjectID is ${BIOProjectID}"
echo "Tissue is ${Tissue}"
echo "Strain is ${Strain}"
echo "Sex is ${Sex}"
echo "unique_name is ${unique_name}"
echo "scratch_dir is ${scratch_dir}"

export LC_ALL=C

cd "$scratch_dir" || { echo "Error: Failed to move into directory $scratch_dir"; exit 1; }
echo "Moved into directory: $scratch_dir"

mkdir -p "${PRJdir}/${geo_accession}/log_files/bedsNjsons"
log_file="${PRJdir}/${geo_accession}/log_files/bedsNjsons/${geo_accession}_GeneBedTPM.out"
exec 1>"$log_file" 2>&1

start_time=$(date +%s)

GBed="/home/wdemos/NCBI_1August2023/mod_genes_sorted.bed"

# === RGB coloring based on TPM ===
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

# === Gene expression processing ===
gene_file="${scratch_dir}/${geo_accession}/${geo_accession}.genes.results"

if [ -f "$gene_file" ]; then
    echo "Gene file exists. Proceeding to create JSON."

    output_file="${scratch_dir}/${geo_accession}/${geo_accession}_Genes_TPMs.txt"
    awk 'BEGIN {OFS="\t"} FNR > 1 {print $1, $6}' "$gene_file" > "$output_file" || {
        echo "Error: Failed to parse $gene_file"
        exit 1
    }

    if [ ! -s "$output_file" ]; then
        echo "Error: $output_file does not exist or is empty"
        exit 1
    fi

    sorted_file="${scratch_dir}/${geo_accession}/${geo_accession}_sorted_Genes_TPMs.txt"
    echo "Sorting $output_file"
    sort -k1,1 "$output_file" > "$sorted_file" 2>"${PRJdir}/${geo_accession}/log_files/bedsNjsons/${geo_accession}_gsort_error.log" || {
        echo "Error: Failed to sort $output_file"
        cat "${PRJdir}/${geo_accession}/log_files/bedsNjsons/${geo_accession}_gsort_error.log"
        exit 1
    }

    echo "Debug: First few lines of sorted file:" >> "$log_file"
    head -n 10 "$sorted_file" >> "$log_file"

    if [ -f "$GBed" ]; then
        merged_file="${scratch_dir}/${geo_accession}/${geo_accession}_merged_GENE.bed"
        echo "Merging $sorted_file with $GBed"
        awk 'BEGIN {OFS="\t"} NR==FNR {a[$1]=$2; next} $4 in a {print $1, $2, $3, $4, a[$4]}' "$sorted_file" "$GBed" > "$merged_file" || {
            echo "Error: Failed to merge $sorted_file with $GBed"
            exit 1
        }
    else
        echo "Error: $GBed does not exist"
        exit 1
    fi

    rawOP="${scratch_dir}/${geo_accession}/${geo_accession}_GENE_rawTPM.bed"
    echo "Generating raw output file: $rawOP"
    while IFS=$'\t' read -r chr start end name score; do
        rgb=$(get_rgb "$score")
        echo -e "$chr\t$start\t$end\t$name\t$score\t.\t$start\t$end\t$rgb"
    done < "$merged_file" > "$rawOP"

    raw_count=$(wc -l < "$rawOP")
    echo "Number of rows in raw output: $raw_count"
    [ "$raw_count" -eq 0 ] && {
        echo "Error: $rawOP has no data"
        exit 1
    }

    chr_filtered="${scratch_dir}/${geo_accession}/${geo_accession}_NWfiltered.bed"
    chr0_filtered="${scratch_dir}/${geo_accession}/${geo_accession}_GENE_TPMfiltchr0.bed"

    echo "Filtering $rawOP to remove NW_ and 0 TPM rows"
    awk '!/^NW_/ && /^chr/' "$rawOP" > "$chr_filtered"
    awk '$5 != "0.00"' "$chr_filtered" > "$chr0_filtered" || {
        echo "Error: Failed to filter $chr_filtered"
        exit 1
    }

    filtered_count=$(wc -l < "$chr0_filtered")
    echo "Filtered row count: $filtered_count"
    [ "$filtered_count" -eq 0 ] && {
        echo "Error: $chr0_filtered has no data after filtering"
        exit 1
    }

    echo "Debug: First few lines of filtered file:" >> "$log_file"
    head -n 10 "$chr0_filtered" >> "$log_file"

    CHRsorted_file="${PRJdir}/${geo_accession}/RNAseq_${unique_name}_geneTPMfinalOP.bed"
    echo "Sorting by chromosome position to $CHRsorted_file"
    sort -k1,1 -k2,2n -k3,3n "$chr0_filtered" > "$CHRsorted_file" 2>>"${PRJdir}/${geo_accession}/log_files/bedsNjsons/${geo_accession}_gCHRsort_error.log" || {
        echo "Error: Sorting failed"
        cat "${PRJdir}/${geo_accession}/log_files/bedsNjsons/${geo_accession}_gCHRsort_error.log"
        exit 1
    }

    if [ ! -s "$CHRsorted_file" ]; then
        echo "Error: Sorted file is empty or missing"
        exit 1
    fi

    echo "Compressing final BED file"
    bgzip "$CHRsorted_file" || {
        echo "Error: bgzip failed"
        exit 1
    }

    tabix -p bed "${CHRsorted_file}.gz" || {
        echo "Error: tabix failed"
        exit 1
    }

    echo "Gene TPM BED processing complete for $geo_accession"
    end_time=$(date +%s)
    runtime=$((end_time - start_time))
    echo "Processing time: $runtime seconds"

else
    echo "Error: Gene results file $gene_file does not exist"
    exit 1
fi
