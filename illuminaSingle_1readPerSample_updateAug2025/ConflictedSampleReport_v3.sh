#!/bin/bash

# Usage: ./script.sh <path/sex_result_file> <path/tpm_matrix_file> <BIOProjectID>
#if [ $# -ne 3 ]; then
#    echo "Usage: $0 <sex_result_file> <tpm_matrix_file> <BIOProjectID>"
if [ $# -ne 1 ]; then
    echo "Usage: $0 <BIOProjectID>"
    exit 1
fi

#sex_result_file="$1"
#tpm_matrix_file="$2"
#BIOProjectID="$3"
BIOProjectID="$1"
output_file=${BIOProjectID}_sex_conflict_report.txt
PRJdir=/home/wdemos/data/expression/GEO/${BIOProjectID}
sex_result_file=${PRJdir}/${BIOProjectID}_sex_result.txt
tpm_matrix_file=${PRJdir}/${BIOProjectID}.genes.TPM.matrix


# Genes of interest
genes=("Xist" "Uty" "Sry" "Ddx3y" "Kdm5d" "Eif2s3y")

# Add header and note to the output file
note="Note: Female samples should have a high TPM for Xist and males high TPM for Uty, Sry, Ddx3y, Kdm5d, and Eif2s3y."
{
    echo "$note"
    echo -e "SampleID\tInputSex\tComputedSex\t${genes[*]}"
} > "${PRJdir}/$output_file"

# Parse the TPM matrix file
awk -v genes="${genes[*]}" '
BEGIN {
    split(genes, gene_list, " ")
    for (i in gene_list) {
        gene_map[gene_list[i]] = 1
    }
}
NR == 1 {
    # Parse header to map sample column indices
    for (i = 2; i <= NF; i++) {
        gsub(/\.genes\.results/, "", $i)
        sample_to_col[$i] = i
    }
    next
}
{
    # For each gene of interest, store TPM values by sample
    gene = $1
    gsub(/"/, "", gene)  # Remove quotes from gene name
    if (gene in gene_map) {
        for (sample in sample_to_col) {
            tpm_values[sample][gene] = $(sample_to_col[sample])
        }
    }
}
END {
    # Store TPM values for quick lookup
    for (sample in tpm_values) {
        for (gene in tpm_values[sample]) {
            print sample, gene, tpm_values[sample][gene]
        }
    }
}
' "$tpm_matrix_file" > ${PRJdir}/tpm_parsed.txt

# Process the sex result file and generate the output
while read -r sample input_sex computed_sex; do
    # Match sample with parsed TPM values
    match=$(grep -w "$sample" ${PRJdir}/tpm_parsed.txt)
    if [ -n "$match" ]; then
        # Extract TPM values for each gene of interest
        tpm_values=()
        for gene in "${genes[@]}"; do
            value=$(echo "$match" | awk -v gene="$gene" '$2 == gene {print $3}')
            tpm_values+=("$value")
        done

        # Append the result to the output file
        echo -e "$sample\t$input_sex\t$computed_sex\t${tpm_values[*]}" >> "${PRJdir}/$output_file"
    fi
done < "$sex_result_file"

# Parse the output file to keep only the header and lines with "Conflict" in the 5th column
#temp_file=$(mktemp)
#{
#    # Include the header
#    head -n 2 "$output_file"
#    # Include only lines with "Conflict" in the 5th column
#    awk -F"\t" '$5 == "Conflict"' "$output_file"
#} > "$temp_file"
#mv "$temp_file" "$output_file"

# Clean up temporary files
rm -f ${PRJdir}/tpm_parsed.txt

echo "Output written to ${PRJdir}/$output_file"

