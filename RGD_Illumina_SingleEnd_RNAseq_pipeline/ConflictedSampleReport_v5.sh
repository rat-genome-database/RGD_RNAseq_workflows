#!/usr/bin/env bash

###############################################################################
# USER CONFIGURATION
# Update the variables below for your system before running.
###############################################################################
# Your home/base directory
myDir="/path/to/home"
###############################################################################


# Usage: ./script.sh <BIOProjectID>
if [ $# -ne 1 ]; then
    echo "Usage: $0 <BIOProjectID>"
    exit 1
fi

# Written by Wendy Demos
# Updated 7 April 2026 (v5): read sex result and TPM matrix from baseDir,
#                            write final conflict report to baseDir,


BIOProjectID="$1"

# Expect PRJdir/baseDir from caller environment
PRJdir="${PRJdir:-}"
baseDir="${baseDir:-}"

if [ -z "$baseDir" ]; then
    echo "ERROR: baseDir environment variable is not set."
    exit 1
fi

output_file="${baseDir}/${BIOProjectID}_sex_conflict_report.txt"
sex_result_file="${baseDir}/${BIOProjectID}_sex_result.txt"
tpm_matrix_file="${baseDir}/${BIOProjectID}.genes.TPM.matrix"
tmp_parsed="${baseDir}/tpm_parsed.txt"

if [ ! -f "$sex_result_file" ]; then
    echo "ERROR: sex result file not found: $sex_result_file"
    exit 1
fi

if [ ! -f "$tpm_matrix_file" ]; then
    echo "ERROR: TPM matrix file not found: $tpm_matrix_file"
    exit 1
fi

# Genes of interest
genes=("Xist" "Uty" "Sry" "Ddx3y" "Kdm5d" "Eif2s3y")

# Add header and note to the output file
note="Note: Female samples should have a high TPM for Xist and males high TPM for Uty, Sry, Ddx3y, Kdm5d, and Eif2s3y."
{
    echo "$note"
    echo -e "SampleID\tInputSex\tComputedSex\tXYRatio\tAgreement\t$(printf "%s\t" "${genes[@]}" | sed 's/\t$//')"
} > "$output_file"

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
    gsub(/"/, "", gene)
    if (gene in gene_map) {
        for (sample in sample_to_col) {
            tpm_values[sample][gene] = $(sample_to_col[sample])
        }
    }
}
END {
    for (sample in tpm_values) {
        for (gene in tpm_values[sample]) {
            print sample, gene, tpm_values[sample][gene]
        }
    }
}
' "$tpm_matrix_file" > "$tmp_parsed"

# Process the sex result file and generate the output
while read -r sample input_sex computed_sex ratio agreement; do
    # Skip header
    if [ "$sample" = "SampleID" ]; then
        continue
    fi

    # Match sample with parsed TPM values
    match=$(grep -w "$sample" "$tmp_parsed" || true)
    if [ -n "$match" ]; then
        # Extract TPM values for each gene of interest
        tpm_values=()
        for gene in "${genes[@]}"; do
            value=$(echo "$match" | awk -v gene="$gene" '$2 == gene {print $3}')
            tpm_values+=("${value:-NA}")
        done

        # Append the result to the output file
        echo -e "$sample\t$input_sex\t$computed_sex\t$ratio\t$agreement\t$(printf "%s\t" "${tpm_values[@]}" | sed 's/\t$//')" >> "$output_file"
    else
        # Still write line even if no TPM values matched
        echo -e "$sample\t$input_sex\t$computed_sex\t$ratio\t$agreement\tNA\tNA\tNA\tNA\tNA\tNA" >> "$output_file"
    fi
done < "$sex_result_file"

# Clean up temporary files
rm -f "$tmp_parsed"

echo "Output written to $output_file"