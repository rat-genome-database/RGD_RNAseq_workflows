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
# Updated 31 Aug 2026 (v6): the sex-gene TPM ratio now makes the final call and
#     overrides the X:Y idxstats call from ComputeSex_v6.sh (unchanged, still
#     runs pre-RSEM, still reported in the XYRatio column).
#         Ymed = median(Uty, Ddx3y, Kdm5d, Eif2s3y)  # Sry excluded, 0.00 in males
#         R    = (Xist + 1) / (Ymed + 1);  R >= 10 -> F,  R <= 1 -> M
#     Ambiguous or no signal -> ComputedSex=Undetermined and Agreement=Conflict,
#     so a sample needing review cannot read as a clean call.
#     The call is written back to the conflict report AND to _sex_result.txt,
#     because BWjson_v7.sh reads column 3 of the latter. The original X:Y calls
#     are kept in ${BIOProjectID}_sex_result.xy.txt. Both formats are unchanged.
# Updated 6 Oct 2026 (v7): strand and library prep columns and the study summary (workflow 2.1.0).
#     Per sample, read from
#     ${baseDir}/reads_fastq/<GSM>/log_files/STAR/<GSM>_library_prep.log (written by libprep_call.sh,
#     run by STAR_SE_v2.sh) and log_files/STARQC/<GSM>_STARLog.final.out; NA if missing. Appended columns:
#       ComputedStrand ComputedLibraryPrep PctMRNA PctIntronic PctIntergenic PctRibosomal
#       UniqPct MultiPct TooShortPct NonPolyA_perM LibraryPrepRule LibraryPrepQC LibraryPrepReason
#     Also writes ${BIOProjectID}_library_prep_summary.txt (StudyCall = a call held by >= 80% of samples
#     with a POLYA/RRNA_DEPLETED call, else MIXED). Sex logic and _sex_result.txt unchanged.


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
libprep_summary="${baseDir}/${BIOProjectID}_library_prep_summary.txt"

# Pipeline version, written at the start of line 1, before the note (update with each release)
WORKFLOW_VERSION="HPC RGD single-end workflow 2.1.0"

# v6 call thresholds
FEMALE_R=10     # R >= this -> F
MALE_R=1        # R <= this -> M
MIN_SIGNAL=1    # max(Xist,Ymed) below this TPM -> no gene call, Undetermined

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

# v6: keep the original X:Y calls once, for rollback
sex_result_backup="${baseDir}/${BIOProjectID}_sex_result.xy.txt"
if [ ! -f "$sex_result_backup" ]; then
    cp "$sex_result_file" "$sex_result_backup"
fi

# Add header and note to the output file
note="Workflow: ${WORKFLOW_VERSION}. Note: ComputedSex is called from sex-gene expression: R = (Xist TPM + 1) / (median TPM of Uty, Ddx3y, Kdm5d and Eif2s3y + 1); R >= ${FEMALE_R} is called F, R <= ${MALE_R} is called M, and other values, or samples where Xist and the Y-gene median are both below ${MIN_SIGNAL} TPM, are Undetermined. Females express Xist and males the Y-linked genes. Sry is reported but not used in the calculation. XYRatio (X:Y read coverage) is shown for reference only. Agreement compares ComputedSex with InputSex (submitted to GEO)."
{
    echo "$note"
#    echo -e "SampleID\tInputSex\tComputedSex\tXYRatio\tAgreement\t$(printf "%s\t" "${genes[@]}" | sed 's/\t$//')"
#v7: strand, library prep and calibration metric columns appended
    echo -e "SampleID\tInputSex\tComputedSex\tXYRatio\tAgreement\t$(printf "%s\t" "${genes[@]}" | sed 's/\t$//')\tComputedStrand\tComputedLibraryPrep\tPctMRNA\tPctIntronic\tPctIntergenic\tPctRibosomal\tUniqPct\tMultiPct\tTooShortPct\tNonPolyA_perM\tLibraryPrepRule\tLibraryPrepQC\tLibraryPrepReason"
} > "$output_file"

# v6: rebuilt _sex_result.txt, same 5 columns ComputeSex_v6.sh writes
updated_result="${baseDir}/${BIOProjectID}_sex_result.tmp"
echo -e "SampleID\tInputSex\tComputedSex\tRatio\tAgreement" > "$updated_result"

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

n_override=0
n_fallback=0
n_legacy=0

# Process the sex result file and generate the output
while read -r sample input_sex computed_sex ratio agreement; do
    # Skip header
    if [ "$sample" = "SampleID" ]; then
        continue
    fi

    # Match sample with parsed TPM values
    match=$(grep -w "$sample" "$tmp_parsed" || true)

    # v7: strand and library prep from the per-sample report written by libprep_call.sh (run by STAR_SE_v2.sh)
    prep_log="${baseDir}/reads_fastq/${sample}/log_files/STAR/${sample}_library_prep.log"
    computed_strand=""
    computed_library_prep=""
    pct_mrna=""; pct_intronic=""; pct_intergenic=""; pct_ribosomal=""
    prep_rule=""; prep_qc=""; prep_reason=""; nonpolyA_perM=""
    if [ -f "$prep_log" ]; then
        computed_strand=$(awk -F': *' '$1 == "Strand used" {print $2}' "$prep_log")
        computed_library_prep=$(awk -F': *' '$1 == "LIBRARY_PREP" {print $2}' "$prep_log")
        prep_rule=$(awk -F': *' '$1 == "Rule version" {print $2}' "$prep_log")
        prep_qc=$(awk -F': *' '$1 == "QC_FLAG" {print $2}' "$prep_log")
        prep_reason=$(awk '/^LIBRARY_PREP_REASON:/ {sub(/^LIBRARY_PREP_REASON: */, ""); print}' "$prep_log")
        nonpolyA_perM=$(awk -F': *' '$1 == "NONPOLYA_PER_M" {print $2}' "$prep_log")
        if [ -z "$prep_rule" ]; then
            prep_rule="LP1 (legacy)"
            prep_reason="made by the retired intronic rule; re-call with libprep_recall.sh"
            n_legacy=$((n_legacy + 1))
        fi
        read -r pct_mrna pct_intronic pct_intergenic pct_ribosomal <<< "$(awk -F': *' '
            $1 == "PCT_MRNA_BASES"       {m = $2}
            $1 == "PCT_INTRONIC_BASES"   {i = $2}
            $1 == "PCT_INTERGENIC_BASES" {g = $2}
            $1 == "PCT_RIBOSOMAL_BASES"  {r = $2}
            END {print (m == "" ? "NA" : m), (i == "" ? "NA" : i), (g == "" ? "NA" : g), (r == "" ? "NA" : r)}' "$prep_log")"
    else
        echo "  WARNING $sample: no library prep report at $prep_log -> ComputedStrand/ComputedLibraryPrep = NA"
    fi
    computed_strand=${computed_strand:-NA}
    computed_library_prep=${computed_library_prep:-NA}
    pct_mrna=${pct_mrna:-NA}; pct_intronic=${pct_intronic:-NA}
    pct_intergenic=${pct_intergenic:-NA}; pct_ribosomal=${pct_ribosomal:-NA}
    prep_rule=${prep_rule:-NA}; prep_qc=${prep_qc:-NA}; prep_reason=${prep_reason:-NA}
    nonpolyA_perM=${nonpolyA_perM:-NA}

    # v7: STAR mapping rates from the per-sample Log.final.out copied by STAR_SE_v2.sh
    star_log="${baseDir}/reads_fastq/${sample}/log_files/STARQC/${sample}_STARLog.final.out"
    uniq_pct="NA"; multi_pct="NA"; short_pct="NA"
    if [ -f "$star_log" ]; then
        read -r uniq_pct multi_pct short_pct <<< "$(awk -F'|' '
            function v(s) { gsub(/[ \t%]/, "", s); return s }
            /Uniquely mapped reads %/            {u = v($2)}
            /% of reads mapped to multiple loci/ {m = v($2)}
            /% of reads unmapped: too short/     {s = v($2)}
            END {print (u == "" ? "NA" : u), (m == "" ? "NA" : m), (s == "" ? "NA" : s)}' "$star_log")"
    else
        echo "  WARNING $sample: no STAR Log.final.out at $star_log -> mapping rates = NA"
    fi
    prep_cols="$computed_strand\t$computed_library_prep\t$pct_mrna\t$pct_intronic\t$pct_intergenic\t$pct_ribosomal\t$uniq_pct\t$multi_pct\t$short_pct\t$nonpolyA_perM\t$prep_rule\t$prep_qc\t$prep_reason"

    if [ -n "$match" ]; then
        # Extract TPM values for each gene of interest
        tpm_values=()
        for gene in "${genes[@]}"; do
            value=$(echo "$match" | awk -v gene="$gene" '$2 == gene {print $3}')
            tpm_values+=("${value:-NA}")
        done

        # v6: call sex from the gene ratio. Median (not sum) of the Y genes so
        # one dropped gene cannot carry the call; Eif2s3y alone runs ~10x Uty.
        # A missing gene comes through as NA, which awk reads as 0.
        gene_call=$(awk -v x="${tpm_values[0]}" -v y1="${tpm_values[1]}" \
                        -v y2="${tpm_values[3]}" -v y3="${tpm_values[4]}" \
                        -v y4="${tpm_values[5]}" \
                        -v fr="$FEMALE_R" -v mr="$MALE_R" -v ms="$MIN_SIGNAL" '
        BEGIN {
            y[1]=y1+0; y[2]=y2+0; y[3]=y3+0; y[4]=y4+0
            for (i=1;i<=4;i++) for (j=i+1;j<=4;j++) if (y[j]<y[i]) { t=y[i]; y[i]=y[j]; y[j]=t }
            ymed = (y[2]+y[3])/2
            R    = (x+1)/(ymed+1)
            call = "NA"
            if (x+0 >= ms || ymed >= ms) { if (R >= fr) call="F"; else if (R <= mr) call="M" }
            printf "%s %.4f\n", call, R
        }')
        gene_sex="${gene_call%% *}"
        gene_R="${gene_call#* }"

        if [ "$gene_sex" != "NA" ]; then
            if [ "$gene_sex" != "$computed_sex" ]; then
                echo "  $sample: X:Y=$ratio said $computed_sex, sex genes (R=$gene_R) say $gene_sex"
                n_override=$((n_override + 1))
            fi
            computed_sex="$gene_sex"

            # v6: agreement is recomputed against the overridden call
            if [ "$input_sex" = "$computed_sex" ]; then
                agreement="Agree"
            else
                agreement="Conflict"
            fi
        else
            # v6: no decisive gene call - flag BOTH columns so the row cannot
            # read as clean. This previously fell back to the X:Y call and could
            # emit "M / Agree" for a sample expressing Xist AND the Y genes, which
            # means contamination, undeclared pooling, or an aneuploid animal.
            echo "  WARNING $sample: no decisive sex-gene signal (R=$gene_R, X:Y said $computed_sex) -> Undetermined/Conflict"
            computed_sex="Undetermined"
            agreement="Conflict"
            n_fallback=$((n_fallback + 1))
        fi

        # Append the result to the output file
#        echo -e "$sample\t$input_sex\t$computed_sex\t$ratio\t$agreement\t$(printf "%s\t" "${tpm_values[@]}" | sed 's/\t$//')" >> "$output_file"
#v7: strand, library prep and calibration metrics appended
        echo -e "$sample\t$input_sex\t$computed_sex\t$ratio\t$agreement\t$(printf "%s\t" "${tpm_values[@]}" | sed 's/\t$//')\t$prep_cols" >> "$output_file"
    else
        # Still write line even if no TPM values matched
#        echo -e "$sample\t$input_sex\t$computed_sex\t$ratio\t$agreement\tNA\tNA\tNA\tNA\tNA\tNA" >> "$output_file"
        echo -e "$sample\t$input_sex\t$computed_sex\t$ratio\t$agreement\tNA\tNA\tNA\tNA\tNA\tNA\t$prep_cols" >> "$output_file"
    fi

    # v6: every row is written here, matched or not, so BWjson can still find it
    printf "%s\t%s\t%s\t%s\t%s\n" "$sample" "$input_sex" "$computed_sex" "$ratio" "$agreement" >> "$updated_result"
done < "$sex_result_file"

# v6: swap the corrected calls into place for BWjson_v7.sh
mv "$updated_result" "$sex_result_file"

# v7: per-study library prep summary.
# Columns are located by header name so this survives future column additions.
awk -F'\t' -v study="$BIOProjectID" '
function median(arr, n,    i, j, t) {
    if (n == 0) return "NA"
    for (i = 1; i <= n; i++) for (j = i + 1; j <= n; j++) if (arr[j] < arr[i]) { t = arr[i]; arr[i] = arr[j]; arr[j] = t }
    return (n % 2) ? arr[(n + 1) / 2] : (arr[n / 2] + arr[n / 2 + 1]) / 2
}
NR == 1 { next }                                   # note line
NR == 2 { for (i = 1; i <= NF; i++) col[$i] = i; next }
{
    n_samples++
    calls[$col["ComputedLibraryPrep"]]++
    rules[$col["LibraryPrepRule"]]++
    split("PctMRNA PctIntronic PctIntergenic PctRibosomal UniqPct MultiPct TooShortPct NonPolyA_perM", keys, " ")
    for (k = 1; k <= 8; k++) {
        v = $col[keys[k]]
        if (v != "NA" && v != "") { cnt[keys[k]]++; vals[keys[k], cnt[keys[k]]] = v + 0 }
    }
}
END {
    printf "Study\tSamples\tPOLYA\tRRNA_DEPLETED\tAMBIGUOUS\tOther"
    split("PctMRNA PctIntronic PctIntergenic PctRibosomal UniqPct MultiPct TooShortPct NonPolyA_perM", keys, " ")
    for (k = 1; k <= 8; k++) printf "\tmedian_%s", keys[k]
    printf "\tStudyCall\tStudyCallBasis\tRules\n"
    other = n_samples - calls["POLYA"] - calls["RRNA_DEPLETED"] - calls["AMBIGUOUS"]
    printf "%s\t%d\t%d\t%d\t%d\t%d", study, n_samples, calls["POLYA"], calls["RRNA_DEPLETED"], calls["AMBIGUOUS"], other
    for (k = 1; k <= 8; k++) {
        n = cnt[keys[k]] + 0
        delete a
        for (i = 1; i <= n; i++) a[i] = vals[keys[k], i]
        printf "\t%s", median(a, n)
    }
    # study call = the call held by >= 80% of samples with a definite call
    p = calls["POLYA"] + 0; r = calls["RRNA_DEPLETED"] + 0; d = p + r
    if (d == 0)               { sc = "UNKNOWN"; basis = "no sample has a POLYA or RRNA_DEPLETED call" }
    else if (p >= 0.8 * d)    { sc = "POLYA"; basis = p " of " d " definite sample calls" }
    else if (r >= 0.8 * d)    { sc = "RRNA_DEPLETED"; basis = r " of " d " definite sample calls" }
    else                      { sc = "MIXED"; basis = p " POLYA, " r " RRNA_DEPLETED - check for mixed preps or sample mix-ups" }
    rl = ""; for (x in rules) rl = rl (rl == "" ? "" : "; ") x " x" rules[x]
    printf "\t%s\t%s\t%s\n", sc, basis, rl
}' "$output_file" > "$libprep_summary"

# Clean up temporary files
rm -f "$tmp_parsed"

echo "Output written to $output_file"
echo "Library prep summary written to $libprep_summary"
echo "Sex calls updated in $sex_result_file"
echo "  $n_override overridden by sex genes, $n_fallback Undetermined"
if [ "$n_legacy" -gt 0 ]; then
    echo "  WARNING: $n_legacy library prep calls were made by the retired LP1 rule - run libprep_recall.sh $BIOProjectID, then rerun this script"
fi
