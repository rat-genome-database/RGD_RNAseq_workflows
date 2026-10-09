#!/usr/bin/env bash
# Per-sample library prep call: POLYA / RRNA_DEPLETED / AMBIGUOUS / UNKNOWN.
# Single source of the call: STAR_SE_v2.sh runs it after Picard for new samples, and
# libprep_recall.sh runs it on already-processed projects. Same inputs -> same call.
#
# Written 30 Sep 2026 by WMD
#
# Usage: libprep_call.sh <GSM> <sample_dir>
#   sample_dir = .../reads_fastq/<GSM>; reads from it:
#     log_files/STAR/<GSM>_STARReadsPerGene.out.tab   (STAR GeneCounts)          required
#     log_files/STAR/<GSM>_rna_metrics.txt            (Picard CollectRnaSeqMetrics) required
#     log_files/STARQC/<GSM>_STARLog.final.out        (STAR mapping rates)       optional, QC flag only
#   writes log_files/STAR/<GSM>_library_prep.log (same keys as before plus the new ones below)
#   env NONPOLYA_GENES overrides the gene list (default: GRCr8 list from build_nonpolyA_genes.sh)
#
# Rule LP2: NonPolyA_perM = reads on non-polyadenylated small-RNA genes per million gene-assigned reads.
# Picard intronic/mRNA are reported as supporting evidence only.

RULE_VERSION="LP2 (2026-09-30)"
POLYA_MAX=350          # POLYA if NonPolyA_perM <= this (calibration polyA max 299)
RRNA_DEPL_MIN=700      # RRNA_DEPLETED if >= this (calibration rRNA-depleted min 865)
MIN_ASSIGNED=500000    # UNKNOWN below this many gene-assigned reads
LOWQ_UNIQ=50      # QC flag only: STAR unique % below this ...
LOWQ_MULTI=30     # ... and multimapped % above this
# Gene list from build_nonpolyA_genes.sh: edit the default for your system, or export NONPOLYA_GENES
NONPOLYA_GENES=${NONPOLYA_GENES:-/path/to/your/GRCr8/reference/GRCr8_nonpolyA_genes.txt}

if [ $# -ne 2 ]; then
    echo "Usage: $0 <GSM> <sample_dir>"
    exit 1
fi
sample=$1
sdir=${2%/}
rpg=$sdir/log_files/STAR/${sample}_STARReadsPerGene.out.tab
metrics=$sdir/log_files/STAR/${sample}_rna_metrics.txt
star_log=$sdir/log_files/STARQC/${sample}_STARLog.final.out
out=$sdir/log_files/STAR/${sample}_library_prep.log

for f in "$NONPOLYA_GENES" "$rpg" "$metrics"; do
    if [ ! -s "$f" ]; then
        echo "ERROR $sample: required input missing or empty: $f"
        exit 1
    fi
done

# Strand Picard actually ran with, from the command line it records in the metrics header.
# Older runs: fall back to the "Strand used" line of the existing report.
strand=$(grep -m1 -oE 'STRAND(_SPECIFICITY)?=[A-Z_]+' "$metrics" | cut -d= -f2)
if [ -z "$strand" ] && [ -f "$out" ]; then
    strand=$(awk -F': *' '$1 == "Strand used" {print $2}' "$out")
fi
case "$strand" in
    FIRST_READ_TRANSCRIPTION_STRAND)  col=3 ;;   # STAR col 3 = htseq-count -s yes
    SECOND_READ_TRANSCRIPTION_STRAND) col=4 ;;   # STAR col 4 = htseq-count -s reverse
    NONE)                             col=2 ;;
    *) echo "ERROR $sample: cannot determine library strand from $metrics"; exit 1 ;;
esac

# Picard percentages, located by column name
read -r pct_ribosomal pct_mrna pct_intronic pct_intergenic < <(awk -F'\t' '
    /^## METRICS CLASS/ {
        getline; for (i = 1; i <= NF; i++) col[$i] = i
        getline
        n = split("PCT_RIBOSOMAL_BASES PCT_MRNA_BASES PCT_INTRONIC_BASES PCT_INTERGENIC_BASES", want, " ")
        for (k = 1; k <= n; k++) {
            v = (want[k] in col) ? $col[want[k]] : ""
            printf "%s%s", (v == "" ? "NA" : v), (k < n ? " " : "\n")
        }
        exit
    }' "$metrics")

# Non-polyA counts per million gene-assigned reads (STAR's first 4 lines are N_ summaries)
read -r assigned nonpolyA_reads < <(awk -v c="$col" '
    NR == FNR { if ($0 !~ /^#/ && NF) g[$1] = 1; next }
    FNR > 4   { t += $c; if ($1 in g) n += $c }
    END       { printf "%d %d\n", t, n }' "$NONPOLYA_GENES" "$rpg")
nonpolyA_perM=$(awk -v n="$nonpolyA_reads" -v t="$assigned" 'BEGIN { if (t > 0) printf "%.0f", n / t * 1e6; else print "NA" }')

# STAR mapping rates -> QC flag (does not change the call: failed-depletion samples were still called correctly)
uniq_pct=NA; multi_pct=NA
if [ -f "$star_log" ]; then
    read -r uniq_pct multi_pct < <(awk -F'|' '
        function v(s) { gsub(/[ \t%]/, "", s); return s }
        /Uniquely mapped reads %/            {u = v($2)}
        /% of reads mapped to multiple loci/ {m = v($2)}
        END {print (u == "" ? "NA" : u), (m == "" ? "NA" : m)}' "$star_log")
fi
qc=$(awk -v u="$uniq_pct" -v m="$multi_pct" -v lu="$LOWQ_UNIQ" -v lm="$LOWQ_MULTI" 'BEGIN {
    if (u == "NA" || m == "NA") print "NA"
    else if (u + 0 < lu + 0 && m + 0 > lm + 0) print "LOW_MAPPING"
    else print "PASS" }')

# The call and a plain-language reason
read -r call < <(awk -v v="$nonpolyA_perM" -v t="$assigned" -v mn="$MIN_ASSIGNED" -v pa="$POLYA_MAX" -v rd="$RRNA_DEPL_MIN" 'BEGIN {
    if (v == "NA" || t + 0 < mn + 0) print "UNKNOWN"
    else if (v + 0 >= rd + 0)        print "RRNA_DEPLETED"
    else if (v + 0 <= pa + 0)        print "POLYA"
    else                             print "AMBIGUOUS" }')
case "$call" in
    RRNA_DEPLETED) reason="NonPolyA_perM $nonpolyA_perM >= $RRNA_DEPL_MIN: non-polyadenylated small RNAs retained, as rRNA depletion does and oligo-dT selection does not" ;;
    POLYA)         reason="NonPolyA_perM $nonpolyA_perM <= $POLYA_MAX: non-polyadenylated small RNAs removed, as oligo-dT (polyA) selection does" ;;
    AMBIGUOUS)     reason="NonPolyA_perM $nonpolyA_perM is between $POLYA_MAX and $RRNA_DEPL_MIN, outside the range of every calibration study; check the study methods" ;;
    UNKNOWN)       reason="only $assigned gene-assigned reads (< $MIN_ASSIGNED) - too few for a reliable estimate" ;;
esac

# Supporting evidence: does Picard intronic agree? (0.15 / 0.25 = the LP1 rule this replaces)
support=$(awk -v i="$pct_intronic" -v c="$call" 'BEGIN {
    if (i == "NA") { print "NA"; exit }
    p = (i + 0 < 0.15) ? "POLYA" : (i + 0 >= 0.25 ? "RRNA_DEPLETED" : "AMBIGUOUS")
    if (c != "POLYA" && c != "RRNA_DEPLETED") print "intronic suggests " p
    else if (p == c) print "AGREES (intronic " i ")"
    else print "DISAGREES (intronic " i " suggests " p "; intronic varies by tissue and prep and is not used for the call)" }')

tmp=$(mktemp "${out}.XXXX")
{
    echo "=== LIBRARY PREP REPORT ==="
    echo "Timestamp: $(date)"
    echo "Rule version:         $RULE_VERSION"
    echo "Picard metrics file:  $metrics"
    echo "Strand used:          $strand"
    echo "PCT_MRNA_BASES:       $pct_mrna"
    echo "PCT_INTRONIC_BASES:   $pct_intronic"
    echo "PCT_INTERGENIC_BASES: $pct_intergenic"
    echo "PCT_RIBOSOMAL_BASES:  $pct_ribosomal"
    echo "GENE_ASSIGNED_READS:  $assigned"
    echo "NONPOLYA_READS:       $nonpolyA_reads"
    echo "NONPOLYA_PER_M:       $nonpolyA_perM"
    echo "NONPOLYA_GENE_LIST:   $NONPOLYA_GENES"
    echo "NONPOLYA_LIST_MD5:    $(md5sum "$NONPOLYA_GENES" | cut -d' ' -f1)"
    echo "STAR_UNIQUE_PCT:      $uniq_pct"
    echo "STAR_MULTI_PCT:       $multi_pct"
    echo "Thresholds: POLYA if NonPolyA_perM <= $POLYA_MAX; RRNA_DEPLETED if >= $RRNA_DEPL_MIN; AMBIGUOUS between; UNKNOWN if < $MIN_ASSIGNED gene-assigned reads"
    echo "LIBRARY_PREP: $call"
    echo "LIBRARY_PREP_REASON: $reason"
    echo "LIBRARY_PREP_SUPPORT: $support"
    echo "QC_FLAG: $qc"
} > "$tmp" && mv "$tmp" "$out"

echo "$sample: $call ($reason) [QC $qc]"
