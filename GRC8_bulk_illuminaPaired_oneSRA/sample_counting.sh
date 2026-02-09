#!/usr/bin/env bash
###############################################################################
# sample_counting.sh
# Test script to verify sample counting from AccList files
###############################################################################

# Example AccList file path (replace with your actual file)
ACCLIST_FILE="$1"

if [[ -z "$ACCLIST_FILE" ]]; then
    echo "Usage: $0 <path_to_acclist_file>"
    echo ""
    echo "Example:"
    echo "  $0 /home/wdemos/STAR_RSEM_pipeline/AccessionLists/Completed/GSE201312_test_AccList.txt"
    exit 1
fi

if [[ ! -f "$ACCLIST_FILE" ]]; then
    echo "ERROR: File not found: $ACCLIST_FILE"
    exit 1
fi

echo "=================================================="
echo "Testing Sample Count for: $ACCLIST_FILE"
echo "=================================================="
echo ""

# Show first few lines of the file
echo "First 5 lines of file:"
echo "--------------------------------------------------"
head -5 "$ACCLIST_FILE"
echo "--------------------------------------------------"
echo ""

# Count total lines (excluding header)
total_lines=$(awk 'NR > 1 && !/^#/ && NF > 0' "$ACCLIST_FILE" | wc -l)
echo "Total data lines (excluding header): $total_lines"
echo ""

# Show unique GSM IDs
echo "Unique GSM IDs (geo_accession in column 2):"
echo "--------------------------------------------------"
awk -F'\t' 'NR > 1 && !/^#/ && $2 != "" {print $2}' "$ACCLIST_FILE" | sort -u
echo "--------------------------------------------------"
echo ""

# Count unique GSM IDs
unique_count=$(awk -F'\t' '
    NR > 1 && !/^#/ && $2 != "" {
        gsm[$2] = 1
    }
    END {
        print length(gsm)
    }
' "$ACCLIST_FILE")

echo "=================================================="
echo "RESULT: $unique_count unique samples (GSM IDs)"
echo "=================================================="
echo ""

# Show distribution
echo "Sample distribution (how many runs per GSM):"
awk -F'\t' 'NR > 1 && !/^#/ && $2 != "" {count[$2]++} END {for (gsm in count) print gsm ": " count[gsm] " runs"}' "$ACCLIST_FILE" | sort

echo ""
echo "This project would be classified as:"
if [[ $unique_count -le 20 ]]; then
    echo " SMALL ($unique_count samples ≤ 20)"
else
    echo " LARGE ($unique_count samples > 20)"
fi