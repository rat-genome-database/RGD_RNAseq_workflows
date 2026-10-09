#!/usr/bin/env bash
# Build the non-polyadenylated small-RNA gene list used by libprep_call.sh.
# Run ONCE per reference annotation; the list is a versioned reference file, not rebuilt per run,
# so every sample called against the same annotation uses the same list.
#
# Written 30 Sep 2026 by WMD
#
# Genes kept: GTF gene_biotype snoRNA, snRNA, ncRNA (scaRNAs on RefSeq), SRP_RNA, RNase_MRP_RNA, telomerase_RNA.
# These RNAs are not polyadenylated, so oligo-dT (polyA) selection removes them while rRNA depletion keeps them.
# misc_RNA is excluded: on GRCr8 it includes mis-annotated protein-coding transcripts.
#
# Usage: build_nonpolyA_genes.sh [GTF] [OUT]

# Defaults below: edit for your system, or pass GTF and OUT as arguments
REF_GTF=${1:-/path/to/your/GRCr8/reference/updated_GCF_036323735.1_CM070413.1_GRCr8_genomic.gtf}
OUT=${2:-/path/to/your/GRCr8/reference/GRCr8_nonpolyA_genes.txt}
BIOTYPES="snoRNA|snRNA|ncRNA|SRP_RNA|RNase_MRP_RNA|telomerase_RNA"

if [ ! -s "$REF_GTF" ]; then
    echo "ERROR: GTF not found or empty: $REF_GTF"
    exit 1
fi
if [ -e "$OUT" ]; then
    echo "ERROR: $OUT already exists. It is a versioned reference file; remove it deliberately or give a new name."
    exit 1
fi

tmp=$(mktemp)
awk -F'\t' '$3 == "gene"' "$REF_GTF" \
  | grep -E "gene_biotype \"($BIOTYPES)\"" \
  | grep -oE 'gene_id "[^"]+"' | cut -d'"' -f2 | sort -u > "$tmp"

n=$(wc -l < "$tmp")
if [ "$n" -eq 0 ]; then
    echo "ERROR: no genes matched biotypes $BIOTYPES in $REF_GTF"
    rm -f "$tmp"
    exit 1
fi

{
    echo "# Non-polyadenylated small-RNA gene list for library prep calls (libprep_call.sh)"
    echo "# Built:     $(date "+%Y-%m-%d %H:%M:%S") by $(whoami)"
    echo "# GTF:       $REF_GTF"
    echo "# GTF md5:   $(md5sum "$REF_GTF" | cut -d' ' -f1)"
    echo "# Biotypes:  ${BIOTYPES//|/, } (misc_RNA excluded)"
    echo "# Genes:     $n"
    cat "$tmp"
} > "$OUT"
rm -f "$tmp"
chmod a-w "$OUT"

echo "Wrote $OUT ($n genes, read-only)"
echo "List md5: $(md5sum "$OUT" | cut -d' ' -f1)"
