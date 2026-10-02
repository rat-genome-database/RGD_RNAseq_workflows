#!/usr/bin/env bash
#SBATCH --job-name=STAR8
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=3gb
#SBATCH --time=08:00:00
#SBATCH --account=your-slurm-account
#SBATCH --partition=normal
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=your@email.edu

###############################################################################
# USER CONFIGURATION
# Update the variables below for your system before running.
###############################################################################
# Scratch area that holds the per-project scratch folders (scratch_dir = SCRATCH_BASE/BIOProjectID)
SCRATCH_BASE="/path/to/your/scratch/mount"
# GRCr8 GTF annotation file
REF_GTF="/path/to/your/GRCr8/reference/updated_GCF_036323735.1_CM070413.1_GRCr8_genomic.gtf"
# GRCr8 genome FASTA file
GENOME_FASTA="/path/to/your/GRCr8/reference/updated_GCF_036323735.1_CM070413.1_GRCr8_genomic.fna"
# Picard CollectRnaSeqMetrics inputs: rRNA interval list and refFlat annotation built from the same GTF
RRNA_INTERVALS="/path/to/your/GRCr8/reference/GRCr8_rRNA.intervals"
REF_FLAT="/path/to/your/GRCr8/reference/updated_GCF_036323735.1_CM070413.1_GRCr8_genomic.gtf.refflat"
###############################################################################

# NOTE: Output and error logs are unified via submission script using --output and --error

# Script details
## Written by Wendy Demos
## Written 16 August 2024
## Modified 12 Feb 2026 by WMD to utilize GRCr8 and run in single project or batch mode
## Also updated to handle AccList file that has expanded Sample characteristics field.
## Modified to properly handle multiple fastq files per sample
## Modified 17 Feb 2026 by WMD: added bamCoverage BigWig generation from sorted STAR genome BAM
## Modified 10 Sep 2026: strand detection added
## Modified 10 Sep 2026: added Picard CollectRnaSeqMetrics call; added polyA vs rRNA-depletion call from its metrics
## Modified 30 Sep 2026 (STAR_bigwig3): library prep call moved to libprep_call.sh (rule LP2, NonPolyA_perM,
##   calibrated on 24 studies). Picard intronic/mRNA thresholds (LP1) retired; Picard still runs and is reported.

# Read sample ID and corresponding reads from command-line arguments
# Arguments match what run_RNApipeline_pairedG8_diskGuard.bash passes:
geo_accession=$1      # Sample ID (e.g., GSM12345)
READ1_FILES=$2        # Comma-separated list of forward reads
READ2_FILES=$3        # Comma-separated list of reverse reads
BIOProjectID=$4       # Project ID
unique_name=$5        # Unique sample identifier

# Enable debugging by printing each command before execution
#set -x
module load star/2.7.10b samtools/1.20 deeptools/3.5.1 picard/2.25.0

# Retrieve environment variables set by the controller script via --export
# These are passed from run_RNApipeline_pairedG8_diskGuard.bash
PRJdir=${PRJdir}      # e.g., /path/to/your/home/data/expression/GEO/PRJNA123/reads_fastq
Logdir=${Logdir}      # e.g., /path/to/your/home/data/expression/GEO/PRJNA123/log_files
baseDir=${baseDir}    # e.g., /path/to/your/home/data/expression/GEO/PRJNA123

# Define scratch and reference paths
scratch_dir="${SCRATCH_BASE}/${BIOProjectID}"
INDEX_DIR="$scratch_dir/RefIndex"
# Library prep call script (shared with libprep_recall.sh). sbatch runs a spooled copy of this script,
# so default to the directory the job was submitted from; override with --export=LIBPREP_CALL=...
LIBPREP_CALL=${LIBPREP_CALL:-${SLURM_SUBMIT_DIR:-$(dirname "$0")}/libprep_call.sh}

# Final output directory for this sample
FinalOPdir="${PRJdir}/${geo_accession}"
mkdir -p "$FinalOPdir"

# Define output paths in scratch directory
OUTPUT_DIR="$scratch_dir/$geo_accession"
mkdir -p "$OUTPUT_DIR"

OUTPUT_PREFIX="${OUTPUT_DIR}/${geo_accession}_STAR"
LOG_OUTPUT="${OUTPUT_PREFIX}Log.final.out"
GENOME_BAM="${OUTPUT_PREFIX}Aligned.out.bam"
TRANSCRIPTOME_BAM="${OUTPUT_PREFIX}Aligned.toTranscriptome.out.bam"
SORT_GENOME_BAM="${OUTPUT_DIR}/${geo_accession}_GENOME_SORT.bam"
bigwig_output="${FinalOPdir}/RNAseq_${unique_name}.bigwig"
GENE_COUNTS_FILE="${FinalOPdir}/log_files/STAR/${geo_accession}_STARReadsPerGene.out.tab"
STRAND_DETECTION_LOG_FILE="${FinalOPdir}/log_files/STAR/${geo_accession}_strand_detection.log"
RNA_METRICS_FILE="${FinalOPdir}/log_files/STAR/${geo_accession}_rna_metrics.txt"
LIBRARY_PREP_LOG_FILE="${FinalOPdir}/log_files/STAR/${geo_accession}_library_prep.log"

echo "============================================"
echo "STAR Alignment Job Started: $(date "+%Y-%m-%d %H:%M:%S")"
echo "============================================"
echo "Sample: $geo_accession"
echo "Unique name: $unique_name"
echo "BioProject ID: $BIOProjectID"
echo "Scratch directory: $scratch_dir"
echo "Final output directory: $FinalOPdir"
echo "GTF reference file: $REF_GTF"
echo "Genome FASTA file: $GENOME_FASTA"
echo "STAR INDEX directory: $INDEX_DIR"
echo ""

# Verify INDEX directory exists
if [ ! -d "$INDEX_DIR" ]; then
    echo "ERROR: STAR index directory does not exist: $INDEX_DIR"
    echo "Please ensure starRef_v4.sh has completed successfully."
    exit 1
fi

# Verify essential index files exist
if [ ! -f "$INDEX_DIR/SA" ] || [ ! -f "$INDEX_DIR/SAindex" ] || [ ! -f "$INDEX_DIR/Genome" ]; then
    echo "ERROR: STAR index files are incomplete in $INDEX_DIR"
    echo "Please check starRef_v4.sh output and regenerate the index."
    exit 1
fi

# Remove trailing commas from input file lists
READ1_FILES="${READ1_FILES%,}"
READ2_FILES="${READ2_FILES%,}"

# Validate input files
if [[ -z "$READ1_FILES" || -z "$READ2_FILES" ]]; then
    echo "ERROR: No valid fastq files provided for sample $geo_accession"
    echo "READ1_FILES: '$READ1_FILES'"
    echo "READ2_FILES: '$READ2_FILES'"
    exit 1
fi

echo "Input FASTQ files:"
echo "READ1: $READ1_FILES"
echo "READ2: $READ2_FILES"
echo ""

# Verify all input files exist
IFS=',' read -ra R1_ARRAY <<< "$READ1_FILES"
IFS=',' read -ra R2_ARRAY <<< "$READ2_FILES"

if [ ${#R1_ARRAY[@]} -ne ${#R2_ARRAY[@]} ]; then
    echo "ERROR: Mismatch in number of R1 and R2 files"
    echo "R1 files: ${#R1_ARRAY[@]}"
    echo "R2 files: ${#R2_ARRAY[@]}"
    exit 1
fi

echo "Verifying ${#R1_ARRAY[@]} pairs of FASTQ files..."
missing_files=0
for i in "${!R1_ARRAY[@]}"; do
    if [ ! -f "${R1_ARRAY[$i]}" ]; then
        echo "ERROR: R1 file not found: ${R1_ARRAY[$i]}"
        missing_files=$((missing_files + 1))
    fi
    if [ ! -f "${R2_ARRAY[$i]}" ]; then
        echo "ERROR: R2 file not found: ${R2_ARRAY[$i]}"
        missing_files=$((missing_files + 1))
    fi
done

if [ $missing_files -gt 0 ]; then
    echo "ERROR: $missing_files file(s) missing. Aborting."
    exit 1
fi

echo "All input files verified successfully."
echo ""

# Check if all expected outputs already exist — skip if so
if [ -f "$LOG_OUTPUT" ] && [ -f "$SORT_GENOME_BAM" ] && [ -f "$bigwig_output" ]; then
    echo "All output files already exist for sample ${geo_accession}:"
    echo "  - $LOG_OUTPUT"
    echo "  - $SORT_GENOME_BAM"
    echo "  - $bigwig_output"
    echo "Skipping STAR alignment and BigWig generation (already processed)."
    exit 0
fi

echo "Output files will be written to: $OUTPUT_DIR"
echo "Output prefix: $OUTPUT_PREFIX"
echo ""

echo "============================================"
echo "Running STAR alignment..."
echo "============================================"

STAR --runThreadN 10 \
    --readFilesCommand zcat \
    --genomeDir "$INDEX_DIR" \
    --outFileNamePrefix "$OUTPUT_PREFIX" \
    --quantMode TranscriptomeSAM GeneCounts \
    --outSAMtype BAM Unsorted \
    --outSAMunmapped Within KeepPairs \
    --outFilterMatchNminOverLread 0.66 \
    --outFilterMultimapNmax 6 \
    --readFilesIn "$READ1_FILES" "$READ2_FILES"

if [ $? -ne 0 ]; then
    echo "ERROR: STAR alignment failed for sample ${geo_accession}"
    echo "Check the STAR log files in: $OUTPUT_DIR"
    exit 1
fi

echo ""
echo "STAR alignment completed successfully."
echo ""

# Verify output files were created
if [ ! -f "$GENOME_BAM" ]; then
    echo "ERROR: Expected genome BAM file not created: $GENOME_BAM"
    exit 1
fi

if [ ! -f "$TRANSCRIPTOME_BAM" ]; then
    echo "WARNING: Transcriptome BAM file not created: $TRANSCRIPTOME_BAM"
    echo "This may affect downstream RSEM analysis."
fi

echo "============================================"
echo "Sorting and indexing genome BAM file..."
echo "============================================"

samtools sort -@ 10 "$GENOME_BAM" -o "$SORT_GENOME_BAM"

if [ $? -ne 0 ]; then
    echo "ERROR: BAM sorting failed for sample ${geo_accession}"
    exit 1
fi

samtools index -b "$SORT_GENOME_BAM"

if [ $? -ne 0 ]; then
    echo "ERROR: BAM indexing failed for sample ${geo_accession}"
    exit 1
fi

echo "Sorted BAM created: $SORT_GENOME_BAM"
echo "BAM index created: ${SORT_GENOME_BAM}.bai"
echo ""

if [ ! -f "$SORT_GENOME_BAM" ]; then
    echo "ERROR: Sorted genome BAM file was not created for sample ${geo_accession}"
    exit 1
fi

echo "============================================"
echo "Copying log files to final output directory..."
echo "============================================"

mkdir -p "$FinalOPdir/log_files/STARQC"

if [ -f "${OUTPUT_PREFIX}Log.final.out" ]; then
    cp "${OUTPUT_PREFIX}Log.final.out" "$FinalOPdir/log_files/STARQC"
    echo "Copied: ${OUTPUT_PREFIX}Log.final.out"
fi

if [ -f "${OUTPUT_PREFIX}Log.out" ]; then
    cp "${OUTPUT_PREFIX}Log.out" "$FinalOPdir/log_files/STAR"
    echo "Copied: ${OUTPUT_PREFIX}Log.out"
fi
if  [ -f "${OUTPUT_PREFIX}ReadsPerGene.out.tab" ]; then
    cp "${OUTPUT_PREFIX}ReadsPerGene.out.tab" "$FinalOPdir/log_files/STAR"
    echo "Copied: ${OUTPUT_PREFIX}ReadsPerGene.out.tab"
fi

echo ""
echo "============================================"
echo "Generating BigWig file from sorted STAR genome BAM..."
echo "Normalized read coverage track (BPM, unique mappers only)"
echo "============================================"

bamCoverage \
    -b "$SORT_GENOME_BAM" \
    -o "$bigwig_output" \
    -p ${SLURM_CPUS_ON_NODE} \
    --normalizeUsing BPM \
    --binSize 10 \
    --minMappingQuality 255

if [ $? -ne 0 ]; then
    echo "ERROR: BigWig file generation failed for sample ${geo_accession}"
    exit 1
fi

echo "BigWig file generated successfully: $bigwig_output"
echo "File size: $(du -h "$bigwig_output" | cut -f1)"
echo ""

echo "================================================"
echo "STAR Processing and bigWig generation Complete!"
echo "================================================"
echo "Sample: $geo_accession"
echo "Unique name: $unique_name"
echo "Processing completed at: $(date "+%Y-%m-%d %H:%M:%S")"
echo ""
echo "Output files in scratch:"
echo "  - Genome BAM (unsorted): $GENOME_BAM"
echo "  - Transcriptome BAM:     $TRANSCRIPTOME_BAM"
echo "  - Sorted genome BAM:     $SORT_GENOME_BAM"
echo "  - STAR log:              $LOG_OUTPUT"
echo ""
echo "Output files in final directory:"
echo "  - BigWig:                $bigwig_output"
echo "  - Log files:             $FinalOPdir/log_files/"
echo "============================================"

echo "========================================================================="
echo "Calculating Strandedness and Library Prep Method - PolyA or rRNA depleted"
echo "========================================================================="

#Script take in the GENE_COUNTS_FILE generated by STAR GeneCounts method
#Determines strandedness of reads to pass to picard for inferring library prep
#Logic to justify  3x threshold:
#1. Antisense Transcription is not clean. Cells naturally express functional antisense RNAs. Real biological antisense 
#transcription would cause automation logic of 100 percent to 0% to falsely conclude that the library is unstranded.
#2. Kit Inefficiencies If a kit degrades slightly or experiences reagent issues, it might drop to 85% or 90% accuracy. The 3x threshold ensures 
# samples aren't misclassified to dub-optimal condition of a kit
#3. Math Breakdown of the Ratio: 
#Unstranded: Expected split is 50% / 50%. Neither side is anywhere near 3x the other. The pipeline safely outputs NONE.
#Stranded: A split of 75% / 25% is exactly a 3x ratio. Even if a sample suffers from severe degradation or high background noise, 
# a 75/25 split is still heavily biased toward one direction and represents a stranded library.
#Stranded: A normal sample split is usually 95% / 5%, which is a 19x ratio
#Calculate strandesness to pass off to picard:

# GENE_COUNTS_FILE and STRAND_DETECTION_LOG_FILE are defined at the top of the script
mkdir -p "$(dirname "$STRAND_DETECTION_LOG_FILE")"

if [ ! -s "$GENE_COUNTS_FILE" ]; then
    echo "ERROR: STAR gene counts file not found or empty: $GENE_COUNTS_FILE"
    echo "Cannot determine strandedness for sample ${geo_accession}"
    exit 1
fi

# 1. Sum up columns 2, 3, and 4 (skipping the header lines starting with N_)
sums=$(awk 'NR > 4 {col2+=$2; col3+=$3; col4+=$4} END {print col2, col3, col4}' "$GENE_COUNTS_FILE")

# 2. Extract individual totals
unstranded_count=$(echo "$sums" | cut -d' ' -f1)
forward_count=$(echo "$sums" | cut -d' ' -f2)
reverse_count=$(echo "$sums" | cut -d' ' -f3)

# Initialize log file
echo "=== STRAND DETECTION REPORT ===" > "$STRAND_DETECTION_LOG_FILE"
echo "Timestamp: $(date)" >> "$STRAND_DETECTION_LOG_FILE"
echo "Unstranded Counts (Col 2): $unstranded_count" >> "$STRAND_DETECTION_LOG_FILE"
echo "Forward Counts    (Col 3): $forward_count" >> "$STRAND_DETECTION_LOG_FILE"
echo "Reverse Counts    (Col 4): $reverse_count" >> "$STRAND_DETECTION_LOG_FILE"

# 3. Handle Edge Case: Empty or ultra-low input
if (( forward_count + reverse_count < 10000 )); then
    PICARD_STRAND="NONE"
    echo "WARNING: Total gene counts are dangerously low (<10,000 reads). Defaulting to NONE." >> "$STRAND_DETECTION_LOG_FILE"
    echo "STATUS: AMBIGUOUS_LOW_COUNTS" >> "$STRAND_DETECTION_LOG_FILE"

# 4. Determine strandedness and flag intermediate gray zones
elif (( forward_count > reverse_count )); then
    # Calculate ratio using awk to handle floating point math safely
    ratio=$(awk -v f="$forward_count" -v r="$reverse_count" 'BEGIN {print (r>0 ? f/r : f)}')
    
    if (( forward_count > reverse_count * 3 )); then
        PICARD_STRAND="FIRST_READ_TRANSCRIPTION_STRAND"
        echo "STATUS: SUCCESS (Forward biased, ratio: ${ratio}x)" >> "$STRAND_DETECTION_LOG_FILE"
    elif (( forward_count > reverse_count * 15 / 10 )); then # 1.5x threshold
        PICARD_STRAND="FIRST_READ_TRANSCRIPTION_STRAND"
        echo "WARNING: Weak forward strand bias detected (ratio: ${ratio}x). Library may be degraded or poorly prepared." >> "$STRAND_DETECTION_LOG_FILE"
        echo "STATUS: WARNING_WEAK_STRAND" >> "$STRAND_DETECTION_LOG_FILE"
    else
        PICARD_STRAND="NONE"
        echo "STATUS: SUCCESS (Classified as Unstranded, ratio: ${ratio}x)" >> "$STRAND_DETECTION_LOG_FILE"
    fi

else
    ratio=$(awk -v f="$forward_count" -v r="$reverse_count" 'BEGIN {print (f>0 ? r/f : r)}')
    
    if (( reverse_count > forward_count * 3 )); then
        PICARD_STRAND="SECOND_READ_TRANSCRIPTION_STRAND"
        echo "STATUS: SUCCESS (Reverse biased, ratio: ${ratio}x)" >> "$STRAND_DETECTION_LOG_FILE"
    elif (( reverse_count > forward_count * 15 / 10 )); then # 1.5x threshold
        PICARD_STRAND="SECOND_READ_TRANSCRIPTION_STRAND"
        echo "WARNING: Weak reverse strand bias detected (ratio: ${ratio}x). Library may be degraded or poorly prepared." >> "$STRAND_DETECTION_LOG_FILE"
        echo "STATUS: WARNING_WEAK_STRAND" >> "$STRAND_DETECTION_LOG_FILE"
    else
        PICARD_STRAND="NONE"
        echo "STATUS: SUCCESS (Classified as Unstranded, ratio: ${ratio}x)" >> "$STRAND_DETECTION_LOG_FILE"
    fi
fi

echo "Pipeline auto-detected strand specification: $PICARD_STRAND"
echo "Strand detection report: $STRAND_DETECTION_LOG_FILE"

#Picard CollectRnaSeqMetrics reports the fraction of aligned bases in mRNA (coding+UTR), intronic, intergenic and rRNA regions.
#These are reported as supporting evidence. The call itself is made by libprep_call.sh from non-polyadenylated small-RNA reads.

for ref_file in "$REF_FLAT" "$RRNA_INTERVALS" "$LIBPREP_CALL"; do
    if [ ! -s "$ref_file" ]; then
        echo "ERROR: Picard reference file not found or empty: $ref_file"
        exit 1
    fi
done

echo "Running Picard CollectRnaSeqMetrics (STRAND=$PICARD_STRAND)..."
java -jar "$PICARD" CollectRnaSeqMetrics \
     I="$SORT_GENOME_BAM" \
     O="$RNA_METRICS_FILE" \
     REF_FLAT="$REF_FLAT" \
     RIBOSOMAL_INTERVALS="$RRNA_INTERVALS" \
     STRAND="$PICARD_STRAND" \
     VALIDATION_STRINGENCY=SILENT

if [ $? -ne 0 ]; then
    echo "ERROR: Picard CollectRnaSeqMetrics failed for sample ${geo_accession}"
    exit 1
fi

# Library prep call (writes $LIBRARY_PREP_LOG_FILE; reads the ReadsPerGene, Picard and STAR logs in $FinalOPdir)
bash "$LIBPREP_CALL" "$geo_accession" "$FinalOPdir"
if [ $? -ne 0 ]; then
    echo "ERROR: library prep call failed for sample ${geo_accession}"
    exit 1
fi
echo "Library prep report: $LIBRARY_PREP_LOG_FILE"
	 
echo "========================================================================="
echo "Calculating Strandedness and Library Prep Method - COMPLETE"
echo "========================================================================="
