#!/bin/bash
#SBATCH --job-name=STAR8
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=3gb
#SBATCH --time=05:00:00
#SBATCH --account=your-slurm-account
#SBATCH --partition=normal
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=your@email.edu

###############################################################################
# USER CONFIGURATION
# Update the variables below for your system before running.
###############################################################################
# Scratch filesystem root; scratch_dir is built as SCRATCH_BASE/BIOProjectID
SCRATCH_BASE="/path/to/your/scratch/mount"
# GRCr8 GTF annotation file
REF_GTF="/path/to/your/GRCr8/reference/updated_GCF_036323735.1_CM070413.1_GRCr8_genomic.gtf"
# GRCr8 genome FASTA file
GENOME_FASTA="/path/to/your/GRCr8/reference/updated_GCF_036323735.1_CM070413.1_GRCr8_genomic.fna"
###############################################################################

# NOTE: Output and error logs are unified via submission script using --output and --error

# Script details
## Written by Wendy Demos
## Written 16 August 2024
## Modified 12 Feb 2026 by WMD to utilize GRCr8 and run in single project or batch mode
## Also updated to handle AccList file that has expanded Sample characteristics field.
## Modified to properly handle multiple fastq files per sample
## Modified 17 Feb 2026 by WMD: added bamCoverage BigWig generation from sorted STAR genome BAM

# Read sample ID and corresponding reads from command-line arguments
# Arguments match what pairedEnd_multi_complete.bash passes:
geo_accession=$1      # Sample ID (e.g., GSM12345)
READ1_FILES=$2        # Comma-separated list of forward reads
READ2_FILES=$3        # Comma-separated list of reverse reads
BIOProjectID=$4       # Project ID
unique_name=$5        # Unique sample identifier

# Enable debugging by printing each command before execution
#set -x
module load star/2.7.10b samtools/1.20 deeptools/3.5.1

# Retrieve environment variables set by the controller script via --export
# These are passed from pairedEnd_multi_complete.bash
PRJdir=${PRJdir}      # e.g., /path/to/your/home/data/expression/GEO/PRJNA123/reads_fastq
Logdir=${Logdir}      # e.g., /path/to/your/home/data/expression/GEO/PRJNA123/log_files
baseDir=${baseDir}    # e.g., /path/to/your/home/data/expression/GEO/PRJNA123

# Define scratch and reference paths
INDEX_DIR="$scratch_dir/RefIndex"

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
    --quantMode TranscriptomeSAM \
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

echo "============================================"
echo "STAR Processing Complete!"
echo "============================================"
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