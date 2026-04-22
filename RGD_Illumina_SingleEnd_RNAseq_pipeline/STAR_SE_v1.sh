#!/usr/bin/env bash
#SBATCH --job-name=STAR_SE
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=3gb
#SBATCH --time=05:00:00
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=normal
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=your@email.edu
# NOTE: Output and error logs are set via submission script using --output and --error

###############################################################################
# USER CONFIGURATION
# Update the variables below for your system before running.
###############################################################################
#path to base directory
myDir="/path/to/home"
# Scratch filesystem root; scratch_dir is built as SCRATCH_BASE/BIOProjectID
SCRATCH_BASE="/path/to/scratch"
# GRCr8 GTF annotation file
REF_GTF="$myDir/GCF_036323735.1_GRCr8/updated_GCF_036323735.1_CM070413.1_GRCr8_genomic.gtf"
# GRCr8 genome FASTA file
GENOME_FASTA="$myDir/GCF_036323735.1_GRCr8/updated_GCF_036323735.1_CM070413.1_GRCr8_genomic.fna"
#mRatBN7.2 for validation testing
#REF_GTF="/path/to/genome/NCBI_1August2023/mod_GCF_015227675.2_mRatBN7.2_genomic.gtf"
#GENOME_FASTA="/path/to/genome/NCBI_1August2023/mod_GCF_015227675.2_mRatBN7.2_genomic.fna"

###############################################################################

# Script details
## Single-end STAR alignment → sorted BAM → BigWig (BPM, unique mappers only)
## Arguments passed from run_RNApipeline_SE_diskGuard.bash:
##   $1  geo_accession   Sample ID (e.g., GSM12345)
##   $2  READ1_FILES     Comma-separated list of single-end FASTQ files
##   $3  BIOProjectID    Project ID
##   $4  unique_name     Unique sample identifier

geo_accession=$1
READ1_FILES=$2
BIOProjectID=$3
unique_name=$4

module load star/2.7.10b samtools/1.20 deeptools/3.5.1

# Retrieve environment variables set by the controller script via --export
PRJdir=${PRJdir}
Logdir=${Logdir}
baseDir=${baseDir}

# Define scratch and reference paths
scratch_dir="${SCRATCH_BASE}/${BIOProjectID}"
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
echo "STAR_SE Alignment Job Started: $(date "+%Y-%m-%d %H:%M:%S")"
echo "============================================"
echo "Sample:          $geo_accession"
echo "Unique name:     $unique_name"
echo "BioProject ID:   $BIOProjectID"
echo "Scratch dir:     $scratch_dir"
echo "Final output:    $FinalOPdir"
echo "GTF reference:   $REF_GTF"
echo "STAR index:      $INDEX_DIR"
echo ""

# Verify STAR index exists
if [ ! -d "$INDEX_DIR" ]; then
    echo "ERROR: STAR index directory does not exist: $INDEX_DIR"
    echo "Please ensure starRef_v4.sh has completed successfully."
    exit 1
fi
if [ ! -f "$INDEX_DIR/SA" ] || [ ! -f "$INDEX_DIR/SAindex" ] || [ ! -f "$INDEX_DIR/Genome" ]; then
    echo "ERROR: STAR index files are incomplete in $INDEX_DIR"
    exit 1
fi

# Remove trailing commas from input file list
READ1_FILES="${READ1_FILES%,}"

if [[ -z "$READ1_FILES" ]]; then
    echo "ERROR: No valid FASTQ files provided for sample $geo_accession"
    exit 1
fi

echo "Input FASTQ files: $READ1_FILES"
echo ""

# Verify all input files exist
IFS=',' read -ra R1_ARRAY <<< "$READ1_FILES"
missing_files=0
for f in "${R1_ARRAY[@]}"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: FASTQ file not found: $f"
        missing_files=$((missing_files + 1))
    fi
done
if [ $missing_files -gt 0 ]; then
    echo "ERROR: $missing_files file(s) missing. Aborting."
    exit 1
fi
echo "All input files verified (${#R1_ARRAY[@]} file(s))."
echo ""

# Skip if all expected outputs already exist
if [ -f "$LOG_OUTPUT" ] && [ -f "$SORT_GENOME_BAM" ] && [ -f "$bigwig_output" ]; then
    echo "All output files already exist for ${geo_accession} — skipping."
    exit 0
fi

echo "============================================"
echo "Running STAR alignment (single-end)..."
echo "============================================"

STAR --runThreadN 10 \
    --readFilesCommand zcat \
    --genomeDir "$INDEX_DIR" \
    --outFileNamePrefix "$OUTPUT_PREFIX" \
    --quantMode TranscriptomeSAM \
    --outSAMtype BAM Unsorted \
    --outFilterMatchNminOverLread 0.66 \
    --outFilterMultimapNmax 6 \
    --readFilesIn "$READ1_FILES"

if [ $? -ne 0 ]; then
    echo "ERROR: STAR alignment failed for sample ${geo_accession}"
    exit 1
fi
echo "STAR alignment completed successfully."
echo ""

# Verify output BAM was created
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

echo "Sorted BAM:  $SORT_GENOME_BAM"
echo "BAM index:   ${SORT_GENOME_BAM}.bai"
echo ""

echo "============================================"
echo "Copying log files to final output directory..."
echo "============================================"
mkdir -p "$FinalOPdir/log_files/STARQC"
mkdir -p "$FinalOPdir/log_files/STAR"

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
echo "Generating BigWig from sorted STAR genome BAM..."
echo "Normalized read coverage (BPM, unique mappers only)"
echo "============================================"

bamCoverage \
    -b "$SORT_GENOME_BAM" \
    -o "$bigwig_output" \
    -p ${SLURM_CPUS_ON_NODE} \
    --normalizeUsing BPM \
    --binSize 10 \
    --minMappingQuality 255

if [ $? -ne 0 ]; then
    echo "ERROR: BigWig generation failed for sample ${geo_accession}"
    exit 1
fi

echo "BigWig generated: $bigwig_output"
echo "File size: $(du -h "$bigwig_output" | cut -f1)"
echo ""

echo "============================================"
echo "STAR_SE Processing Complete!"
echo "============================================"
echo "Sample:      $geo_accession"
echo "Unique name: $unique_name"
echo "Completed:   $(date "+%Y-%m-%d %H:%M:%S")"
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
