#!/usr/bin/env bash
#SBATCH --job-name=STAR_SE
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=3gb
#SBATCH --time=08:00:00
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
# v2: Picard CollectRnaSeqMetrics reference files
RRNA_INTERVALS="$myDir/GCF_036323735.1_GRCr8/GRCr8_rRNA.intervals"
REF_FLAT="$myDir/GCF_036323735.1_GRCr8/updated_GCF_036323735.1_CM070413.1_GRCr8_genomic.gtf.refflat"

###############################################################################

# Library prep call script (shared with libprep_recall.sh). sbatch runs a spooled copy of this script,
# so default to the directory the job was submitted from; override with --export=LIBPREP_CALL=...
LIBPREP_CALL=${LIBPREP_CALL:-${SLURM_SUBMIT_DIR:-$(dirname "$0")}/libprep_call.sh}

# Script details
## v2 (6 Oct 2026): STAR GeneCounts, strand detection, Picard CollectRnaSeqMetrics and the library prep
##   call (libprep_call.sh, rule LP2) added. Wall time 8 h.
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

module load star/2.7.10b samtools/1.20 deeptools/3.5.1 picard/2.25.0

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
GENE_COUNTS_FILE="${FinalOPdir}/log_files/STAR/${geo_accession}_STARReadsPerGene.out.tab"
STRAND_DETECTION_LOG_FILE="${FinalOPdir}/log_files/STAR/${geo_accession}_strand_detection.log"
RNA_METRICS_FILE="${FinalOPdir}/log_files/STAR/${geo_accession}_rna_metrics.txt"
LIBRARY_PREP_LOG_FILE="${FinalOPdir}/log_files/STAR/${geo_accession}_library_prep.log"

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
    --quantMode TranscriptomeSAM GeneCounts \
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
if  [ -f "${OUTPUT_PREFIX}ReadsPerGene.out.tab" ]; then
    cp "${OUTPUT_PREFIX}ReadsPerGene.out.tab" "$FinalOPdir/log_files/STAR"
    echo "Copied: ${OUTPUT_PREFIX}ReadsPerGene.out.tab"
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
