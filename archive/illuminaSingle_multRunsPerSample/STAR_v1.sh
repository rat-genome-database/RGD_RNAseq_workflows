#!/bin/sh
#SBATCH --job-name=mStar
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=3gb
#SBATCH --time=01:30:00
#SBATCH --account=akwitek
#SBATCH --partition=normal
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wdemos@mcw.edu

# Load modules
module load star/2.7.10b samtools/1.20 deeptools/3.5.1

# Read sample ID and corresponding reads from command-line arguments
SAMPLE_ID=$1
READ1_FILES=$2  # Comma-separated list of forward reads
BIOProjectID=$3  # Project ID
unique_name=$4

scratch_dir="/scratch/g/akwitek/wdemos/${BIOProjectID}"
wrkDir="${scratch_dir}/${SAMPLE_ID}"
OUTPUT_PREFIX="${wrkDir}/${SAMPLE_ID}_STAR"
GENOME_BAM="${OUTPUT_PREFIX}Aligned.out.bam"
SORT_GENOME_BAM="${wrkDir}/${SAMPLE_ID}_GENOME_SORT.bam"
baseDir="${baseDir}"
Logdir="${Logdir}" #/log_files
PRJdir="${PRJdir}" #/data/expression/GEO/${BIOProjectID}/reads_fastq
FinalOPdir="${PRJdir}/${SAMPLE_ID}"
GENOME_DIR="${scratch_dir}/RefIndex"
REF_GTF="/home/wdemos/NCBI_1August2023/mod_GCF_015227675.2_mRatBN7.2_genomic.gtf"
GENOME_FASTA="/home/wdemos/NCBI_1August2023/mod_GCF_015227675.2_mRatBN7.2_genomic.fna"
INDEX_DIR="$scratch_dir/RefIndex"
BIGWIG_FILE="${FinalOPdir}/RNAseq_${unique_name}.bigwig"
# Ensure directories exist
mkdir -p "$scratch_dir"
mkdir -p "$wrkDir"
mkdir -p "$FinalOPdir"

# Remove trailing commas
READ1_FILES="${READ1_FILES%,}"

# Check if we have reads to process
if [[ -z "$READ1_FILES" ]]; then
    echo "Error: No valid fastq files found for sample $SAMPLE_ID."
    exit 1
fi

echo "READ1_FILES: \"${READ1_FILES}\""

# Run STAR alignment once per sample
echo "Processing sample: $sample JOB START $(date "+%Y-%m-%d %H:%M:%S")"

#outFilterMatchNminOverLread 0.66 same as outFilterMatchNmin, but normalized to read length if 100 paired reads n matches = 0.66*200=132 bases must match
#outFilterMultimapNmax default = 20 max number of multiple alignments allowed for a read: if exceeded, the read is considered unmapped,
#a value of 6 Allows a read to map to at most 6 locations, reducing ambiguity albeit strict for rRNA or paralogous genes
STAR --runThreadN 10 \
     --readFilesCommand zcat \
     --genomeDir "$INDEX_DIR" \
     --outFileNamePrefix "$OUTPUT_PREFIX" \
     --quantMode TranscriptomeSAM \
     --outSAMtype BAM Unsorted \
     --outFilterMatchNminOverLread 0.66 \
     --outFilterMultimapNmax 6 \
     --readFilesIn "$READ1_FILES"


# Check if STAR completed successfully
if [ $? -ne 0 ]; then
    echo "Error: STAR alignment failed for sample ${SAMPLE_ID}"
    exit 1
fi

echo "STAR alignment for ${SAMPLE_ID} completed."

# Sort and index BAM file for BigWig generation
samtools sort -@ 10 "$GENOME_BAM" -o "$SORT_GENOME_BAM"
samtools index -b "$SORT_GENOME_BAM"

# Check if sorting was successful
if [ ! -f "$SORT_GENOME_BAM" ]; then
    echo "Error: Genome sorted BAM file was not created for sample ${SAMPLE_ID}."
    exit 1
fi

# Generate BigWig files
bamCoverage -b "$SORT_GENOME_BAM" -o "$BIGWIG_FILE" -p ${SLURM_CPUS_ON_NODE} --normalizeUsing BPM --binSize 50

# Check for successful BigWig generation
if [ $? -ne 0 ]; then
    echo "Error: bamCoverage failed for sample ${SAMPLE_ID}"
    exit 1
fi
#Added 15 April 2025 keep the STAR Log file with the data
cp  "$wrkDir/${SAMPLE_ID}_STARLog.final.out" "$FinalOPdir/log_files"

# Final confirmation
echo "Final output directory: $FinalOPdir"
echo "$BIGWIG_FILE has been generated."

