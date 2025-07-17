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

# Set input parameters
#GENOME_DIR="/scratch/g/akwitek/wdemos/GSE219045/RefIndex"  # Path to STAR genome index
#REF_GTF="/home/wdemos/NCBI_1August2023/mod_GCF_015227675.2_mRatBN7.2_genomic.gtf"
#GENOME_FASTA="/home/wdemos/NCBI_1August2023/mod_GCF_015227675.2_mRatBN7.2_genomic.fna"

# Read sample ID and corresponding reads from command-line arguments
SAMPLE_ID=$1
READ1_FILES=$2  # Comma-separated list of forward reads
READ2_FILES=$3  # Comma-separated list of reverse reads
BIOProjectID=$4  # Project ID
#Run=$5  # Run ID
unique_name=$5

scratch_dir="/scratch/g/akwitek/wdemos/${BIOProjectID}"
wrkDir="/scratch/g/akwitek/wdemos/${BIOProjectID}/${SAMPLE_ID}"
OUTPUT_PREFIX="${wrkDir}/${SAMPLE_ID}_STAR"
GENOME_BAM="${OUTPUT_PREFIX}Aligned.out.bam"
SORT_GENOME_BAM="${wrkDir}/${SAMPLE_ID}_GENOME_SORT.bam"
PRJdir="/home/wdemos/data/expression/GEO/${BIOProjectID}/reads_fastq"
FinalOPdir="${PRJdir}/${SAMPLE_ID}"
#unique_name="${unique_name}"
GENOME_DIR="${scratch_dir}/RefIndex"
REF_GTF="/home/wdemos/NCBI_1August2023/mod_GCF_015227675.2_mRatBN7.2_genomic.gtf"
GENOME_FASTA="/home/wdemos/NCBI_1August2023/mod_GCF_015227675.2_mRatBN7.2_genomic.fna"

# Ensure directories exist
mkdir -p "$scratch_dir"
mkdir -p "$wrkDir"
mkdir -p "$FinalOPdir"

if [[ -z "$READ1_FILES" || -z "$READ2_FILES" ]]; then
echo "Error: No valid fastq files found for sample $SAMPLE_ID."
exit 1
fi

#echo "Troubleshoot: Checking if READ1_FILES and READ2_FILES are passed correctly" >> "$log_file"
#echo "READ1_FILES: \"${READ1_FILES}\"" >> "$log_file"
#echo "READ2_FILES: \"${READ2_FILES}\"" >> "$log_file"

# Remove trailing commas
READ1_FILES="${READ1_FILES%,}"
READ2_FILES="${READ2_FILES%,}"

# Check if we have reads to process
if [[ -z "$READ1_FILES" || -z "$READ2_FILES" ]]; then
echo "Error: No valid fastq files found for sample $SAMPLE_ID."
exit 1
fi
#echo "Running STAR with inputs: ${READ1_FILES} ${READ2_FILES}" >> "$log_file"

# Run STAR alignment once per sample
STAR --runThreadN 10 \
--readFilesCommand zcat \
--genomeDir "$GENOME_DIR" \
--outFileNamePrefix "$OUTPUT_PREFIX" \
--quantMode TranscriptomeSAM \
--outSAMtype BAM Unsorted \
--outSAMunmapped Within KeepPairs \
--outFilterMatchNminOverLread 0.66 \
--outFilterMultimapNmax 6 \
--readFilesIn "${READ1_FILES}" "${READ2_FILES}"

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

# Check if BigWig file already exists
BIGWIG_FILE="$FinalOPdir/RNAseq_${unique_name}.bigwig"
if [ -f "$BIGWIG_FILE" ]; then
echo "Skipping sample ${SAMPLE_ID} as the bigWig file has already been generated."
exit 0
fi

# Generate BigWig files
bamCoverage -b "$SORT_GENOME_BAM" -o "$BIGWIG_FILE" -p ${SLURM_CPUS_ON_NODE} --normalizeUsing BPM --binSize 50

# Check for successful BigWig generation
if [ $? -ne 0 ]; then
echo "Error: bamCoverage failed for sample ${SAMPLE_ID}"
exit 1
fi

# Final confirmation
echo "Final output directory: $FinalOPdir"
echo "$BIGWIG_FILE has been generated."
