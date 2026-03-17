#!/bin/bash
#SBATCH --job-name=STAR72
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=3gb
#SBATCH --time=05:00:00
#SBATCH --account=akwitek
#SBATCH --partition=normal
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wdemos@mcw.edu
# NOTE: Output and error logs are unified via submission script using --output and --error

# Script details
## Written by Wendy Demos
## Written 16 August 2024
## Last Modified 6 Nov 2024 by WMD to run on RCC cluster
#18 June 2025 removed form header to reduce log files: #SBATCH --output=%x-%j.out

# Enable debugging by printing each command before execution
#set -x
module load star/2.7.10b  samtools/1.20

# Retrieve environment variables
Run=${Run}
geo_accession=${geo_accession}
PRJdir=${PRJdir} ###"$myDir/data/expression/GEO/$BIOProjectID/reads_fastq"
LENGTH=${length}
BIOProjectID=${BIOProjectID}
unique_name="${unique_name}"
scratch_dir="/scratch/g/akwitek/wdemos/$BIOProjectID"
FinalOPdir=${PRJdir}/${geo_accession}
echo "Logdir is set to: $Logdir" #"$myDir/data/expression/GEO/$BIOProjectID"
echo " FinalOPdir= $PRJdir/$geo_accession"
mkdir -p "$scratch_dir"
FASTQ_DIR="$scratch_dir/$Run"
REF_GTF="/home/wdemos/GCF_036323735.1_GRCr8/updated_GCF_036323735.1_CM070413.1_GRCr8_genomic.gtf"
GENOME_FASTA="/home/wdemos/GCF_036323735.1_GRCr8/updated_GCF_036323735.1_CM070413.1_GRCr8_genomic.fna"
INDEX_DIR="$scratch_dir/RefIndex"

# Verify project directory exists
if [ ! -d "$scratch_dir/$Run" ]; then
    echo "Error: Directory $scratch_dir/$Run does not exist. Please check the project name."
    exit 1
fi

echo "Project directory: $scratch_dir/$Run"
echo "GTF reference file: $REF_GTF"
echo "Genome FASTA file: $GENOME_FASTA"
echo "INDEX_DIR: $INDEX_DIR"


# Define sample ID
sample="$geo_accession"
echo "Processing sample for project: $Run"

# Define paths to the FASTQ files
R1=$(find "$scratch_dir/$Run/" -type f -name "${sample}_1.fastq.gz" | head -n 1)
R2=$(find "$scratch_dir/$Run/" -type f -name "${sample}_2.fastq.gz" | head -n 1)

if [ -z "$R1" ] || [ -z "$R2" ]; then
    echo "Error: fastq.gz files not found for sample ${sample}. Skipping."
    exit 1
fi

OUTPUT_PREFIX="$scratch_dir/$Run/${sample}_STAR"
LOG_OUTPUT="${OUTPUT_PREFIX}Log.final.out"
SORT_GENOME_BAM="$scratch_dir/$Run/${sample}_GENOME_SORT.bam"
GENOME_BAM="${OUTPUT_PREFIX}Aligned.out.bam"

# Check if the output BAM files and final log already exists
if [ -f "$LOG_OUTPUT" ] && [ -f "$SORT_GENOME_BAM" ]; then
    echo "Skipping sample ${sample} as it has already been processed."
    exit 0
fi

echo "Processing sample: $sample JOB START $(date "+%Y-%m-%d %H:%M:%S")"

#updated 3 Jan 2025 by WD
#outFilterMatchNminOverLread 0.66 same as outFilterMatchNmin, but normalized to read length if 100 paired reads n matches = 0.66*200=132 bases must match
#outFilterMultimapNmax default = 20 max number of multiple alignments allowed for a read: if exceeded, the read is considered unmapped,
#a value of 6 Allows a read to map to at most 6 locations, reducing ambiguity albeit strict for rRNA or paralogous genes
STAR --runThreadN 10 \
     --readFilesCommand zcat \
     --genomeDir "$INDEX_DIR" \
     --outFileNamePrefix "$OUTPUT_PREFIX" \
     --quantMode TranscriptomeSAM \
     --outSAMtype BAM Unsorted \
     --outSAMunmapped Within KeepPairs \
     --outFilterMatchNminOverLread 0.66 \
     --outFilterMultimapNmax 6 \
     --readFilesIn "$R1" "$R2"

if [ $? -ne 0 ]; then
    echo "Error: STAR alignment failed for sample ${sample}"
    exit 1
fi

# Sort and index Genome BAM file for computation of sample sex
samtools sort -@ 10 "$GENOME_BAM" -o "$scratch_dir/$Run/${sample}_GENOME_SORT.bam"
samtools index -b "$scratch_dir/$Run/${sample}_GENOME_SORT.bam"

# Check for successful BAM sorting
if [ ! -f "$scratch_dir/$Run/${sample}_GENOME_SORT.bam" ]; then
    echo "Error: Genome sorted BAM file was not created for sample ${sample}."
    exit 1
fi

#Added 15 April 2025 keep the STAR Log file with the data
cp  "$scratch_dir/$Run/${sample}_STARLog.final.out" "$FinalOPdir/log_files"
echo "unique_name: $unique_name"
echo "Final output directory: $FinalOPdir"


