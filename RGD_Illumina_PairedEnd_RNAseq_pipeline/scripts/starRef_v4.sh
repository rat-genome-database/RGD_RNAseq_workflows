#!/bin/bash
#SBATCH --job-name=STARref
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=7gb
#SBATCH --time=02:00:00
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


# Script details
## Written by Wendy Demos
## Written 16 August 2024
## Last Modified 5 Feb 2026 by WMD to run on RCC cluster with GRCr8 reference

# Enable debugging by printing each command before execution
set -x
module load star/2.7.10b  samtools/1.20 deeptools/3.5.1

# Retrieve environment variables
Run=${Run}
geo_accession=${geo_accession}
PRJdir=${PRJdir}
LENGTH=${length}
BIOProjectID=${BIOProjectID}
unique_name="${unique_name}"
FinalOPdir=${Logdir}
echo "Logdir is set to: $Logdir"
mkdir -p "$scratch_dir"
FASTQ_DIR="$scratch_dir/$Run"

mkdir -p "${SCRATCH_BASE}/${BIOProjectID}/RefIndex/"
INDEX_DIR="$scratch_dir/RefIndex"

# Verify project directory exists
if [ ! -d "$scratch_dir" ]; then
    echo "Error: Directory $scratch_dir does not exist. Please check the project name.(This is the first catch)"
    exit 1
fi

# Check if the reference index directory exists, if not, create it
if [ ! -d "$INDEX_DIR" ]; then
    echo "Directory $INDEX_DIR does not exist. Creating it."
    mkdir -p "$INDEX_DIR"
fi

echo "Project directory: $scratch_dir"
echo "GTF reference file: $REF_GTF"
echo "Genome FASTA file: $GENOME_FASTA"
echo "INDEX_DIR: $INDEX_DIR"

# Check if the index directory exists and contains key files
if [ -d "$INDEX_DIR" ] && [ -f "$INDEX_DIR/SA" ] && [ -f "$INDEX_DIR/SAindex" ] && [ -f "$INDEX_DIR/Genome" ] && [ -f "$INDEX_DIR/chrNameLength.txt" ] && [ -f "$INDEX_DIR/chrName.txt" ]; then
    echo "Genome index already exists. Skipping genomeGenerate."
#    exit 0 # explicitly return success to avoid error with dependency for next step
else
    echo "Genome index not found. Running genomeGenerate."
# Calculate sjdbOverhang as LENGTH - 1
SJDB_OVERHANG=$((LENGTH - 1))
    STAR --runMode genomeGenerate \
         --genomeDir "$INDEX_DIR" \
         --genomeFastaFiles "$GENOME_FASTA" \
         --sjdbGTFfile "$REF_GTF" \
         --sjdbOverhang "$SJDB_OVERHANG" \
         --outFileNamePrefix "$INDEX_DIR/" \
         --runThreadN 10 \
         --limitGenomeGenerateRAM 50000000000
fi
