#!/bin/bash
#SBATCH --job-name=STARref
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=00:30:00
#SBATCH --account=akwitek
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --partition=normal
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wdemos@mcw.edu

# Script details
## Written by Wendy Demos
## Written 16 August 2024
## Last Modified 28 Oct 2024 by WMD to run on RCC cluster
# 18June2025 removed form header to reduce redundante log files #SBATCH --output=%x-%j.out

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
scratch_dir="/scratch/g/akwitek/wdemos/$BIOProjectID"
FinalOPdir=${Logdir}
echo "Logdir is set to: $Logdir"
mkdir -p "$scratch_dir"
FASTQ_DIR="$scratch_dir/$Run"
REF_GTF="/home/wdemos/NCBI_1August2023/mod_GCF_015227675.2_mRatBN7.2_genomic.gtf"
GENOME_FASTA="/home/wdemos/NCBI_1August2023/mod_GCF_015227675.2_mRatBN7.2_genomic.fna"
mkdir -p /scratch/g/akwitek/wdemos/$BIOProjectID/RefIndex/
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
         --runThreadN 8
fi
