#!/bin/bash
#SBATCH --job-name=RSEMref
#SBATCH --ntasks=1
#SBATCH --mem=2gb
#SBATCH --account=your-slurm-account
#SBATCH --time=00:10:00
#SBATCH --output=RSEMref-%j.log
#SBATCH --error=RSEMref-%j.log

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


# Software Modules
module load rsem/1.3.3 samtools/1.20

# Variables
BIOProjectID=${BIOProjectID}
# Reference Files (GRCr8)
rsemref="${scratch_dir}/rsemref"

echo "RSEM reference will be generated at: $rsemref"

# Check if RSEM reference already exists
if [ -f "$rsemref/rsemref.transcripts.fa" ]; then
    echo "RSEM reference already exists. Skipping generation."
    exit 0
fi

# Generate RSEM reference
echo "RSEM reference not found. Generating now..."
rsem-prepare-reference --gtf "$REF_GTF" "$GENOME_FASTA" "$rsemref"

if [ $? -eq 0 ]; then
    echo "RSEM reference generated successfully."
else
    echo "Error generating RSEM reference."
    exit 1
fi
