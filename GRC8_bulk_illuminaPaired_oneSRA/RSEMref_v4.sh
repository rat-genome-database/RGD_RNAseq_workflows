#!/bin/bash
#SBATCH --job-name=RSEMref
#SBATCH --ntasks=1
#SBATCH --mem=2gb
#SBATCH --account=akwitek
#SBATCH --time=00:10:00
#SBATCH --output=RSEMref-%j.log
#SBATCH --error=RSEMref-%j.log

# Software Modules
module load rsem/1.3.3 samtools/1.20

# Variables
BIOProjectID=${BIOProjectID}
scratch_dir="/scratch/g/akwitek/wdemos/$BIOProjectID"
# Reference Files (GRCr8)
REFdir="/home/wdemos/GCF_036323735.1_GRCr8"
REF_GTF="/home/wdemos/GCF_036323735.1_GRCr8/updated_GCF_036323735.1_CM070413.1_GRCr8_genomic.gtf"
GENOME_FASTA="/home/wdemos/GCF_036323735.1_GRCr8/updated_GCF_036323735.1_CM070413.1_GRCr8_genomic.fna"
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
