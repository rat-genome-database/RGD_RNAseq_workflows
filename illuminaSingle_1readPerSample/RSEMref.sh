#!/bin/bash
#SBATCH --job-name=RSEMref
#SBATCH --ntasks=1
#SBATCH --mem=2gb
#SBATCH --time=00:10:00
#SBATCH --output=RSEMref-%j.out
#SBATCH --error=RSEMref-%j.err

# Software Modules
module load rsem/1.3.3 samtools/1.20

# Arguments
BIOProjectID=${BIOProjectID}
scratch_dir="/scratch/g/akwitek/wdemos/$BIOProjectID"
REFdir="/home/wdemos/NCBI_1August2023"
#REF_GTF="/home/wdemos/NCBI_1August2023/rn7chr.gtf"
REF_GTF="/home/wdemos/NCBI_1August2023/mod_GCF_015227675.2_mRatBN7.2_genomic.gtf"
GENOME_FASTA="/home/wdemos/NCBI_1August2023/mod_GCF_015227675.2_mRatBN7.2_genomic.fna"
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
