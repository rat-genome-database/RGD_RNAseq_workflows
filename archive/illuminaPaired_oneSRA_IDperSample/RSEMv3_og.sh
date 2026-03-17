#!/bin/sh
#SBATCH --job-name=RSEM72
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=03:00:00
#SBATCH --account=akwitek
#SBATCH --partition=normal
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wdemos@mcw.edu
# NOTE: Output and error logs are unified via submission script using --output and --error

# ---------------------------------------------
# Script: RSEM.sh
# Last modified: June 18, 2025 by WMD
# Purpose: Estimate expression using RSEM after STAR alignment
# ---------------------------------------------

# Load required software modules
module load rsem/1.3.3 samtools/1.20

# Retrieve environment variables
Run=${Run}
geo_accession=${geo_accession}
PRJdir=${PRJdir}
BIOProjectID=${BIOProjectID}
Logdir=${Logdir}

echo "PRJdir is set to $PRJdir"
scratch_dir="/scratch/g/akwitek/wdemos/$BIOProjectID"
echo "Scratch working directory is set to $scratch_dir"
tempOPdir="${scratch_dir}/${Run}"

# Reference Files (mRatBN7.2)
REFdir="/home/wdemos/NCBI_1August2023/"
REF_GTF="/home/wdemos/NCBI_1August2023/mod_GCF_015227675.2_mRatBN7.2_genomic.gtf"
GENOME_FASTA="/home/wdemos/NCBI_1August2023/mod_GCF_015227675.2_mRatBN7.2_genomic.fna"
rsemref="${scratch_dir}/rsemref"

echo "GTF reference file: $REF_GTF"
echo "Genome FASTA file: $GENOME_FASTA"

# Ensure project directory exists
if [ ! -d "$PRJdir" ]; then
    echo "Creating project directory: $PRJdir"
    mkdir -p "$PRJdir" || { echo "Error: Failed to create directory $PRJdir"; exit 1; }
else
    echo "Project directory already exists: $PRJdir"
fi

# Capture start time
start_time=$(date +%s)

# Move into the working directory
echo "Moving into the working directory ${scratch_dir}"
cd "$scratch_dir" || { echo "Error: Failed to move into directory $scratch_dir"; exit 1; }

# Define sample ID
sample="$geo_accession"
echo "Processing sample $sample (Run: $Run)"

# Define path to BAM & verify it's existence
inbam="${scratch_dir}/${Run}/${geo_accession}_STARAligned.toTranscriptome.out.bam"
if [ ! -f "$inbam" ]; then
    echo "Error: Input BAM file not found for sample ${sample}. Skipping."
    exit 1
fi

# Run RSEM
echo "Estimating expression for sample: $geo_accession JOB START $(date "+%Y-%m-%d %H:%M:%S")"
echo "Running RSEM on $inbam with reference $rsemref"

rsem-calculate-expression \
    --sort-bam-by-coordinate \
    --paired-end \
    --alignments -p 8 "$inbam" "$rsemref" "$tempOPdir/$geo_accession"

if [ $? -ne 0 ]; then
    echo "Error: RSEM expression estimation failed for sample ${geo_accession}"
    exit 1
fi

# Validate output files
if [ ! -f "$tempOPdir/$geo_accession.genes.results" ] || [ ! -f "$tempOPdir/$geo_accession.isoforms.results" ]; then
    echo "Error: One or more RSEM output files missing for ${geo_accession}"
    exit 1
fi

# Cleanup input BAM
echo "RSEM completed successfully for ${geo_accession}. Removing input BAM ${inbam}"
#rm -f "$inbam"

# Wrap up
echo "Processing completed for sample $geo_accession on $(date "+%Y-%m-%d %H:%M:%S")"
echo "Output directory: ${Logdir}/reads_fastq/${geo_accession}/"

# Report elapsed time
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
echo "Total elapsed time: $(date -d@$elapsed_time -u +%H:%M:%S)"
