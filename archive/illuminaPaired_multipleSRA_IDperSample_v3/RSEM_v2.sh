#!/bin/bash
#SBATCH --job-name=RSEM72
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2gb
#SBATCH --time=03:00:00
#SBATCH --account=akwitek
#SBATCH --partition=normal
#SBATCH --output=%x-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wdemos@mcw.edu
#SBATCH --error=%x-%j.err

# Written by Wendy Demos
# Last modified on 28 Feb 2025 by WMD

# Enable debugging by printing each command before execution
#set -x

# Software Modules
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
echo "Logdir is set to: $Logdir"

mkdir -p "$scratch_dir"
tempOPdir="${scratch_dir}/${geo_accession}"

# mRatBN7.2 Reference Files
REFdir="/home/wdemos/NCBI_1August2023/"
REF_GTF="/home/wdemos/NCBI_1August2023/mod_GCF_015227675.2_mRatBN7.2_genomic.gtf"
GENOME_FASTA="/home/wdemos/NCBI_1August2023/mod_GCF_015227675.2_mRatBN7.2_genomic.fna"
rsemref="${scratch_dir}/rsemref"

echo "GTF reference file: $REF_GTF"
echo "Genome FASTA file: $GENOME_FASTA"

# Ensure project directory and log directory exist
if [ ! -d "$PRJdir" ]; then
    echo "Creating project directory: $PRJdir"
    mkdir -p "$PRJdir" || { echo "Error: Failed to create directory $PRJdir"; exit 1; }
else
    echo "Project directory already exists: $PRJdir"
fi

mkdir -p "${PRJdir}/${geo_accession}/log_files/RSEM"

# Redirect output to the log file in the project directory
log_file="${PRJdir}/${geo_accession}/log_files/RSEM/${geo_accession}_RSEM.out"
exec 1>"$log_file" 2>&1

# Capture start time
start_time=$(date +%s)

# Move into the working directory
echo "Moving into the working directory ${scratch_dir}"
cd "$scratch_dir" || { echo "Error: Failed to move into directory $scratch_dir"; exit 1; }

# Define sample ID
sample="$geo_accession"
echo "Processing sample for project: $geo_accession"

# Define path to the _STARAligned.toTranscriptome.out.bam file
inbam="${scratch_dir}/${geo_accession}/${geo_accession}_STARAligned.toTranscriptome.out.bam"
if [ ! -f "$inbam" ]; then
    echo "Error: Input BAM file not found for sample ${sample}. Skipping."
    exit 1
fi

echo "Estimating expression for sample: $geo_accession JOB START $(date "+%Y-%m-%d %H:%M:%S")"

# Estimate expression with RSEM
rsem-calculate-expression \
    --sort-bam-by-coordinate \
    --paired-end \
    --alignments -p 8 "$inbam" "$rsemref" "$tempOPdir/$geo_accession"

if [ $? -ne 0 ]; then
    echo "Error: RSEM expression estimation failed for sample ${geo_accession}"
    exit 1
fi

echo "Processing completed for sample $geo_accession on $(date "+%Y-%m-%d %H:%M:%S")"

# Copy isoforms.results file to transcripts.results in scratch directory
echo "Copy isoforms.results file to transcripts.results in scratch directory"
cp "$tempOPdir/$geo_accession.isoforms.results" "${scratch_dir}/${geo_accession}/${geo_accession}.transcripts.results"

# Capture end time
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
echo "Total elapsed time: $(date -d@$elapsed_time -u +%H:%M:%S)"
