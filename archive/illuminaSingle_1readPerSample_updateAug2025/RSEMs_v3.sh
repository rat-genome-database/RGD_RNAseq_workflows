#!/bin/bash
#SBATCH --job-name=RSEM72s
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2gb
#SBATCH --time=02:00:00
#SBATCH --account=akwitek
#SBATCH --partition=normal
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wdemos@mcw.edu

# Written by Wendy Demos
# Last modified on 29 Oct 2024 by WMD
# Enable debugging by printing each command before execution
#set -x
# Software Modules
module load rsem/1.3.3 samtools/1.20

# Retrieve environment variables
Run=${Run}
geo_accession=${geo_accession}
PRJdir=${PRJdir}
BIOProjectID=${BIOProjectID}
echo "PRJdir is set to $PRJdir"
scratch_dir="/scratch/g/akwitek/wdemos/$BIOProjectID"
echo "Scratch working directory is set to $scratch_dir"
FinalOPdir=${Logdir}
echo "Logdir is set to: $Logdir"
mkdir -p "$scratch_dir"
tempOPdir="${scratch_dir}/${Run}"


# mRatBN7.2 Reference Files
REFdir="/home/wdemos/NCBI_1August2023/"
REF_GTF="/home/wdemos/NCBI_1August2023/mod_GCF_015227675.2_mRatBN7.2_genomic.gtf"
GENOME_FASTA="/home/wdemos/NCBI_1August2023/mod_GCF_015227675.2_mRatBN7.2_genomic.fna"
rsemref="${scratch_dir}/rsemref"
#depfiles="/home/wdemos/NCBI_1August2023/RSEMfiles"
echo "GTF reference file: $REF_GTF"
echo "Genome FASTA file: $GENOME_FASTA"

# Ensure project directory exists
if [ ! -d "$PRJdir" ]; then
    echo "Creating project directory: $PRJdir"
    mkdir -p "$PRJdir" || { echo "Error: Failed to create directory $PRJdir"; exit 1; }
else
    echo "Project directory already exists: $PRJdir"
fi

# Redirect output to the log file in the project directory
log_file="${PRJdir}/${geo_accession}/log_files/${geo_accession}_RSEM.out"
mkdir -p "$(dirname "$log_file")"
exec 1>"$log_file" 2>&1

# Capture start time
start_time=$(date +%s)

# Move into the working directory
echo "Moving into the working directory ${scratch_dir}"
cd "$scratch_dir" || { echo "Error: Failed to move into directory $scratch_dir"; exit 1; }

# Define sample ID
sample="$geo_accession"
echo "Processing sample for project: $Run"

# Define path to the _STARAligned.toTranscriptome.out.bam file
inbam="${scratch_dir}/${Run}/${geo_accession}_STARAligned.toTranscriptome.out.bam"
if [ ! -f "$inbam" ]; then
    echo "Error: Input BAM file not found for sample ${sample}. Skipping."
    exit 1
fi

# Calculate fragment length mean and standard deviation using samtools stats
read MEAN SD < <(samtools stats "$inbam" | awk '
    /average length/ { mean=$4 }
    /insert size standard deviation/ { sd=$6 }
    END { print mean, sd }')
echo "Extracted (average length)mean: $MEAN, insert size standard deviation: $SD"
# Ensure values are valid
if [[ -z "$MEAN" || -z "$SD" || "$MEAN" == "nan" || "$SD" == "nan" ]]; then
    echo "Error: Unable to estimate fragment length for sample ${geo_accession}. Check BAM file."
    exit 1
fi

echo "Estimating expression for sample: $geo_accession JOB START $(date "+%Y-%m-%d %H:%M:%S")"

# Estimate expression with RSEM for single-end reads
rsem-calculate-expression \
    --sort-bam-by-coordinate \
    --fragment-length-mean "$MEAN" \
    --fragment-length-sd "$SD" \
    --alignments -p 8 "$inbam" "$rsemref" "$tempOPdir/$geo_accession" 2>> "$log_file"

# Check for errors
if [ $? -ne 0 ]; then
    echo "Error: RSEM expression estimation failed for sample ${geo_accession}"
    exit 1
fi

#copy isoforms.results file to handle renaming to transcripts.results
cp "$tempOPdir/$geo_accession.isoforms.results" $tempOPdir/$geo_accession.transcripts.results

# Remove intermediate transcript bam files only upon successful processing
# rm "${scratch_dir}/${Run}/${geo_accession}.transcript.bam"
# rm "${scratch_dir}/${Run}/${geo_accession}_STARAligned.toTranscriptome.out.bam"
# cp "$tempOPdir/$geo_accession.isoforms.results" $FinalOPdir
# cp "$tempOPdir/$geo_accession.genes.results" $FinalOPdir

# Capture end time
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
echo "Total elapsed time: $(date -d@$elapsed_time -u +%H:%M:%S)"
module unload rsem/1.3.3 samtools/1.20
