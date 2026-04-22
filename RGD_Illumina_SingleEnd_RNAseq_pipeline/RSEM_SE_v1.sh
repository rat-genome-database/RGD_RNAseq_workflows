#!/usr/bin/env bash
#SBATCH --job-name=RSEM_SE
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=10:00:00
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --partition=normal
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=your@email.edu
# NOTE: Output and error logs are set via submission script using --output and --error

###############################################################################
# USER CONFIGURATION
# Update the variables below for your system before running.
###############################################################################
# Scratch filesystem root; scratch_dir is built as SCRATCH_BASE/BIOProjectID
SCRATCH_BASE="/path/to/scratch"
# GRCr8 GTF annotation file
REF_GTF="/path/to/genome/GCF_036323735.1_GRCr8/updated_GCF_036323735.1_CM070413.1_GRCr8_genomic.gtf"
# GRCr8 genome FASTA file
GENOME_FASTA="/path/to/genome/GCF_036323735.1_GRCr8/updated_GCF_036323735.1_CM070413.1_GRCr8_genomic.fna"
#mRatBN7.2 for validation testing
#REF_GTF="/path/to/genome/NCBI_1August2023/mod_GCF_015227675.2_mRatBN7.2_genomic.gtf"
#GENOME_FASTA="/path/to/genome/NCBI_1August2023/mod_GCF_015227675.2_mRatBN7.2_genomic.fna"

###############################################################################

# ---------------------------------------------
# Script: RSEM_SE_v1.sh
# Purpose: Estimate expression using RSEM after single-end STAR alignment.
#
# Key difference from paired-end:
#   - No --paired-end flag
#   - Requires --fragment-length-mean and --fragment-length-sd,
#     estimated from samtools stats on the transcriptome BAM
# ---------------------------------------------

module load rsem/1.3.3 samtools/1.20

# Retrieve environment variables
Run=${Run}
geo_accession=${geo_accession}
PRJdir=${PRJdir}
BIOProjectID=${BIOProjectID}
Logdir=${Logdir}
unique_name="${unique_name}"

echo "PRJdir is set to $PRJdir"
scratch_dir="${SCRATCH_BASE}/${BIOProjectID}"
echo "Scratch working directory is set to $scratch_dir"
tempOPdir="${scratch_dir}/${geo_accession}"
echo "tempOPdir is set to $tempOPdir"
echo "unique_name is set to ${unique_name}"

rsemref="${scratch_dir}/rsemref"

echo "GTF reference file: $REF_GTF"
echo "Genome FASTA file:  $GENOME_FASTA"

# Ensure project directory exists
if [ ! -d "$PRJdir" ]; then
    echo "Creating project directory: $PRJdir"
    mkdir -p "$PRJdir" || { echo "Error: Failed to create directory $PRJdir"; exit 1; }
else
    echo "Project directory already exists: $PRJdir"
fi

# Redirect output to log file
log_file="${PRJdir}/${geo_accession}/log_files/RSEM/${geo_accession}_RSEM.out"
mkdir -p "$(dirname "$log_file")"
exec 1>"$log_file" 2>&1

# Capture start time
start_time=$(date +%s)

echo "Moving into the working directory ${scratch_dir}"
cd "$scratch_dir" || { echo "Error: Failed to move into directory $scratch_dir"; exit 1; }

sample="$geo_accession"
echo "Processing sample: $sample"

# Define path to transcriptome BAM
inbam="${scratch_dir}/${geo_accession}/${geo_accession}_STARAligned.toTranscriptome.out.bam"
if [ ! -f "$inbam" ]; then
    echo "Error: Input BAM file not found for sample ${sample}."
    echo "Expected: $inbam"
    echo "Please verify STAR alignment completed successfully."
    ls -lh "${scratch_dir}/${geo_accession}/" || echo "Directory does not exist"
    exit 1
fi

echo "Transcriptome BAM found: $inbam"
echo "File size: $(du -h "$inbam" | cut -f1)"
echo ""

# ---------------------------------------------------------------------------
# Estimate fragment length mean and SD from samtools stats
# Required for single-end RSEM (no insert size information from read pairs)
# ---------------------------------------------------------------------------
echo "Estimating fragment length parameters from BAM stats..."
read MEAN SD < <(samtools stats "$inbam" | awk '
    /average length/                 { mean=$4 }
    /insert size standard deviation/ { sd=$6   }
    END { print mean, sd }')

echo "Fragment length mean (average read length): $MEAN"
echo "Fragment length SD:                         $SD"

if [[ -z "$MEAN" || -z "$SD" || "$MEAN" == "nan" || "$SD" == "nan" ]]; then
    echo "Error: Unable to estimate fragment length for sample ${geo_accession}. Check BAM file."
    exit 1
fi
echo ""

echo "Estimating expression for sample: $geo_accession JOB START $(date "+%Y-%m-%d %H:%M:%S")"

# ---------------------------------------------------------------------------
# Run RSEM — single-end mode
# No --paired-end flag; fragment length params required
# ---------------------------------------------------------------------------
rsem-calculate-expression \
    --sort-bam-by-coordinate \
    --fragment-length-mean "$MEAN" \
    --fragment-length-sd "$SD" \
    --alignments -p 8 "$inbam" "$rsemref" "$tempOPdir/$sample" 2>> "$log_file"

if [ $? -ne 0 ]; then
    echo "Error: RSEM expression estimation failed for sample ${sample}"
    exit 1
fi
echo "RSEM completed successfully"
echo ""

# Validate output files
if [ ! -f "$tempOPdir/${sample}.genes.results" ]; then
    echo "ERROR: RSEM genes.results file missing for ${sample}"
    exit 1
fi
if [ ! -f "$tempOPdir/${sample}.isoforms.results" ]; then
    echo "ERROR: RSEM isoforms.results file missing for ${sample}"
    exit 1
fi

echo "RSEM output files validated:"
echo "  - genes.results:    $(wc -l < "$tempOPdir/${sample}.genes.results") lines"
echo "  - isoforms.results: $(wc -l < "$tempOPdir/${sample}.isoforms.results") lines"
echo ""

echo "Copying isoforms.results to transcripts.results"
cp "$tempOPdir/${sample}.isoforms.results" "$tempOPdir/${sample}.transcripts.results"
echo ""

echo "============================================"
echo "RSEM_SE Processing Complete!"
echo "============================================"
echo "Sample:      $geo_accession"
echo "Unique name: $unique_name"
echo "Completed:   $(date "+%Y-%m-%d %H:%M:%S")"
echo ""
echo "Output files in sample scratch directory ($tempOPdir):"
echo "  - ${sample}.genes.results"
echo "  - ${sample}.transcripts.results"
echo "  - ${sample}.isoforms.results"
echo ""

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
echo "Total elapsed time: $(date -d@$elapsed_time -u +%H:%M:%S)"
echo "============================================"
