#!/bin/sh
#SBATCH --job-name=RSEM8
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=10:00:00
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

# NOTE: Output and error logs are unified via submission script using --output and --error

# ---------------------------------------------
# Script: RSEM_noBW.bash
# Last modified: 17 Feb 2026 by WMD
# Purpose: Estimate expression using RSEM after STAR alignment
# BigWig generation has been moved to STAR_bigwig.sh
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
echo "Scratch working directory is set to $scratch_dir"
tempOPdir="${scratch_dir}/${geo_accession}"
echo "tempOPdir is set to $tempOPdir"
unique_name="${unique_name}"
echo "unique_name is set to ${unique_name}"

# Reference Files (GRCr8)
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
echo "Processing sample $sample"

# Define path to BAM & verify its existence
inbam="${scratch_dir}/${geo_accession}/${geo_accession}_STARAligned.toTranscriptome.out.bam"
if [ ! -f "$inbam" ]; then
    echo "Error: Input BAM file not found for sample ${sample}."
    echo "Expected location: $inbam"
    echo "Please verify STAR alignment completed successfully."
    ls -lh "${scratch_dir}/${geo_accession}/" || echo "Directory does not exist"
    exit 1
fi

echo "Transcriptome BAM found: $inbam"
echo "File size: $(du -h "$inbam" | cut -f1)"
echo ""

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

echo "RSEM completed successfully"
echo ""

# Validate RSEM output files
if [ ! -f "$tempOPdir/$geo_accession.genes.results" ]; then
    echo "ERROR: RSEM genes.results file missing for ${geo_accession}"
    echo "Expected: $tempOPdir/$geo_accession.genes.results"
    exit 1
fi

if [ ! -f "$tempOPdir/$geo_accession.isoforms.results" ]; then
    echo "ERROR: RSEM isoforms.results file missing for ${geo_accession}"
    echo "Expected: $tempOPdir/$geo_accession.isoforms.results"
    exit 1
fi

echo "RSEM output files validated:"
echo "  - genes.results: $(wc -l < "$tempOPdir/$geo_accession.genes.results") lines"
echo "  - isoforms.results: $(wc -l < "$tempOPdir/$geo_accession.isoforms.results") lines"
echo ""

echo "Copying isoforms.results to transcripts.results"
cp "$tempOPdir/$geo_accession.isoforms.results" "$tempOPdir/$geo_accession.transcripts.results"
echo ""

echo "============================================"
echo "RSEM Processing Complete!"
echo "============================================"
echo "Sample: $geo_accession"
echo "Unique name: $unique_name"
echo "Processing completed: $(date "+%Y-%m-%d %H:%M:%S")"
echo ""
echo "Output files in sample scratch directory ($tempOPdir):"
echo "  - ${geo_accession}.genes.results"
echo "  - ${geo_accession}.transcripts.results"
echo "  - ${geo_accession}.isoforms.results"
echo ""

# Report elapsed time
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
echo "Total elapsed time: $(date -d@$elapsed_time -u +%H:%M:%S)"
echo "============================================"