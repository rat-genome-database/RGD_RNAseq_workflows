#!/bin/sh
#SBATCH --job-name=RSEM8
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=10:00:00
#SBATCH --account=YOUR_SLURM_ACCOUNT       # <-- replace with your SLURM account
#SBATCH --partition=normal
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=YOUR_EMAIL@example.com  # <-- replace with your email
# NOTE: Output and error logs are unified via submission script using --output and --error

# ---------------------------------------------
# Script: RSEM_noBW.bash
# Purpose: Estimate expression using RSEM after STAR alignment.
#          BigWig generation is handled separately by STAR_bigwig2.sh.
#
# Change log:
#   17 Feb 2026 WMD — BigWig generation moved to STAR_bigwig2.sh
#   [current]   WMD — Removed --sort-bam-by-coordinate; the coordinate-sorted
#                     BAM produced by that flag is not used by any downstream
#                     step. Sorting is already performed on the STAR genome BAM
#                     in STAR_bigwig2.sh. Removing this flag eliminates a
#                     samtools sort pass that was unnecessary and slow for
#                     high-read-count samples.
# ---------------------------------------------

# ---------------------------------------------------------------------------
# CONFIGURATION — edit these variables before running
# ---------------------------------------------------------------------------
# Root directory under which all BioProject data directories live
DATA_BASE="/path/to/your/expression/data"   # <-- replace

# Scratch filesystem base — per-project temp files go here
SCRATCH_BASE="/scratch/your/scratch/area"   # <-- replace

# Reference genome directory and files
REFdir="/path/to/your/reference/genome"     # <-- replace
REF_GTF="${REFdir}/your_genome.gtf"         # <-- replace
GENOME_FASTA="${REFdir}/your_genome.fna"    # <-- replace
# ---------------------------------------------------------------------------

# Load required software modules — adjust versions for your environment
module load rsem/1.3.3 samtools/1.20

# Retrieve environment variables exported by the parent script
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

rsemref="${scratch_dir}/rsemref"

echo "unique_name is set to ${unique_name}"
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
# Note: --sort-bam-by-coordinate is intentionally omitted. The coordinate-sorted
# genome BAM used for BigWig generation is produced by STAR_bigwig2.sh. RSEM's
# own sort output is not consumed by any downstream step in this pipeline.
echo "Estimating expression for sample: $geo_accession JOB START $(date "+%Y-%m-%d %H:%M:%S")"
echo "Running RSEM on $inbam with reference $rsemref"

rsem-calculate-expression \
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
