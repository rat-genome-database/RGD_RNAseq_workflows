#!/bin/sh
#SBATCH --job-name=RSEM72
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=05:30:00
#SBATCH --account=akwitek
#SBATCH --partition=normal
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wdemos@mcw.edu
# NOTE: Output and error logs are unified via submission script using --output and --error

# ---------------------------------------------
# Script: RSEM.sh
# Last modified: 4 Feb 2026 by WMD
# Purpose: Estimate expression using RSEM after STAR alignment
# ---------------------------------------------

# Load required software modules
module load rsem/1.3.3 samtools/1.20 deeptools/3.5.1

# Retrieve environment variables
Run=${Run}
geo_accession=${geo_accession}
PRJdir=${PRJdir}
BIOProjectID=${BIOProjectID}
Logdir=${Logdir}
FinalOPdir="${PRJdir}/${geo_accession}"
echo "PRJdir is set to $PRJdir"
echo "FinalOPdir is set to $FinalOPdir"
scratch_dir="/scratch/g/akwitek/wdemos/$BIOProjectID"
echo "Scratch working directory is set to $scratch_dir"
tempOPdir="${scratch_dir}/${Run}"
unique_name="${unique_name}"
# Reference Files (GRCr8)
REFdir="/home/wdemos/GCF_036323735.1_GRCr8/"
REF_GTF="/home/wdemos/GCF_036323735.1_GRCr8/updated_GCF_036323735.1_CM070413.1_GRCr8_genomic.gtf"
GENOME_FASTA="/home/wdemos/GCF_036323735.1_GRCr8/updated_GCF_036323735.1_CM070413.1_GRCr8_genomic.fna"
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
    --output-genome-bam \
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
echo "RSEM completed successfully for ${geo_accession}. Removing input BAM ${inbam} and copying isoform result file to transcript result file"
#rm -f "$inbam"
cp "$tempOPdir/$geo_accession.isoforms.results" $tempOPdir/$geo_accession.transcripts.results


# Check that sorted RSEM BAM exists
if [ ! -f "$scratch_dir/$Run/${sample}.genome.sorted.bam" ]; then
    echo "Error: The RSEM Genome sorted BAM file was not created for sample ${sample}."
    exit 1
fi

###echo "Generate bigWig files with command bamCoverage -b $scratch_dir/$Run/${sample}.genome.sorted.bam" -o "$FinalOPdir/RNAseq_${unique_name}.bigwig" -p ${SLURM_CPUS_ON_NOD$
###here I need other variables to properly name the output file RNAseq_${tissue}_${strain}_${sex}_${sample}.bigwig
###update 20 November 2024,normalize with BPM scale the coverage within each bin so that it accounts for differences in sequencing depth across samples.
### By normalizing by BPM, you are adjusting for sequencing depth across different samples, so that each sample's
### coverage is scaled by the number of reads (in millions)
###binSize of 50 is used because expecting paired reads of length 50-150. This size is ideal to  provide smoother coverage
###because they aggregate coverage over a #larger region, which is generally better at capturing overall transcriptional activity.
###This is less sensitive to noise in low-coverage regions and often provides a good balance between resolution and smoothing,
### for high coverage studies a bin size of 100 may be more appropriate, whereas
###if its both highcoverage and more granualar info is needed, use a smaller bin size of 10-25
echo "Generate bigWig files with command bamCoverage -b "$scratch_dir/$Run/${sample}.genome.sorted.bam" -o "$FinalOPdir/RNAseq_${unique_name}.bigwig" -p ${SLURM_CPUS_ON_NODE}--normalizeUsing BPM --binSize 50"
bamCoverage -b "$scratch_dir/$Run/${sample}.genome.sorted.bam" -o "$FinalOPdir/RNAseq_${unique_name}.bigwig" -p ${SLURM_CPUS_ON_NODE} --normalizeUsing BPM --binSize 50

if [ $? -ne 0 ]; then
    echo "Error: expression profile bigwig file generation failed for sample ${sample}"
    exit 1
fi
#Added 15 April 2025 keep the STAR Log file with the data
echo "unique_name: $unique_name"
echo "Final output directory: $FinalOPdir"
echo "$FinalOPdir/RNAseq_${unique_name}.bigwig has been generated."

# Wrap up
echo "Processing completed for sample $geo_accession on $(date "+%Y-%m-%d %H:%M:%S")"
echo "Output directory: ${tempOPdir}"

# Report elapsed time
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
echo "Total elapsed time: $(date -d@$elapsed_time -u +%H:%M:%S)"
