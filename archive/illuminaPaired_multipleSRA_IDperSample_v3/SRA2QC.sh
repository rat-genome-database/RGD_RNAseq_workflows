#!/bin/bash
#SBATCH --job-name=SRA2QC
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=1gb
#SBATCH --time=01:35:00
#SBATCH --account=akwitek
#SBATCH --partition=normal
#SBATCH --output=%x-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wdemos@mcw.edu
#SBATCH --error=SRA2QC.err
# Enable debugging by printing each command before execution
set -x
# Load required modules
module load sratoolkit/3.1.1 fastq-screen/0.15.2 fastqc/0.11.9 bowtie/2.5 python/3.9.1


# Retrieve environment variables
Run=${Run}
geo_accession=${geo_accession}
PRJdir=${PRJdir}
screenconfig=${screenconfig}
BIOProjectID=${BIOProjectID}
echo "Run: $Run"
echo "geo_accession: $geo_accession"
echo "PRJdir: $PRJdir"
echo "screenconfig: $screenconfig"
#scratch_dir="/scratch/g/akwitek/wdemos/$Run"
scratch_dir="/scratch/g/akwitek/wdemos/$BIOProjectID"

# Create a scratch directory for the run
mkdir -p "$scratch_dir"
echo "Processing sample: $Run"

# Prefetch data
prefetch -O "$scratch_dir" "$Run" || { echo "Prefetch failed for $Run"; exit 1; }

# Validate data
vdb-validate "$scratch_dir/$Run" || { echo "Validation failed for $Run"; exit 1; }

# Attempt to run fasterq-dump up to 3 times if needed
attempt=1
max_attempts=3

while [ $attempt -le $max_attempts ]; do
    echo "Attempt $attempt: Running fasterq-dump for $Run"
    if fasterq-dump.3.1.1 "$scratch_dir/$Run" --threads 12 --outdir "$scratch_dir/$Run" -t /scratch/g/akwitek/wdemos/$BIOProjectID -x; then
        fastq1="${scratch_dir}/${Run}/${Run}_1.fastq"
        fastq2="${scratch_dir}/${Run}/${Run}_2.fastq"
        fastq3="${scratch_dir}/${Run}/${Run}.fastq"

        if [ -f "$fastq1" ] && [ -f "$fastq2" ]; then 
            echo "Renaming fastq files to GEO accession: $geo_accession"
            new_fastq1="${scratch_dir}/${Run}/${geo_accession}_${Run}_1.fastq"
            new_fastq2="${scratch_dir}/${Run}/${geo_accession}_${Run}_2.fastq"
            mv "$fastq1" "$new_fastq1"
            mv "$fastq2" "$new_fastq2"

            if [ -f "$fastq3" ]; then
                new_fastq3="${scratch_dir}/${Run}/${geo_accession}_${Run}.fastq"
                mv "$fastq3" "$new_fastq3"
                echo "Renamed $fastq3 to $new_fastq3"
            fi

            echo " Run FastQC and FastQ-Screen in the scratch directory"
            fastqc "$new_fastq1" "$new_fastq2" --outdir "$scratch_dir/${Run}/" || { echo "FastQC failed"; exit 1; }
            fastq_screen --outdir "$scratch_dir/${Run}/" -conf "$screenconfig" "$new_fastq1" "$new_fastq2" || { echo "FastQ-Screen failed"; exit 1; }

            echo " Compress fastq files"
            pigz "$new_fastq1" "$new_fastq2"
            echo "Processed and compressed $geo_accession"

            # If R3 exists, compress it
            if [ -f "$new_fastq3" ]; then
                echo "Compressing FASTQ files for unmated reads for sample $geo_accession START: $(date "+%Y-%m-%d %H:%M:%S")"
                pigz "$new_fastq3"
            fi

            #echo "Moving results to PRJdir"
            #cp "${scratch_dir}/${Run}/${geo_accession}*.fastq.gz" "$PRJdir/"
            #cp "${scratch_dir}/${Run}/*fastqc.html" "$PRJdir/"
            #cp "${scratch_dir}/${Run}/*screen.html" "$PRJdir/"
            #echo "Moved results to $PRJdir"

             echo "Remove .SRA file"
             sra_file="${scratch_dir}/${Run}/${Run}.sra"
             echo "Run directory: ${scratch_dir}/${Run}"
             echo "SRA file: $sra_file "
             if [ -f "$sra_file" ]; then
                 rm "$sra_file"
             fi

            exit 0
        else
            echo "FastQ files not found for $Run on attempt $attempt"
        fi
    else
        echo "fasterq-dump failed for $Run on attempt $attempt"
    fi
    attempt=$((attempt + 1))
    sleep 10  # Optional: wait before retrying
done

echo "Failed to process $Run after $max_attempts attempts"

