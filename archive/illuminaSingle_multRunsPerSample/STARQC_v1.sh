#!/usr/bin/env bash
#SBATCH --job-name=STARsum
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=512mb
#SBATCH --time=00:05:00
#SBATCH --account=akwitek
#SBATCH --partition=normal
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wdemos@mcw.edu

# Usage: STAR_unmapped_summary.sh <UniqueAccList.tsv> <BioProject_ID>
# UniqueAccList.tsv must keep the header and have geo_accession in column 2.

set -uo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <UniqueAccList.tsv> <BioProject_ID>" >&2
  exit 1
fi

uniqueAccList="$1"
BIOProjectID="$2"

# Match your workflow paths
myDir="/home/wdemos"
baseDir="${myDir}/data/expression/GEO/${BIOProjectID}"
Logdir="${baseDir}/log_files/STARQC"
scratch_dir="/scratch/g/akwitek/wdemos/${BIOProjectID}"

mkdir -p "$baseDir" "$Logdir"

#OUT="${baseDir}/${BIOProjectID}_STAR_Align_sum.txt"
OUT="${Logdir}/${BIOProjectID}_STAR_Align_sum.txt"
TMP="${OUT}.tmp"

# Helper: extract numeric field from STAR Log.final.out
# STAR writes lines like: "Number of input reads | <tab> 24776293"
get_value () {
  local pattern="$1" file="$2"
  awk -F '|' -v p="$pattern" '
    $1 ~ p {
      v=$2
      gsub(/^[ \t]+|[ \t]+$/, "", v)
      gsub(/,/, "", v)
      print v
      exit
    }
  ' "$file"
}

echo -e "SampleID\tinput_reads\tunaligned_reads\tUnmapped_Rate\tStatus" > "$TMP"

# Iterate deduplicated list; skip header
tail -n +2 "$uniqueAccList" | while IFS=$'\t' read -r _run geo_accession _rest; do
  [[ -z "$geo_accession" ]] && continue

  wrkDir="${scratch_dir}/${geo_accession}"
  # With outFileNamePrefix="${wrkDir}/${geo_accession}_STAR", STAR writes "${prefix}Log.final.out"
  star_log="${wrkDir}/${geo_accession}_STARLog.final.out"

  if [[ ! -f "$star_log" ]]; then
    echo -e "${geo_accession}\tNA\tNA\tNA\tNO_LOG" >> "$TMP"
    continue
  fi

  input_reads="$(get_value 'Number of input reads' "$star_log")"
  unm_mm="$(get_value 'Number of reads unmapped: too many mismatches' "$star_log")"
  unm_short="$(get_value 'Number of reads unmapped: too short' "$star_log")"
  unm_other="$(get_value 'Number of reads unmapped: other' "$star_log")"

  # Validate input_reads
  if [[ -z "$input_reads" || "$input_reads" == "0" ]]; then
    echo -e "${geo_accession}\t${input_reads:-NA}\tNA\tNA\tINVALID_LOG" >> "$TMP"
    continue
  fi

  # Totals and percent (use awk for safe float math)
  unmapped_total="$(awk -v a="$unm_mm" -v b="$unm_short" -v c="$unm_other" 'BEGIN{printf "%.0f", a+b+c}')"
  unmapped_pct="$(awk -v u="$unmapped_total" -v n="$input_reads" 'BEGIN{printf "%.2f", (u/n)*100.0}')"

  # PASS if mapped > 50% => unmapped < 50%
  status="$(awk -v p="$unmapped_pct" 'BEGIN{print (p<50.0) ? "PASS" : "FAIL"}')"

  echo -e "${geo_accession}\t${input_reads}\t${unmapped_total}\t${unmapped_pct}\t${status}" >> "$TMP"
done

mv "$TMP" "$OUT"
echo "STAR alignment summary written to: $OUT"
#Copy to Logdir:
# cp -f "$OUT" "$Logdir/"
