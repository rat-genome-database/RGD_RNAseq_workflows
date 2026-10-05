#!/bin/sh
#SBATCH --job-name=bwjson
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5000kb
#SBATCH --time=00:01:00
#SBATCH --account=your-slurm-account
#SBATCH --partition=normal
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=your@email.edu

#this script is meant to replace BWjsonV4. It has updated code to handle special characters in the AccList file that previously caused json generation issues
###############################################################################
# SAFETY SETTINGS
###############################################################################
set -euo pipefail

# Enable debugging only if DEBUG=1
[ "${DEBUG:-0}" = "1" ] && set -x

###############################################################################
# REQUIRED ENVIRONMENT VARIABLES (from wrapper / pipeline)
###############################################################################
Run=${Run}
geo_accession=${geo_accession}
BIOProjectID=${BIOProjectID}
unique_name=${unique_name}
PRJdir=${PRJdir}
scratch_dir=${scratch_dir}
Logdir=${Logdir}
baseDir=${baseDir}

# Optional metadata variables - set defaults to avoid unbound variable errors
tissue="${tissue:-}"
strain="${strain:-}"
sex="${sex:-}"
title="${title:-}"
Sample_characteristics="${Sample_characteristics:-}"
StrainInfo="${StrainInfo:-}"
PMID="${PMID:-}"
GEOpath="${GEOpath:-}"

###############################################################################
# BASIC VALIDATION
###############################################################################
for v in Run geo_accession BIOProjectID unique_name PRJdir scratch_dir baseDir; do
    eval "val=\${$v:-}"
    if [ -z "$val" ]; then
        echo "ERROR: Required variable '$v' is not set"
        exit 1
    fi
done

###############################################################################
# JSON ESCAPING (CRITICAL FIX)
###############################################################################
json_escape() {
    printf '%s' "$1" | sed \
        -e 's/\\/\\\\/g' \
        -e 's/"/\\"/g' \
        -e 's/\t/\\t/g' \
        -e 's/\r/\\r/g' \
        -e 's/\n/\\n/g'
}

###############################################################################
# MOVE TO SCRATCH
###############################################################################
cd "$scratch_dir" || {
    echo "ERROR: Cannot cd to scratch directory: $scratch_dir"
    exit 1
}

###############################################################################
# COMPUTED SEX HANDLING
###############################################################################
sex_results_file="${baseDir}/${BIOProjectID}_sex_result.txt"

get_calc_sex() {
    sample_id="$1"
    awk -F'\t' -v s="$sample_id" '
        $1 == s { print $3; found=1 }
        END { if (!found) exit 1 }
    ' "$sex_results_file"
}

CalcSex="$(get_calc_sex "$geo_accession" 2>/dev/null || echo "Unknown")"

###############################################################################
# COMPUTED STRANDEDNESS AND RNA SELECTION (from the libprep_call.sh report)
###############################################################################
prep_log="${PRJdir}/${geo_accession}/log_files/STAR/${geo_accession}_library_prep.log"

get_prep_field() {
    awk -F': *' -v k="$1" '
        $1 == k { print $2; found=1 }
        END { if (!found) exit 1 }
    ' "$prep_log"
}

CalcStrand_raw="$(get_prep_field "Strand used" 2>/dev/null || echo "")"
CalcPrep_raw="$(get_prep_field "LIBRARY_PREP" 2>/dev/null || echo "")"

case "$CalcStrand_raw" in
    FIRST_READ_TRANSCRIPTION_STRAND)  CalcStrand="Stranded (forward)" ;;
    SECOND_READ_TRANSCRIPTION_STRAND) CalcStrand="Stranded (reverse)" ;;
    NONE)                             CalcStrand="Unstranded" ;;
    *)                                CalcStrand="Unknown (sample data couldn't be measured)" ;;
esac

case "$CalcPrep_raw" in
    POLYA)         CalcPrep="Poly(A) selection" ;;
    RRNA_DEPLETED) CalcPrep="rRNA depletion" ;;
    AMBIGUOUS)     CalcPrep="Ambiguous (sample data doesn't definitively fit either group)" ;;
    *)             CalcPrep="Unknown (sample data couldn't be measured)" ;;
esac

###############################################################################
# ESCAPE ALL METADATA FOR JSON
###############################################################################
Title_esc=$(json_escape "$title")
Tissue_esc=$(json_escape "$tissue")
Strain_esc=$(json_escape "$strain")
Sex_esc=$(json_escape "$sex")
CalcSex_esc=$(json_escape "$CalcSex")
CalcStrand_esc=$(json_escape "$CalcStrand")
CalcPrep_esc=$(json_escape "$CalcPrep")
Sample_characteristics_esc=$(json_escape "$Sample_characteristics")
StrainInfo_esc=$(json_escape "$StrainInfo")
PMID_esc=$(json_escape "$PMID")
GEOpath_esc=$(json_escape "$GEOpath")
BIOProjectID_esc=$(json_escape "$BIOProjectID")
Sample_esc=$(json_escape "$geo_accession")

###############################################################################
# BIGWIG CHECK
###############################################################################
bw_file="${PRJdir}/${geo_accession}/RNAseq_${unique_name}.bigwig"

if [ ! -f "$bw_file" ]; then
    echo "ERROR: BigWig file not found: $bw_file"
    exit 1
fi

###############################################################################
# JSON GENERATION
###############################################################################
json_file="${PRJdir}/${geo_accession}/RNAseq_${unique_name}.json"

cat > "$json_file" <<EOF
{
  "type": "FeatureTrack",
  "trackId": "RNAseq_${unique_name}",
  "name": "RNAseq_${unique_name}",
  "category": [
    "RNA-Seq",
    "${Tissue_esc}",
    "${Strain_esc}"
  ],
  "assemblyNames": [
    "GRCr8"
  ],
  "metadata": {
    "Sample Characteristic": "${Sample_characteristics_esc}",
    "Tissue": "${Tissue_esc}",
    "Strain": "${Strain_esc}",
    "RGD Strain Report": "${StrainInfo_esc}",
    "Sex": "${Sex_esc}",
    "Computed Sex": "${CalcSex_esc}",
    "Computed Strandedness": "${CalcStrand_esc}",
    "Computed RNA Selection": "${CalcPrep_esc}",
    "RGD Metadata Report": "https://rgd.mcw.edu/rgdweb/report/expressionStudy/main.html?geoAcc=${BIOProjectID_esc}",
    "Project Title": "${Title_esc}",
    "Project Repository Link": "${GEOpath_esc}",
    "Project Accession ID": "${BIOProjectID_esc}",
    "Sample Accession ID": "${Sample_esc}",
    "PubMed ID": "PMID:${PMID_esc}",
    "Data Processing": "HPC RGD workflow",
    "Read alignment": "STAR v2.7.10b",
    "Genome version": "GCF_036323735.1 GRCr8",
    "Expression Quantification": "RSEM v1.3.1"
  },
  "adapter": {
    "type": "BigWigAdapter",
    "bigWigLocation": {
      "locationType": "UriLocation",
      "uri": "RNAseq_${unique_name}.bigwig"
    }
  },
  "displays": [
    {
      "type": "LinearWiggleDisplay",
      "displayId": "RNAseq_${unique_name}-LinearWiggleDisplay"
    }
  ]
}
EOF

###############################################################################
# JSON VALIDATION (FAIL FAST)
###############################################################################
# Make python3 available on compute nodes; adjust the module name for your cluster
module load python/3.9.1 2>/dev/null || module load python3 2>/dev/null || true
python3 - <<EOF
import json
with open("$json_file") as f:
    json.load(f)
EOF

###############################################################################
# FINAL STATUS
###############################################################################
echo "SUCCESS: JSON generated and validated:"
echo "  $json_file"
