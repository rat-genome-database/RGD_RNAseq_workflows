# RGD Illumina Single-End RNAseq Pipeline

Workflow for processing GEO single-end Illumina RNA-seq data using STAR + RSEM on a SLURM HPC cluster.

## Overview

This pipeline downloads SRA data, performs QC, aligns reads with STAR, and quantifies expression with RSEM. It is designed for bulk RNA-seq datasets from NCBI GEO.

For each sample it also infers library strandedness and library preparation (poly(A) selection vs. rRNA depletion) from the data, and calls sample sex from sex-specific gene expression. See `CHANGELOG.md` for version history.

### Pipeline Steps

| Step | Script | Description |
|------|--------|-------------|
| 1 | `run_SRA2QC_SE_v1.bash` | Download SRA reads and run FastQC / FastQ Screen |
| 2 | `run_RNApipeline_SE_diskGuard_v2.bash` | STAR alignment, strand and library prep calls, RSEM quantification, reports and JBrowse files |
| — | `bulk_orchestrator_production_diskGuard.bash` | Orchestrates both steps across multiple projects |

### Supporting Scripts

| Script | Description |
|--------|-------------|
| `SRA2QC_SE_v1.sh` | Core SRA download + QC logic |
| `STAR_SE_v2.sh` | STAR alignment (single-end, with gene counts) → sorted BAM → BigWig; strand detection; Picard CollectRnaSeqMetrics; calls `libprep_call.sh` |
| `libprep_call.sh` | Per-sample library prep call (POLYA / RRNA_DEPLETED / AMBIGUOUS / UNKNOWN) with reason, rule version and QC flag |
| `libprep_recall.sh` | Utility: re-call library prep for an already-processed project without realignment |
| `build_nonpolyA_genes.sh` | One-time setup: builds the gene list used by the library prep call from the GTF |
| `RSEM_SE_v1.sh` | RSEM quantification (single-end) |
| `RSEMref_v4.sh` | Build RSEM reference |
| `RSEMmatrix_v7.sh` | Generate count/TPM expression matrices, final MultiQC report, and the sample report |
| `starRef_v4.sh` | Build STAR genome index |
| `pSTARQC_v1.sh` | Parse STAR alignment QC logs |
| `ComputeSex_v6.sh` | Estimate sex from chrX/chrY read depth ratio |
| `ConflictedSampleReport_v7.sh` | Final sex call from sex-specific gene TPMs (overrides the chrX/chrY call); adds strand, library prep metrics and calls; writes the per-study library prep summary |
| `BWjson_v7.sh` | Generate BigWig track JSON for JBrowse (includes computed sex, strandedness and library capture) |
| `JBrowseSession_v2.sh` | Create JBrowse2 session files (STAR QC PASS samples only) |
| `make_jbrowse_session_for_bioproject.py` | Python helper for JBrowse session generation |
| `lib_v10.sh` | Shared library functions |
| `sample_counting.sh` | Count samples per project |

## Requirements

- SLURM workload manager
- [SRA Toolkit](https://github.com/ncbi/sra-tools)
- [STAR](https://github.com/alexdobin/STAR)
- [RSEM](https://github.com/deweylab/RSEM)
- [samtools](https://www.htslib.org/)
- [Picard](https://broadinstitute.github.io/picard/) (CollectRnaSeqMetrics; version 2.25.0 used)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/)
- [deepTools](https://deeptools.readthedocs.io/) (for BigWig generation)
- Python 3

## Configuration

Before running, update the `USER CONFIGURATION` block at the top of each script:

```bash
SCRIPT_DIR="/path/to/RGD_Illumina_SingleEnd_RNAseq_pipeline"   # Location of pipeline scripts
myDir="/path/to/home"                                            # Your home/base directory
SCRATCH_BASE="/path/to/scratch"                                  # Scratch filesystem for temp files
REF_GTF="/path/to/genome/GCF_036323735.1_GRCr8/..."            # GTF annotation file
GENOME_FASTA="/path/to/genome/GCF_036323735.1_GRCr8/..."       # Genome FASTA file
RRNA_INTERVALS=".../GRCr8_rRNA.intervals"                        # STAR_SE_v2.sh: Picard rRNA interval list
REF_FLAT=".../updated_GCF_036323735.1_CM070413.1_GRCr8_genomic.gtf.refflat"   # STAR_SE_v2.sh: Picard refFlat
```

In `libprep_call.sh`, set the default `NONPOLYA_GENES` path to the gene list built by `build_nonpolyA_genes.sh` (or export `NONPOLYA_GENES` before running).

Also update SLURM directives at the top of each script:
```bash
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --mail-user=your@email.edu
```

## Reference Files for Strand and Library Prep Calls

Three files are built once from the same GTF used for the STAR index.

### Library prep gene list

```bash
bash build_nonpolyA_genes.sh /path/to/GRCr8.gtf /path/to/GRCr8_nonpolyA_genes.txt
```

The list (snoRNA, snRNA, scaRNA, SRP RNA, RNase MRP RNA and telomerase RNA genes, which lack stable poly(A) tails) records its source GTF and MD5 checksum and is made read-only; every library prep call records the list's checksum.

### Picard reference files (CollectRnaSeqMetrics)

`STAR_SE_v2.sh` runs Picard `CollectRnaSeqMetrics` (module `picard/2.25.0`, called as
`java -jar "$PICARD"`) on each sample's sorted genome BAM. Picard needs two reference
files, both built from the same GTF used to build the STAR index:

| File | Picard argument |
|---|---|
| `updated_GCF_036323735.1_CM070413.1_GRCr8_genomic.gtf.refflat` | `REF_FLAT` |
| `GRCr8_rRNA.intervals` | `RIBOSOMAL_INTERVALS` |

Both live in your reference directory.
Chromosome names in both must match the BAM (`chr1` ... `chrMT`, the STAR index naming).

Requirements: UCSC `gtfToGenePred`, `samtools`, `awk`.

```bash
cd /path/to/your/GRCr8/reference
gtf=updated_GCF_036323735.1_CM070413.1_GRCr8_genomic.gtf
```

#### 1. refFlat

```bash
gtfToGenePred -genePredExt -ignoreGroupsWithoutExons $gtf tmp.genePred
awk 'BEGIN{OFS="\t"} {print $12,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' tmp.genePred > $gtf.refflat
rm tmp.genePred
```

- `-ignoreGroupsWithoutExons` is required for NCBI GTFs: their `gene` feature lines have an
  empty transcript ID and no exons, and `gtfToGenePred` otherwise stops with
  `no exons defined for group , feature gene`. Skipping them loses no transcripts.
- Do **not** add `-geneNameAsName2`: NCBI GTFs keep the gene symbol in `gene_id`, and that
  flag can leave the gene-name column blank.
- The awk step puts genePred column 12 (gene name) first, which is the refFlat layout:
  gene name, transcript, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount,
  exonStarts, exonEnds. Printing column **13** instead (the CDS-status flag) produces a
  file that passes a simple 11-column check but has `cmpl` / `incmpl` / `none` / `unk` as
  gene names; Picard then loads ~4 "genes" and reports every base as intergenic.

Check:

```bash
head -3 $gtf.refflat                       # column 1 = gene symbols (e.g. LOC134485287)
cut -f1 $gtf.refflat | sort -u | wc -l     # ~38,673 for this GTF
awk -F'\t' '
    NF != 11 { print "Error: line " NR " does not have 11 columns"; err=1 }
    $1 == "" || $2 == "" { print "Error: line " NR " missing gene or transcript ID"; err=1 }
    $1 ~ /^(cmpl|incmpl|none|unk)$/ { flag++ }
    { genes[$1] = 1 }
    END {
        for (g in genes) n++
        print n " unique gene names in column 1"
        if (flag) { print "Error: " flag " lines have a CDS status flag in column 1"; err=1 }
        if (!err) print "Success: 11 columns, gene names look valid."
    }' $gtf.refflat
```

Picard's log should then show `Loaded 38673 genes` (STAR's `ReadsPerGene.out.tab` has 38,674).

#### 2. rRNA interval list

A Picard interval list is a SAM-style header (`@HD` + `@SQ` lines) followed by
tab-separated `chrom  start  end  strand  name` rows (1-based, closed). The header must
match the BAM's sequence names **and lengths** exactly, so it is copied from a BAM
aligned against this STAR index (any sample's `_GENOME_SORT.bam` works).

```bash
grep -c 'gene_biotype "rRNA"' $gtf        # 159 in this GTF; must be > 0

bam=<any sample>_GENOME_SORT.bam
samtools view -H $bam | grep -E '^@(HD|SQ)' > GRCr8_rRNA.intervals
awk -F'\t' 'BEGIN{OFS="\t"} $3=="gene" && $9 ~ /gene_biotype "rRNA"/ {
    match($9, /gene_id "[^"]+"/)
    print $1, $4, $5, $7, substr($9, RSTART+9, RLENGTH-10) }' $gtf >> GRCr8_rRNA.intervals
```

Check:

```bash
grep -vc '^@' GRCr8_rRNA.intervals                                  # 159, same as the GTF count
grep -v '^@' GRCr8_rRNA.intervals | cut -f1 | sort | uniq -c        # chrMT should be present
grep -v '^@' GRCr8_rRNA.intervals | grep -Ei 'rn45s|rn18s|rn28s'     # nuclear rRNAs annotated
```

Rn45s, Rn18s, Rn28s, Rn5-8s and Rn5s are annotated in this GTF, so the major nuclear
rRNAs are covered. If Picard logs `The RIBOSOMAL_INTERVALS file ... does not contain
intervals`, the file has a header but no usable rows; rebuild it as above.

#### Notes

- If the BAM header changes (new STAR index, different assembly), rebuild
  `GRCr8_rRNA.intervals` from a BAM made with the new index.
- A refFlat that gives `PCT_MRNA_BASES = 0`, `PCT_INTRONIC_BASES = 0` and
  `PCT_INTERGENIC_BASES = 1` for every sample is a broken reference, not a library-prep
  result.
- The library prep call uses reads on the gene list above from STAR `ReadsPerGene.out.tab`
  (rule LP2, `libprep_call.sh`); Picard output is kept as supporting evidence.

## Input Files

Place the following in your working directory before running:

- **AccList file**: Tab-separated, one row per SRA run, with a header line: `Run  geo_accession  Tissue  Strain  Sex  PMID  GEOpath  Title  Sample_characteristics  StrainInfo` (see `docs/example_AccList.txt`)
- **Project list file** (orchestrator only): one project per line, `<path_to_AccList>  <BIOProjectID>  <read_length>`; lines starting with `#` are ignored (see `docs/example_project_list.txt`)

## Usage

**Single project (two-step):**
```bash
# Step 1: Download and QC
sbatch run_SRA2QC_SE_v1.bash <AccList_file> <BIOProjectID>

# Step 2: Align and quantify (read length in bp; STAR uses sjdbOverhang = read length - 1)
sbatch run_RNApipeline_SE_diskGuard_v2.bash <AccList_file> <BIOProjectID> <read_length>
```

Submit Step 2 from the pipeline directory: `STAR_SE_v2.sh` finds `libprep_call.sh` in the directory the job was submitted from (or set `--export=LIBPREP_CALL=/path/to/libprep_call.sh`).

**Multiple projects (orchestrated):**
```bash
bash bulk_orchestrator_production_diskGuard.bash project_list.txt
```

**Re-call library prep for an already-processed project** (no realignment; requires the STAR gene counts and Picard metrics written by `STAR_SE_v2.sh`), then refresh its report:
```bash
bash libprep_recall.sh <BIOProjectID>
baseDir=/path/to/data/expression/GEO/<BIOProjectID> bash ConflictedSampleReport_v7.sh <BIOProjectID>
```

## Output Files

All permanent outputs go to `<myDir>/data/expression/GEO/<BIOProjectID>/`:

```
<BIOProjectID>/
├── reads_fastq/<GSM_ID>/
│   ├── RNAseq_<unique_name>.bigwig
│   ├── RNAseq_<unique_name>.json
│   └── log_files/
│       └── STAR/
│           ├── <GSM_ID>_STARReadsPerGene.out.tab   STAR gene counts
│           ├── <GSM_ID>_strand_detection.log       inferred strand
│           ├── <GSM_ID>_rna_metrics.txt            Picard CollectRnaSeqMetrics
│           └── <GSM_ID>_library_prep.log           library prep call, reason, rule version, QC flag
├── log_files/
├── <BIOProjectID>.genes.TPM.matrix
├── <BIOProjectID>.genes.counts.matrix
├── <BIOProjectID>.transcripts.TPM.matrix
├── <BIOProjectID>.transcripts.counts.matrix
├── <BIOProjectID>_sex_result.txt
├── <BIOProjectID>_sex_conflict_report.txt      workflow version (line 1), sex call, strand, library prep per sample
├── <BIOProjectID>_library_prep_summary.txt     per-study call counts, metric medians, StudyCall
├── <BIOProjectID>_final_multiQC_report.html
└── <BIOProjectID>_jbrowse_session_GRCr8.json
```

Samples that fail STAR alignment QC keep their BigWig in `reads_fastq/<GSM_ID>/` but have no track JSON and are not included in the JBrowse session, the TPM/count matrices, or the sex call.

## Reference Genome

Scripts are configured for **GRCr8** (rat genome GCF_036323735.1). Commented-out lines for **mRatBN7.2** are retained for reference.

## License

GPL-3.0 — see `dependencies/LICENSE`
