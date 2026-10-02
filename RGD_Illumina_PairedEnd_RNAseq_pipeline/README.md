# RGD Bulk RNA-seq Pipeline (GRCr8)

A SLURM-based pipeline for processing paired-end bulk RNA-seq data from NCBI SRA. Designed for HPC clusters with SLURM job scheduling. Produces STAR-aligned BAM files, BigWig coverage tracks, RSEM expression matrices, sex estimation, library strandedness and library preparation (poly(A) selection vs rRNA depletion) inferred from the data, and JBrowse2 visualization sessions.

**Genome:** GCF_036323735.1 GRCr8 (rat)  
**Developed at:** Rat Genome Database (RGD), Medical College of Wisconsin

---

## Pipeline Overview

```
Input: Accession list (TSV) + BioProject ID
       │
       ▼
STEP 1: SRA Download & QC
  ├── prefetch (up to 8 attempts with exponential backoff)
  ├── vdb-validate
  ├── fasterq-dump → paired FASTQ (3 retries; exits code 2 if single-end detected)
  ├── FastQC + FastQ-Screen (per sample)
  └── MultiQC report
       │
       ▼
STEP 2: Alignment, Quantification & Visualization
  ├── STAR genome index (GRCr8)
  ├── STAR alignment (with GeneCounts) → sorted BAM + BigWig (BPM normalized)
  ├── Strand detection from STAR gene counts → Picard CollectRnaSeqMetrics
  ├── Library prep call: poly(A) vs rRNA-depleted (libprep_call.sh, rule LP2)
  ├── STARQC alignment rate summary → PASS/FAIL filter (threshold: <50% unmapped)
  ├── Sex estimation (chrX/Y read depth ratio)
  ├── RSEM reference + per-sample quantification (PASS samples only)
  ├── RSEM expression matrices (TPM + counts, genes + transcripts)
  ├── MultiQC final report
  ├── Sex and library prep report: final sex call from sex-gene TPMs (Xist vs
  │     median of Uty, Ddx3y, Kdm5d, Eif2s3y), strand, library prep, study-level call
  ├── BigWig JSON track configs (for JBrowse2, PASS samples only)
  └── JBrowse2 session JSON (PASS samples only)
```

---

## Requirements

### Software (via `module load` on your cluster)

| Tool | Version used |
|------|-------------|
| STAR | 2.7.10b |
| RSEM | 1.3.3 |
| samtools | 1.20 |
| deeptools | 3.5.1 |
| Picard (Java) | 2.25.0 |
| SRA Toolkit | 3.1.1 |
| FastQC | 0.11.9 |
| FastQ-Screen | 0.15.2 |
| Bowtie2 | 2.5 |
| MultiQC | 1.18 |
| Python | 3.9+ |
| Perl | 5.x |
| pigz | any |

### Reference Files (GRCr8)
Download from NCBI accession `GCF_036323735.1`:
- Genome FASTA: `GCF_036323735.1_GRCr8_genomic.fna`
- Genome GTF: `GCF_036323735.1_GRCr8_genomic.gtf`

Also needed, built from the same GTF:
- rRNA interval list and refFlat annotation, used by Picard CollectRnaSeqMetrics (`RRNA_INTERVALS`, `REF_FLAT` in `STAR_bigwig3.sh`); see [Picard reference files](#picard-reference-files-collectrnaseqmetrics) below
- Non-polyadenylated small-RNA gene list for the library prep call. Build it once per annotation:
  ```bash
  bash scripts/build_nonpolyA_genes.sh /path/to/GRCr8.gtf /path/to/GRCr8_nonpolyA_genes.txt
  ```
  The list records its source GTF and MD5 checksum and is made read-only; every library prep call records the list's checksum.

### Picard reference files (CollectRnaSeqMetrics)

`STAR_bigwig3.sh` runs Picard `CollectRnaSeqMetrics` (module `picard/2.25.0`, called as
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
- chrMT in this reference has no RefSeq gene models; the GTF carries a single `MT_mask`
  feature instead (6-22% of assigned counts in our data). If Picard's intergenic fraction
  looks inflated, check `grep -c MT_mask $gtf.refflat`.
- A refFlat that gives `PCT_MRNA_BASES = 0`, `PCT_INTRONIC_BASES = 0` and
  `PCT_INTERGENIC_BASES = 1` for every sample is a broken reference, not a library-prep
  result.
- Picard's intronic fraction did not separate polyA from rRNA-depleted libraries across
  our control studies; the library-prep call uses non-polyadenylated small-RNA reads from
  STAR `ReadsPerGene.out.tab` (rule LP2, `libprep_call.sh`), with Picard output kept as
  supporting evidence.

### FastQ-Screen Genomes
Download the standard genome set using FastQ-Screen's built-in tool:
```bash
fastq_screen --get_genomes --outdir /path/to/your/FastQ_Screen_Genomes
```
Then update the paths in `dependencies/fastq_screen.conf`.

---

## Repository Structure

```
RGD_Illumina_PairedEnd_RNAseq_pipeline/
├── scripts/                  Pipeline scripts (see table below)
├── dependencies/
│   ├── rsem-generate-data-matrix        Custom Perl script (outputs TPM matrix)
│   ├── rsem-generate-data-matrix-counts Custom Perl script (outputs counts matrix)
│   └── fastq_screen.conf               FastQ-Screen config template
├── docs/
│   ├── example_AccList.txt             Example accession list format
│   └── example_project_list.txt        Example orchestrator project list
├── README.md
├── CHANGELOG.md
├── CONTRIBUTING.md
├── LICENSE
└── .gitignore
```

---

## Configuration

Before running, update the `USER CONFIGURATION` block at the top of each script. The variables to set are:

| Variable | Description |
|----------|-------------|
| `SCRIPT_DIR` | Full path to the `scripts/` directory |
| `myDir` | Your home/base directory |
| `SCRATCH_BASE` | Scratch area that holds the per-project scratch folders (each script uses `SCRATCH_BASE/<BIOProjectID>`). **Set the same value in every script.** |
| `REF_GTF` | Path to GRCr8 GTF file |
| `GENOME_FASTA` | Path to GRCr8 genome FASTA |
| `RRNA_INTERVALS`, `REF_FLAT` | Picard rRNA interval list and refFlat file (`STAR_bigwig3.sh`) |
| `NONPOLYA_GENES` | Non-polyA gene list from `build_nonpolyA_genes.sh` (default in `libprep_call.sh`, or export it) |
| `--account` | Your SLURM account name (in `#SBATCH` headers) |
| `--mail-user` | Your email for SLURM failure notifications |

Also update `dependencies/fastq_screen.conf` with paths to your local FastQ-Screen genome indices.

---

## Input Format

Tab-delimited accession list with a header row:

```
Run	geo_accession	Tissue	Strain	Sex	PMID	GEOpath	Title	Sample_characteristics	StrainInfo
SRR12345678	GSM1234567	Liver	BN/NHsdMcwi	M	12345678	https://...	Study Title	age: 12 weeks	https://...
```

Multiple SRA runs per GSM sample are supported — they are automatically grouped and merged for alignment.

---

## Usage

### Recommended: Batch Mode (Orchestrator)

```bash
# project_list.txt: <acclist_path>  <BioProjectID>  <read_length>
bash scripts/bulk_orchestrator_production_diskGuard.bash project_list.txt
```

Run it on the login node inside `screen` or `tmux`; it submits STEP 1 and STEP 2 as SLURM jobs and waits on them.

The orchestrator enforces concurrent project limits, disk space guards, and handles large vs. small project scheduling automatically.

### Single Project

```bash
# Step 1: Download and QC
sbatch scripts/run_SRA2QC_diskGuard.bash /path/to/AccList.txt PRJNA123456

# Step 2: Alignment, quantification, visualization
sbatch scripts/run_RNApipeline_pairedG8_diskGuard.bash /path/to/AccList.txt PRJNA123456 150
```

### Utility Scripts

```bash
# Count unique samples in an accession list
bash scripts/sample_counting.sh /path/to/AccList.txt

# Re-call library prep for an already-processed project (no realignment), then refresh its report
bash scripts/libprep_recall.sh PRJNA123456
bash scripts/ConflictedSampleReport_v8.sh PRJNA123456
```

Submit Step 2 from the `scripts/` directory: `STAR_bigwig3.sh` finds `libprep_call.sh` in the directory the job was submitted from (or set `--export=LIBPREP_CALL=/path/to/libprep_call.sh`).

---

## Script Reference

| Script | Role |
|--------|------|
| `bulk_orchestrator_production_diskGuard.bash` | Batch orchestrator with hybrid scheduling and disk guard |
| `run_SRA2QC_diskGuard.bash` | Step 1 controller: SRA download and QC |
| `run_RNApipeline_pairedG8_diskGuard.bash` | Step 2 controller: alignment through visualization |
| `SRA2QC_production.sh` | Per-sample SLURM job: SRA → FASTQ → FastQC/FastQ-Screen |
| `starRef_v4.sh` | Generates STAR genome index |
| `STAR_bigwig3.sh` | STAR alignment (with GeneCounts) → sorted BAM → BigWig (BPM); strand detection; Picard CollectRnaSeqMetrics; calls `libprep_call.sh` |
| `libprep_call.sh` | Per-sample library prep call (POLYA / RRNA_DEPLETED / AMBIGUOUS / UNKNOWN) with reason, rule version and QC flag |
| `libprep_recall.sh` | Utility: re-call library prep for an already-processed project without realignment |
| `build_nonpolyA_genes.sh` | One-time setup: builds the non-polyA small-RNA gene list from the GTF |
| `pSTARQC_v1.sh` | Parses STAR logs; generates PASS/FAIL alignment summary |
| `ComputeSex_v5.sh` | Estimates sex from chrX/Y read depth ratio |
| `ConflictedSampleReport_v8.sh` | Final sex call from sex-linked gene TPMs (overrides the chrX/Y call); adds strand, library prep metrics and calls; writes the per-study library prep summary |
| `RSEMref_v4.sh` | Generates RSEM reference |
| `RSEM_noBW.bash` | Per-sample RSEM expression quantification |
| `RSEMmatrix_v5.sh` | Combines RSEM results into project-level matrices + final MultiQC |
| `BWjson_v7.sh` | Generates per-sample JBrowse2 BigWig track JSON |
| `JBrowseSession_v2.sh` | SLURM wrapper for the JBrowse session builder; passes the STARQC PASS list |
| `make_jbrowse_session_for_bioproject.py` | Builds combined JBrowse2 session JSON with color grouping (PASS samples only) |
| `sample_counting.sh` | Utility: count unique samples in an accession list |
| `lib_v10.sh` | Shared helper functions (logging, job submission, status checking) |

---

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
│   ├── STARQC/<BIOProjectID>_STAR_Align_sum.txt
│   └── ...
├── <BIOProjectID>.genes.TPM.matrix
├── <BIOProjectID>.genes.counts.matrix
├── <BIOProjectID>.transcripts.TPM.matrix
├── <BIOProjectID>.transcripts.counts.matrix
├── <BIOProjectID>_sex_result.txt
├── <BIOProjectID>_sex_conflict_report.txt      sex call, strand, library prep per sample
├── <BIOProjectID>_library_prep_summary.txt     per-study call counts, metric medians, StudyCall
├── <BIOProjectID>_fastq_multiQC_report.html
├── <BIOProjectID>_final_multiQC_report.html
└── <BIOProjectID>_jbrowse_session_GRCr8.json
```

> **Samples that fail STAR alignment QC:** a BigWig is still generated and kept in the sample's
> `reads_fastq/<GSM_ID>/` folder, so the alignment remains available for download. These samples
> have no track JSON and are not included in the JBrowse session, the TPM/count matrices, or the
> sex call. Each sample's PASS/FAIL status is listed in `log_files/STARQC/<BIOProjectID>_STAR_Align_sum.txt`.

---

## Custom Dependencies

The `dependencies/` directory contains two custom Perl scripts derived from RSEM that output **TPM values** and **expected counts** respectively, with filenames extracted using `File::Basename` (added August 2024) for cleaner matrix column headers.

---

## Citation

If you use this pipeline, please cite:
- Dobin et al. (2013) STAR. *Bioinformatics* 29(1):15–21
- Li & Dewey (2011) RSEM. *BMC Bioinformatics* 12:323
- Andrews S. (2010) FastQC. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- Ewels et al. (2016) MultiQC. *Bioinformatics* 32(19):3047–3048
- Ramírez et al. (2016) deepTools2. *Nucleic Acids Research* 44(W1):W160–W165

---

## License

Pipeline scripts: MIT — see [LICENSE](LICENSE)  
RSEM-derived scripts in `dependencies/`: GPL-3.0 — see [dependencies/LICENSE](dependencies/LICENSE)
