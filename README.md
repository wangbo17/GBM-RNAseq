# GBM RNA-Seq Processing Pipeline

<pre>
=============================================================================================
  
    ____ ____  __  __   ____        _ _      ____  _   _    _        ____             
   / ___| __ )|  \/  | | __ ) _   _| | | __ |  _ \| \ | |  / \      / ___|  ___  __ _ 
  | |  _|  _ \| |\/| | |  _ \| | | | | |/ / | |_) |  \| | / _ \ ____\___ \ / _ \/ _` |
  | |_| | |_) | |  | | | |_) | |_| | |   <  |  _ <| |\  |/ ___ |________) |  __| (_| |
   \____|____/|_|  |_| |____/ \__,_|_|_|\_\ |_| \_|_| \_/_/   \_\   |____/ \___|\__, |
                                                                                   |_|
      
                                                             Author: Bo Wang | Version: Beta
=============================================================================================
</pre>

## Overview
GBM Bulk RNA-Seq is a bioinformatics pipeline designed for the analysis of glioblastoma (GBM) bulk RNA sequencing data. This pipeline performs comprehensive preprocessing, alignment, quantification, and quality control of RNA-seq data using modular Nextflow processes. It supports both paired-end and single-end reads and is containerized for reproducibility. The workflow integrates industry-standard tools including **FastQC**, **Trim Galore**, **SortMeRNA**, **STAR**, **RSEM**, and **RSeQC**, with summary reporting via **MultiQC**.

---

## Step 1: Raw Read Quality Assessment using FastQC

Raw sequencing data quality was assessed using **FastQC (v0.12.1)** to ensure the reliability of downstream analyses. For each sample, FastQC was executed on the untrimmed FASTQ files (both forward and reverse reads for paired-end libraries, or single files for single-end libraries). The process generates comprehensive quality control (QC) reports in both HTML and compressed formats.

The FastQC report evaluates multiple key metrics, including:
- **Per-base sequence quality**: to detect low-quality cycles.
- **Per-sequence quality scores**: to identify overall poor reads.
- **Per-base GC content**: to check for unexpected biases.
- **Adapter content**: to screen for the presence of sequencing adapters.
- **Overrepresented sequences**: which may indicate contamination or technical artifacts.

These outputs were used to screen for any samples exhibiting poor sequence quality, severe adapter contamination, or other sequencing anomalies that could interfere with alignment or quantification. The resulting reports were archived and later integrated into a unified **MultiQC** report for batch-level inspection.

---

## Step 2: Adapter Removal and Quality Trimming

Following the initial quality check, sequencing adapters and low-quality bases were removed using **Trim Galore (v0.6.10)**, a wrapper around **Cutadapt** with optional FastQC execution.

For paired-end libraries, Trim Galore was applied to both forward and reverse reads in paired mode. For single-end libraries, the single FASTQ file was processed accordingly. The tool automatically detects common adapter sequences (e.g., Illumina) and removes them, along with trimming low-quality bases (Phred score < 20) from both ends of the reads.

Additionally, Trim Galore produces:
- Adapter-trimmed FASTQ files;
- Trimming summary reports detailing read retention and base quality before/after trimming;
- Post-trimming FastQC reports to evaluate improvements in data quality.

These cleaned and trimmed reads were used for downstream rRNA removal and alignment. The trimming step is essential to eliminate technical artifacts that may otherwise reduce mapping accuracy and introduce biases in quantification.

---

## Step 3: Removal of Ribosomal RNA Contaminants

To minimize the influence of highly abundant but non-informative ribosomal RNA (rRNA) sequences, we employed **SortMeRNA (v4.3.7)** for computational rRNA removal at the FASTQ level. This step helps to enrich for informative mRNA-derived reads and improves downstream alignment and quantification accuracy.

For paired-end libraries, both read files were input into SortMeRNA with the `--paired_in` option. For single-end libraries, only one read file was used. The process involved aligning reads against a curated rRNA reference database, including:
- **SILVA 18S and 28S rRNA** (eukaryotic large and small subunits),
- **Rfam 5S and 5.8S rRNA** databases.

Reads aligning to rRNA sequences were separated from the rest. The **non-rRNA reads** were retained for downstream STAR alignment, while log files were generated for quality tracking.

This step is particularly important in experiments where no rRNA depletion protocol was applied during library preparation, or where residual contamination remains. Removing these reads improves signal-to-noise ratio and ensures that the computational resources focus on informative sequences.

---

## Step 4: Splice-Aware Alignment with STAR

Cleaned, non-rRNA reads were aligned to the reference genome using **STAR (Spliced Transcripts Alignment to a Reference, v2.7.11b)**, a splice-aware aligner optimized for RNA-seq data. The genome index was generated dynamically within the workflow using the provided reference FASTA and GTF files, with `sjdbOverhang` automatically determined based on the read length.

For paired-end samples, both reads were input in `--readFilesIn` with `--paired_in`. For single-end samples, only one file was used. STAR was executed in **two-pass mode**, which improves splice junction detection and overall alignment accuracy.

Each alignment produced:
- A **sorted BAM file** (`*_Aligned.sortedByCoord.out.bam`) for genome-level visualization and QC;
- A **transcriptome-aligned BAM file** (`*_Aligned.toTranscriptome.out.bam`) used for downstream quantification by RSEM;
- A **log file** (`*_Log.final.out`) summarizing alignment statistics such as mapping rate, mismatch rate, and number of uniquely mapped reads.

STAR was configured with additional parameters to ensure optimal performance and specificity, including:
- Minimum intron length (`--alignIntronMin`),
- Mismatch tolerance (`--outFilterMismatchNoverReadLmax`),
- Read attribute tagging (`--outSAMattributes All`).

This step bridges the transition from raw sequencing data to structured, mappable transcriptomic information, setting the foundation for reliable expression quantification.

---

## Step 5: Library Strandedness Inference (RSeQC)

To ensure accurate expression quantification, the library strandedness was first inferred using **RSeQC (v5.0.4)** via the `infer_experiment.py` module. Each genome-aligned BAM file was evaluated alongside a filtered BED annotation derived from the GTF. The tool quantifies the proportion of reads mapping in the sense versus antisense orientation relative to gene annotation.

Based on empirical thresholds:
- If the strand bias difference was < 0.1, the library was classified as **unstranded**;
- If ≥ 80% of reads aligned in one orientation, the library was classified as **forward** or **reverse** accordingly.

The inferred strandedness was saved in a plain text file for each sample and directly used as a parameter in downstream RSEM quantification.

---

## Step 6: Gene and Isoform Quantification with RSEM

Transcriptome-aligned BAM files (`*_Aligned.toTranscriptome.out.bam`) were used to estimate gene- and isoform-level expression using **RSEM (v1.3.3)**.

RSEM was executed with the `--alignments` and `--strandedness` options, utilizing the reference generated earlier from the genome and GTF. It outputs normalized expression values as:
- **Expected counts**: raw transcript abundance (non-normalized);
- **TPM (Transcripts Per Million)**: normalized for sequencing depth and transcript length;
- **FPKM (Fragments Per Kilobase of transcript per Million mapped reads)**: traditional normalization metric.

Separate files were generated for **genes** and **isoforms**, including detailed statistics and metadata for quality control.

Together, **RSeQC** and **RSEM** provide a robust and biologically coherent framework for RNA-seq expression quantification, tailored to each sample’s orientation characteristics.

---

## Step 7: Expression Matrix Consolidation

To facilitate downstream statistical analyses and visualization, all RSEM quantification results were consolidated into unified expression matrices using a custom R script.

Specifically, the `.genes.results` and `.isoforms.results` files from all samples were parsed and merged into the following six tab-delimited files:
- `gene_expected_count_matrix.tsv`
- `gene_TPM_matrix.tsv`
- `gene_FPKM_matrix.tsv`
- `isoform_expected_count_matrix.tsv`
- `isoform_TPM_matrix.tsv`
- `isoform_FPKM_matrix.tsv`

Each matrix includes samples as columns and genes/transcripts as rows, allowing immediate use in tools such as **DESeq2**, **edgeR**, or custom visualization workflows. This step standardizes the data structure and supports reproducible downstream analyses.

---

## Step 8: Integrated Quality Control Reporting with MultiQC

All quality control outputs generated during the pipeline were aggregated using **MultiQC (v1.28)** into a single interactive HTML report.

MultiQC parses and summarizes results from:
- **FastQC** (raw and trimmed reads),
- **Trim Galore** (adapter removal summaries),
- **SortMeRNA** (rRNA removal logs),
- **STAR** (alignment metrics),
- **RSeQC** (strandedness, read distribution, BAM stats),
- **RSEM** (alignment and multimapping summaries).

This report enables comprehensive inspection of sequencing quality, sample consistency, and pipeline performance at a glance. It serves as both a documentation artifact and a valuable diagnostic tool for identifying outlier samples or processing issues.

---

## Reproducibility and Modularity
All processes are containerized (via Docker or Apptainer) to ensure version control and reproducibility. The pipeline is written in Nextflow DSL2 modular format, enabling scalable, maintainable, and extensible RNA-seq data processing.

---

