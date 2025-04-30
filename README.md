# GBM RNA-Seq Processing Pipeline

<pre>
==============================================================================================
  
     ____ ____  __  __   ____        _ _      ____  _   _    _        ____             
    / ___| __ )|  \/  | | __ ) _   _| | | __ |  _ \| \ | |  / \      / ___|  ___  __ _ 
   | |  _|  _ \| |\/| | |  _ \| | | | | |/ / | |_) |  \| | / _ \ ____\___ \ / _ \/ _` |
   | |_| | |_) | |  | | | |_) | |_| | |   <  |  _ <| |\  |/ ___ |________) |  __| (_| |
    \____|____/|_|  |_| |____/ \__,_|_|_|\_\ |_| \_|_| \_/_/   \_\   |____/ \___|\__, |
                                                                                    |_|
      
                                                              Author: Bo Wang | Version: Beta
===============================================================================================
</pre>

## Overview
This pipeline performs comprehensive preprocessing, alignment, quantification, and quality control of RNA-seq data using modular Nextflow processes. It supports both paired-end and single-end reads and is containerized for reproducibility. The workflow integrates industry-standard tools including **FastQC**, **Trim Galore**, **SortMeRNA**, **STAR**, **RSEM**, and **RSeQC**, with summary reporting via **MultiQC**.

---

## Steps and Rationale

### 1. Initial Quality Control (FastQC)
- **Tool**: FastQC  
- **Purpose**: Assess raw sequencing quality (e.g., per-base quality, adapter contamination, GC content).  
- **Rationale**: Ensures that input data meet minimum quality standards before downstream processing.

### 2. Adapter Trimming (Trim Galore)
- **Tool**: Trim Galore (Cutadapt-based)  
- **Purpose**: Remove adapter sequences and low-quality bases.  
- **Rationale**: Improves read mappability and reduces bias in alignment and quantification.

### 3. rRNA Read Removal (SortMeRNA)
- **Tool**: SortMeRNA  
- **Purpose**: Filter out reads mapping to ribosomal RNA (rRNA) reference sequences.  
- **Rationale**: rRNA is highly abundant and non-informative in mRNA-focused studies; removal enhances transcriptome coverage.

### 4. GTF Filtering and BED Generation
- **Tool**: Custom AWK-based script  
- **Purpose**: Remove non-relevant gene biotypes and convert GTF to BED for RSeQC.  
- **Rationale**: Retains biologically meaningful annotations (e.g., `protein_coding`, `lncRNA`); improves downstream QC accuracy.

### 5. Genome Index Generation (STAR & RSEM)
- **Tool**: STAR, RSEM  
- **Purpose**: Build genome indices for read alignment and quantification.  
- **Rationale**: STAR index is tailored to read length (SJDB overhang); RSEM index ensures compatibility for expression estimation.

### 6. Read Alignment (STAR)
- **Tool**: STAR  
- **Purpose**: Align reads to the reference genome and transcriptome.  
- **Rationale**: Provides both genome-level BAM for QC and transcriptome-level BAM for expression quantification. Supports accurate spliced alignment.

### 7. Strandedness Inference (RSeQC)
- **Tool**: `infer_experiment.py` (RSeQC)  
- **Purpose**: Determine RNA-seq library strandedness.  
- **Rationale**: Accurate strandedness improves quantification and annotation-aware processing.

### 8. Alignment Quality Metrics (RSeQC)
- **Tools**: `bam_stat.py`, `read_distribution.py` (RSeQC)  
- **Purpose**: Evaluate alignment statistics and genomic feature distribution.  
- **Rationale**: Enables identification of potential biases or experimental artifacts.

### 9. Gene and Isoform Quantification (RSEM)
- **Tool**: RSEM  
- **Purpose**: Estimate transcript- and gene-level expression (expected counts, TPM, FPKM).  
- **Rationale**: EM-based model accounts for multi-mapping and isoform uncertainty; supports differential expression analysis.

### 10. Expression Matrix Generation
- **Tool**: Custom R script  
- **Purpose**: Merge individual RSEM results into consolidated matrices for downstream analysis.  
- **Rationale**: Enables cohort-level expression profiling and statistical testing.

### 11. Comprehensive Report (MultiQC)
- **Tool**: MultiQC  
- **Purpose**: Aggregate quality control metrics from all tools into a single report.  
- **Rationale**: Facilitates rapid review of data quality and processing success.

---

## Output Summary
- **Quality Reports**: FastQC HTMLs, trimming and rRNA removal logs, STAR alignment logs  
- **Alignment Files**: Sorted genome BAMs, transcriptome BAMs  
- **Quantification Files**: `*.genes.results`, `*.isoforms.results`, merged expression matrices  
- **QC Metrics**: RSeQC strandedness, read distribution, alignment statistics  
- **Final Summary**: MultiQC HTML report

---

## Reproducibility and Modularity
All processes are containerized (via Docker or Apptainer) to ensure version control and reproducibility. The pipeline is written in Nextflow DSL2 modular format, enabling scalable, maintainable, and extensible RNA-seq data processing.

---

