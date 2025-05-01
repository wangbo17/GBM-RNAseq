#!/usr/bin/env nextflow

// Module INCLUDE statements
include { PREPARE_GTF } from './modules/prepare_gtf.nf'
include { FASTQC_PE } from './modules/fastqc_pe.nf'
include { TRIM_GALORE_PE } from './modules/trim_galore_pe.nf'
include { SORTMERNA_PE } from './modules/sortmerna_pe.nf'
include { STAR_INDEX } from './modules/star_index.nf'
include { STAR_ALIGN_PE } from './modules/star_align_pe.nf'
include { RSEQC_STRANDEDNESS } from './modules/rseqc_strandedness.nf'
include { RSEQC_BAMSTAT } from './modules/rseqc_bamstat.nf'
include { RSEQC_DISTRIBUTION } from './modules/rseqc_distribution.nf'
include { RSEM_REFERENCE } from './modules/rsem_reference.nf'
include { RSEM_QUANT } from './modules/rsem_quant.nf'
include { RSEM_MERGE } from './modules/rsem_merge.nf'
include { MULTIQC } from './modules/multiqc.nf'

/*
 * Pipeline parameters
 */

// Primary input
params.input_csv = "data/paired-end.csv"
params.fasta = "data/reference/genome.fa"
params.gtf = "data/reference/genes.gtf"
params.report_id = "all_paired-end"

log.info """\
    ==================================================================================================
          ____ ____  __  __   ____        _ _      ____  _   _    _        ____             
         / ___| __ )|  \\/  | | __ ) _   _| | | __ |  _ \\| \\ | |  / \\      / ___|  ___  __ _ 
        | |  _|  _ \\| |\\/| | |  _ \\| | | | | |/ / | |_) |  \\| | / _ \\ ____\\___ \\ / _ \\/ _` |
        | |_| | |_) | |  | | | |_) | |_| | |   <  |  _ <| |\\  |/ ___ |________) |  __| (_| |
         \\____|____/|_|  |_| |____/ \\__,_|_|_|\\_\\ |_| \\_|_| \\_/_/   \\_\\   |____/ \\___|\\__, |
                                                                                         |_|
    ==================================================================================================

      ▸ Workflow Summary:

        ▹ Quality Control
           - FASTQC: Evaluate raw read quality.
           - Trim Galore: Remove adapters and low-quality bases.

        ▹ rRNA Removal
           - SortMeRNA: Remove rRNA reads at the FASTQ level to reduce non-informative signal.

        ▹ Reference Preparation
           - GTF Filtering: Filter GTF annotation for relevant biotypes and generate BED for RSeQC.
           - STAR Index: Generate STAR genome index tailored to read length.
           - RSEM Reference: Prepare RSEM reference for quantification.

        ▹ Alignment & Quantification
           - STAR: Align reads to genome and transcriptome using STAR.
           - RSEM: Quantify gene and isoform expression using RSEM.

        ▹ Quality Assessment
           - RSeQC: Assess strandedness, alignment quality, and read distribution.

        ▹ Results Summary
           - RSEM Merge: Combine gene and isoform expression estimates across all samples.
           - MultiQC: Summarize QC metrics into a single HTML report.

    --------------------------------------------------------------------------------------------------

      ▸ Key Parameters:

        ▹ Executed By       : ${System.getProperty('user.name')}
        ▹ Run Date          : ${workflow.start.format("yyyy-MM-dd HH:mm 'UTC'")}

        ▹ Config Profile    : ${workflow.profile}
        ▹ Read Type         : ${params.report_id}
        ▹ Sample Sheet      : ${params.input_csv}
        ▹ Reference FASTA   : ${params.fasta}
        ▹ GTF Annotation    : ${params.gtf}
                                                                     Author: Bo Wang | Version: Beta
    ==================================================================================================
    """
    .stripIndent(true)


workflow {

    // Create input channel
    read_ch = Channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> tuple(row.sample_id, [file(row.fastq_1), file(row.fastq_2)]) }

    // Call processes
    PREPARE_GTF(file(params.gtf))

    FASTQC_PE(read_ch)

    TRIM_GALORE_PE(read_ch)

    SORTMERNA_PE(
        TRIM_GALORE_PE.out.trimmed_reads.map { read1, read2 ->
            tuple(
            read1,
            read2,
            [
                file("data/rrna_db/rfam-5.8s-database-id98.fasta"),
                file("data/rrna_db/rfam-5s-database-id98.fasta"),
                file("data/rrna_db/silva-euk-18s-id95.fasta"),
                file("data/rrna_db/silva-euk-28s-id98.fasta")
            ]
            )
        }
    )

    STAR_INDEX(
        read_ch.map { it[1] }.first(),
        file(params.fasta),
        PREPARE_GTF.out.filtered_gtf
    )

    STAR_ALIGN_PE(SORTMERNA_PE.out.cleaned_reads, STAR_INDEX.out.index_zip)

    RSEQC_STRANDEDNESS(STAR_ALIGN_PE.out.genome_bam, PREPARE_GTF.out.filtered_bed)
    RSEQC_BAMSTAT(STAR_ALIGN_PE.out.genome_bam)
    RSEQC_DISTRIBUTION(STAR_ALIGN_PE.out.genome_bam, PREPARE_GTF.out.filtered_bed)

    RSEM_REFERENCE(file(params.fasta), PREPARE_GTF.out.filtered_gtf)

    RSEM_QUANT(
        STAR_ALIGN_PE.out.transcriptome_bam,
        RSEM_REFERENCE.out.rsem_rindex,
        RSEQC_STRANDEDNESS.out.strandedness_report.map { file -> file.text.trim() }
    )

    RSEM_MERGE(
        RSEM_QUANT.out.gene_result.collect(),
        RSEM_QUANT.out.isoform_result.collect()
    )

    MULTIQC(
        FASTQC_PE.out.zip.mix(
            FASTQC_PE.out.html,
            TRIM_GALORE_PE.out.trimming_reports,
            TRIM_GALORE_PE.out.fastqc_reports_1,
            TRIM_GALORE_PE.out.fastqc_reports_2,
            SORTMERNA_PE.out.logs,
            STAR_ALIGN_PE.out.log,
            RSEQC_STRANDEDNESS.out.infer_raw,
            RSEQC_BAMSTAT.out.bamstat_report,
            RSEQC_DISTRIBUTION.out.read_distribution,
            RSEM_QUANT.out.cnt_report
        ).collect(),
        params.report_id
    )
}
