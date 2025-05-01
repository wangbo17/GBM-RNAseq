#!/usr/bin/env nextflow

process RSEM_QUANT {

    container "community.wave.seqera.io/library/rsem_star:8a0fe0c5e7aa01d5"
    publishDir "results/rsem_expression", mode: 'copy'

    input:
    path transcriptome_bam
    path rsem_index
    val strandedness

    output:
    path "*.genes.results", emit: gene_result
    path "*.isoforms.results", emit: isoform_result
    path "*.cnt", emit: cnt_report

    script:
    def prefix = transcriptome_bam.simpleName
                .replaceFirst(/_Aligned\.toTranscriptome\.out$/, '')
                
    """
    rsem-calculate-expression -p 8 \\
      --alignments \\
      --strandedness ${strandedness} \\
      ${transcriptome_bam} \\
      ${rsem_index}/rsem_ref \\
      ${prefix}
    
    cp ${prefix}.stat/${prefix}.cnt ${prefix}.cnt
    """
}
