#!/usr/bin/env nextflow

process RSEM_QUANT_PE {
    label 'process_medium'
    
    container 'containers/rsem_1.3.3_star_2.7.11b.sif'
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
    rsem-calculate-expression --num-threads $task.cpus \\
      --alignments \\
      --paired-end \\
      --strandedness ${strandedness} \\
      ${transcriptome_bam} \\
      ${rsem_index}/rsem_ref \\
      ${prefix}
    
    cp ${prefix}.stat/${prefix}.cnt ${prefix}.cnt
    """
}
