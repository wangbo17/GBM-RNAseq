#!/usr/bin/env nextflow

process RSEM_REFERENCE {
    label 'process_high'
    
    container 'containers/rsem_1.3.3_star_2.7.11b.sif'
    publishDir "results/rsem_reference", mode: 'copy'

    input:
    path fasta
    path gtf

    output:
    path "rsem_ref", emit: rsem_rindex

    script:
    """
    mkdir rsem_ref
    rsem-prepare-reference \\
      --num-threads $task.cpus \\
      --gtf $gtf \\
      --star \\
      $fasta \\
      rsem_ref/rsem_ref
    """
}
