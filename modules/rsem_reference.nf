#!/usr/bin/env nextflow

process RSEM_REFERENCE {
    label 'process_high'
    
    container 'oras://community.wave.seqera.io/library/rsem_star:8a0fe0c5e7aa01d5'
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
