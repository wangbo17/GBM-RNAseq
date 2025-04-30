#!/usr/bin/env nextflow

process RSEM_REFERENCE {

    container "community.wave.seqera.io/library/rsem_star:8a0fe0c5e7aa01d5"
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
      --gtf $gtf \\
      --star \\
      $fasta \\
      rsem_ref/rsem_ref
    """
}
