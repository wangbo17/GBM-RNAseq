#!/usr/bin/env nextflow

process STAR_INDEX {
    label 'process_high'
    
    container "containers/star_2.7.11b.sif"
    publishDir "results/star_index", mode: 'copy'

    input:
    path reads
    path fasta
    path gtf

    output:
    path "star_index.tar.gz", emit: index_zip

    script:
    def memory      = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
    """
    READ_LENGTH=\$(zcat $reads | awk 'NR % 4 == 2 {print length(\$0)}' | head -n 100000 | sort -nr | head -n 1)
    SJDB_OVERHANG=\$((READ_LENGTH - 1))

    mkdir -p star_index

    STAR --runThreadN $task.cpus $memory --runMode genomeGenerate \
         --genomeDir star_index \
         --genomeFastaFiles $fasta \
         --sjdbGTFfile $gtf \
         --sjdbOverhang \$SJDB_OVERHANG

    tar -czvf star_index.tar.gz -C star_index .
    """
}
