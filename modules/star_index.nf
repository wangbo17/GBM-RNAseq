#!/usr/bin/env nextflow

process STAR_INDEX {

    container "community.wave.seqera.io/library/star:2.7.11b--84fcc19fdfab53a4"
    publishDir "results/star_index", mode: 'copy'

    input:
    path reads
    path fasta
    path gtf

    output:
    path "star_index.tar.gz", emit: index_zip

    script:
    """
    READ_LENGTH=\$(zcat $reads | awk 'NR % 4 == 2 {print length(\$0)}' | head -n 100000 | sort -nr | head -n 1)
    SJDB_OVERHANG=\$((READ_LENGTH - 1))

    mkdir -p star_index

    STAR --runThreadN 6 --runMode genomeGenerate \
         --genomeDir star_index \
         --genomeFastaFiles $fasta \
         --sjdbGTFfile $gtf \
         --sjdbOverhang \$SJDB_OVERHANG

    tar -czvf star_index.tar.gz -C star_index .
    """
}
