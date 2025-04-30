#!/usr/bin/env nextflow

process STAR_ALIGN_SE {

    container "community.wave.seqera.io/library/star:2.7.11b--84fcc19fdfab53a4"
    publishDir "results/star_align", mode: 'copy'

    input:
    path reads
    path index_zip

    output:
    path "*_Aligned.sortedByCoord.out.bam", emit: genome_bam
    path "*_Aligned.toTranscriptome.out.bam", emit: transcriptome_bam
    path "*_Log.final.out", emit: log

    script:
    def prefix = reads.simpleName
                .replaceFirst(/raw$/, '')
                .replaceFirst(/\.non_rRNA/, '')
                
    """
    mkdir star_index
    tar -xzf $index_zip -C star_index

    STAR --runThreadN 6 --runMode alignReads \\
        --genomeDir star_index \\
        --readFilesIn $reads \\
        --readFilesCommand zcat \\
        --twopassMode Basic \\
        --alignIntronMin 30 \\
        --outFilterType BySJout \\
        --outFilterScoreMinOverLread 0.1 \\
        --outFilterMatchNminOverLread 0.1 \\
        --outFilterMismatchNoverReadLmax 0.04 \\
        --outSAMtype BAM SortedByCoordinate \\
        --outFileNamePrefix ${prefix}_ \\
        --outSAMattributes All \\
        --outSAMattrIHstart 0 \\
        --outMultimapperOrder Random \\
        --outReadsUnmapped Fastx \\
        --quantMode TranscriptomeSAM GeneCounts
    """
}

