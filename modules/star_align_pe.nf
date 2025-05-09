#!/usr/bin/env nextflow

process STAR_ALIGN_PE {
    label 'process_high'
    
    container "containers/star_2.7.11b.sif"
    publishDir "results/star_align", mode: 'copy'

    input:
    tuple path(read1), path(read2)
    path index_zip

    output:
    path "*_Aligned.sortedByCoord.out.bam", emit: genome_bam
    path "*_Aligned.toTranscriptome.out.bam", emit: transcriptome_bam
    path "*_Log.final.out", emit: log

    script:
    def prefix = read1.simpleName
                .replaceFirst(/\.non_rRNA$/, '')
                .replaceFirst(/_1$/, '')

    """
    mkdir star_index
    tar -xzf $index_zip -C star_index

    STAR --runThreadN $task.cpus --runMode alignReads \
         --genomeDir star_index \
         --readFilesIn ${read1} ${read2} \
         --readFilesCommand zcat \
         --twopassMode Basic \
         --alignIntronMin 30 \
         --outFilterType BySJout \
         --outFilterScoreMinOverLread 0.1 \
         --outFilterMatchNminOverLread 0.1 \
         --outFilterMismatchNoverReadLmax 0.04 \
         --outSAMtype BAM SortedByCoordinate \
         --outFileNamePrefix ${prefix}_ \\
         --outSAMattributes All \
         --outSAMattrIHstart 0 \
         --outMultimapperOrder Random \
         --outReadsUnmapped Fastx \
         --quantMode TranscriptomeSAM GeneCounts
    """
}
