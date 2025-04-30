#!/usr/bin/env nextflow

process SORTMERNA_SE {

    container './containers/sortmerna_4.3.7--5e99cfae38f7f1eb.sif'
    publishDir "results/sortmerna", mode: 'copy'

    input:
    tuple path(read), val(fasta_refs)

    output:
    path "*.non_rRNA.fastq.gz", emit: cleaned_reads
    path "*.sortmerna.log", emit: logs

    script:
    def prefix = read.simpleName
                .replaceFirst(/_val_trimmed$/, '')
    def ref_string = fasta_refs.collect { "--ref ${it}" }.join(" \\\n    ")

    """
    sortmerna \\
        ${ref_string} \\
        --reads ${read} \\
        --aligned rRNA_reads --fastx --other non_rRNA_reads \\
        --workdir . \\
        --threads ${task.cpus}

    mv non_rRNA_reads.fq.gz ${prefix}.non_rRNA.fastq.gz
    mv rRNA_reads.log ${prefix}.sortmerna.log

    sed -i "s/Reads file: .*_val_[12]\\.fq\\.gz/Reads file: ${prefix}/g" ${prefix}.sortmerna.log
    """
}
