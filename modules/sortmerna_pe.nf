#!/usr/bin/env nextflow

process SORTMERNA_PE {

    container './containers/sortmerna_4.3.7--5e99cfae38f7f1eb.sif'
    publishDir "results/sortmerna", mode: 'copy'

    input:
    tuple path(read1), path(read2), val(fasta_refs)

    output:
    tuple path("*_1.non_rRNA.fastq.gz"), path("*_2.non_rRNA.fastq.gz"), emit: cleaned_reads
    path "*.sortmerna.log", emit: logs

    script:
    def prefix = read1.simpleName
                .replaceFirst(/_val_1$/, '')
    def ref_string = fasta_refs.collect { "--ref ${it}" }.join(" \\\n    ")

    """
    sortmerna \\
        ${ref_string} \\
        --reads ${read1} \\
        --reads ${read2} \\
        --paired_in \\
        --aligned rRNA_reads --fastx --other non_rRNA_reads --out2 \\
        --workdir . \\
        --threads ${task.cpus}

    mv non_rRNA_reads_fwd.f*q.gz ${prefix}_1.non_rRNA.fastq.gz
    mv non_rRNA_reads_rev.f*q.gz ${prefix}_2.non_rRNA.fastq.gz
    mv rRNA_reads.log ${prefix}.sortmerna.log

    sed -i "s/Reads file: .*_val_[12]\\.fq\\.gz/Reads file: ${prefix}/g" ${prefix}.sortmerna.log
    """
}