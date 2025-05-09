#!/usr/bin/env nextflow

process FASTQC_PE {
    label 'process_single'

    container 'containers/fastqc_0.12.1.sif'
    publishDir "results/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*_fastqc.zip", emit: zip
    path "*_fastqc.html", emit: html

    script:
    def memory_in_mb = MemoryUnit.of("${task.memory}").toUnit('MB') / task.cpus
    def fastqc_memory = memory_in_mb > 10000 ? 10000 : (memory_in_mb < 100 ? 100 : memory_in_mb)

    """
    read1_name=\$(basename ${reads[0]})
    read2_name=\$(basename ${reads[1]})

    target1=${sample_id}_raw_1.fastq.gz
    target2=${sample_id}_raw_2.fastq.gz

    if [ "\$read1_name" != "\$target1" ]; then
        ln -s \$read1_name \$target1
    fi

    if [ "\$read2_name" != "\$target2" ]; then
        ln -s \$read2_name \$target2
    fi

    fastqc --threads $task.cpus --memory $fastqc_memory \$target1 \$target2
    """
}