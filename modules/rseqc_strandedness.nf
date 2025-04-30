#!/usr/bin/env nextflow

process RSEQC_STRANDEDNESS {

    container "community.wave.seqera.io/library/rseqc:5.0.4--8c855f0b915334d2" 
    publishDir "results/infer_strandedness", mode: 'copy'

    input:
    path bam_file
    path bed_file

    output:
    path "${bam_file.simpleName}.infer_experiment.txt", emit: infer_raw
    path "${bam_file.simpleName}.strandedness.txt", emit: strandedness_report

    script:
    """
    infer_experiment.py -r $bed_file -i $bam_file > ${bam_file.simpleName}.infer_experiment.txt

    FORWARD=\$(grep 'Fraction of reads explained by "++,--"' ${bam_file.simpleName}.infer_experiment.txt | awk '{print \$NF}')
    REVERSE=\$(grep 'Fraction of reads explained by "+-,-+"' ${bam_file.simpleName}.infer_experiment.txt | awk '{print \$NF}')


    FORWARD_NUM=\$(echo "\$FORWARD" | awk '{printf "%f", \$1}')
    REVERSE_NUM=\$(echo "\$REVERSE" | awk '{printf "%f", \$1}')

    DIFF=\$(echo "\$FORWARD_NUM - \$REVERSE_NUM" | bc -l)
    DIFF_ABS=\$(echo "\$DIFF" | awk '{print (\$1<0)?-\$1:\$1}')

    if (( \$(echo "\$DIFF_ABS < 0.1" | bc -l) )); then
        echo "none" > ${bam_file.simpleName}.strandedness.txt
    elif (( \$(echo "\$FORWARD_NUM >= 0.8" | bc -l) )); then
        echo "forward" > ${bam_file.simpleName}.strandedness.txt
    elif (( \$(echo "\$REVERSE_NUM >= 0.8" | bc -l) )); then
        echo "reverse" > ${bam_file.simpleName}.strandedness.txt
    else
        echo "none" > ${bam_file.simpleName}.strandedness.txt
    fi
    """
}
