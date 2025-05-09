#!/usr/bin/env nextflow

process RSEQC_STRANDEDNESS_SE {
    label 'process_low'
    
    container 'containers/rseqc_5.0.4.sif'
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

    awk -v f=\$FORWARD -v r=\$REVERSE -v out=${bam_file.simpleName}.strandedness.txt '
    BEGIN {
        diff = f - r
        diff_abs = (diff < 0) ? -diff : diff

        if (diff_abs < 0.1)
            result = "none"
        else if (f >= 0.8)
            result = "forward"
        else if (r >= 0.8)
            result = "reverse"
        else
            result = "none"

        print result > out
    }
    '
    """
}