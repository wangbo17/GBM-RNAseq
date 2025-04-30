#!/usr/bin/env nextflow

process PREPARE_GTF {

    input:
    path gtf_file

    output:
    path "genes.filtered.gtf", emit: filtered_gtf
    path "genes.filtered.bed", emit: filtered_bed

    script:
    """
    # Remove the following 15 gene biotypes that are not relevant or may interfere with downstream analysis:
    # - "rRNA", "rRNA_pseudogene", "Mt_rRNA": Ribosomal RNAs; highly expressed but non-informative, often dominate expression profiles.
    # - "snRNA", "snoRNA", "scaRNA", "scRNA", "vault_RNA": Small nuclear/nucleolar RNAs; non-coding and not the focus of this study.
    # - "miRNA", "misc_RNA", "sRNA", "ncRNA": Small non-coding RNAs; typically low in abundance and difficult to resolve in bulk RNA-seq data.
    # - "TEC": Transcripts of uncertain function ("To be Experimentally Confirmed"); often unstable or non-functional.
    # - "tRNA", "Mt_tRNA": Transfer RNAs; not involved in mRNA transcriptional regulation and typically excluded from differential expression analysis.
    awk '
    BEGIN {
      FS = OFS = "\t";
      keep["protein_coding"];
      keep["lncRNA"];
      keep["pseudogene"];
      keep["processed_pseudogene"];
      keep["transcribed_processed_pseudogene"];
      keep["transcribed_unitary_pseudogene"];
      keep["transcribed_unprocessed_pseudogene"];
      keep["translated_processed_pseudogene"];
      keep["unitary_pseudogene"];
      keep["unprocessed_pseudogene"];
      keep["IG_C_gene"]; keep["IG_C_pseudogene"]; keep["IG_D_gene"];
      keep["IG_J_gene"]; keep["IG_J_pseudogene"]; keep["IG_pseudogene"];
      keep["IG_V_gene"]; keep["IG_V_pseudogene"];
      keep["TR_C_gene"]; keep["TR_D_gene"]; keep["TR_J_gene"];
      keep["TR_J_pseudogene"]; keep["TR_V_gene"]; keep["TR_V_pseudogene"];
    }
    {
      if (\$0 ~ /^#/) { print; next }
      match(\$0, /gene_biotype "([^"]+)"/, a)
      if (a[1] in keep) print
    }' ${gtf_file} > genes.filtered.gtf

    awk '\$3 == "exon" {
        print \$1"\\t"(\$4-1)"\\t"\$5"\\texon\\t0\\t"\$7"\\t"\$4"\\t"\$5"\\t0\\t1\\t"(\$5-\$4+1)"\\t0"
    }' genes.filtered.gtf > genes.filtered.bed
    """
}
