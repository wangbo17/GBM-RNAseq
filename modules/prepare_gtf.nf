#!/usr/bin/env nextflow

process PREPARE_GTF {
    label 'process_low'

    input:
    path gtf_file

    output:
    path "excluded.gtf", emit: excluded_gtf
    path "filtered.gtf", emit: filtered_gtf
    path "filtered.bed", emit: filtered_bed
    path "gtf_report.txt", emit: gtf_report

    script:
    """
    # GTF Filtering Strategy for Bulk RNA-seq Analysis
    #
    # Purpose:
    # Remove low-confidence or undesired transcript biotypes to improve the accuracy of alignment,
    # quantification, and downstream interpretation.
    #
    # Biotypes to exclude:
    # rRNA, rRNA_pseudogene, Mt_rRNA, Mt_tRNA, miRNA, sRNA, snoRNA, scaRNA, scRNA, snRNA,
    # misc_RNA, vault_RNA, TEC, non_stop_decay, nonsense_mediated_decay, retained_intron,
    # artifact, pseudogene, ribozyme
    #
    # Filtering logic:
    # 1. Parse all transcript lines from the input GTF.
    # 2. Retain transcript_ids whose transcript_biotype is in the curated 'keep' list,
    #    which includes: protein_coding, lncRNA, processed_transcript, and selected pseudogenes.
    # 3. Also record their corresponding gene_ids.
    # 4. Retain a GTF line if:
    #    - It contains a transcript_id present in the keep list, OR
    #    - It lacks a transcript_id but has a gene_id present in the keep list.
    #    → This ensures all relevant annotations (gene, exon, CDS, UTR, etc.) are preserved.
    # 5. Lines not meeting the criteria are written to excluded.gtf for auditing.
    #
    # Outputs:
    # - filtered.gtf:     Final cleaned GTF used for STAR and RSEM.
    # - filtered.bed:     BED file derived from filtered.gtf, used for RSeQC.
    # - excluded.gtf:     GTF entries filtered out due to undesired biotypes.
    # - gtf_report.txt:   Summary of feature types and biotype distributions (before and after filtering).

    awk '
    BEGIN {
      FS = OFS = "\t";
      keep["protein_coding"];
      keep["lncRNA"];
      keep["processed_transcript"];
      keep["protein_coding_CDS_not_defined"];
      keep["protein_coding_LoF"];
      keep["translated_processed_pseudogene"];

      keep["processed_pseudogene"];
      keep["transcribed_unitary_pseudogene"];
      keep["transcribed_processed_pseudogene"];
      keep["transcribed_unprocessed_pseudogene"];
      keep["unprocessed_pseudogene"];
      keep["unitary_pseudogene"];

      keep["IG_V_gene"]; keep["IG_D_gene"]; keep["IG_J_gene"]; keep["IG_C_gene"];
      keep["IG_V_pseudogene"]; keep["IG_J_pseudogene"]; keep["IG_C_pseudogene"];
      keep["IG_pseudogene"];

      keep["TR_V_gene"]; keep["TR_D_gene"]; keep["TR_J_gene"]; keep["TR_C_gene"];
      keep["TR_V_pseudogene"]; keep["TR_J_pseudogene"];
    }
    \$3 == "transcript" {
      match(\$0, /transcript_id \"([^\"]+)\"/, t);
      match(\$0, /transcript_biotype \"([^\"]+)\"/, b);
      match(\$0, /gene_id \"([^\"]+)\"/, g);
      if (t[1] != "" && b[1] in keep) {
        print t[1] > "keep.transcript_ids.txt";
        print g[1] > "keep.gene_ids.txt";
      }
    }
    ' ${gtf_file}

    sort -u keep.transcript_ids.txt > tmp.tids && mv tmp.tids keep.transcript_ids.txt
    sort -u keep.gene_ids.txt > tmp.gids && mv tmp.gids keep.gene_ids.txt

    awk '
    BEGIN { FS = OFS = "\t" }
    ARGIND == 1 { keep_tid[\$1]; next }
    ARGIND == 2 { keep_gid[\$1]; next }
    \$0 ~ /^#/ { print; next }
    {
      match(\$0, /transcript_id "([^"]+)"/, t)
      match(\$0, /gene_id "([^"]+)"/, g)

      tid = (t[1] != "") ? t[1] : ""
      gid = (g[1] != "") ? g[1] : ""

      if ((tid != "" && tid in keep_tid) || (tid == "" && gid in keep_gid))
        print > "filtered.gtf"
      else
        print > "excluded.gtf"
    }
    ' keep.transcript_ids.txt keep.gene_ids.txt ${gtf_file}

    awk '
    \$3 == "exon" {
        print \$1"\t"(\$4-1)"\t"\$5"\texon\t0\t"\$7"\t"\$4"\t"\$5"\t0\t1\t"(\$5-\$4+1)"\t0"
    }
    \$3 == "CDS" {
        print \$1"\t"(\$4-1)"\t"\$5"\tCDS\t0\t"\$7"\t"\$4"\t"\$5"\t0\t1\t"(\$5-\$4+1)"\t0"
    }
    \$3 == "five_prime_utr" {
        print \$1"\t"(\$4-1)"\t"\$5"\tfive_prime_utr\t0\t"\$7"\t"\$4"\t"\$5"\t0\t1\t"(\$5-\$4+1)"\t0"
    }
    \$3 == "three_prime_utr" {
        print \$1"\t"(\$4-1)"\t"\$5"\tthree_prime_utr\t0\t"\$7"\t"\$4"\t"\$5"\t0\t1\t"(\$5-\$4+1)"\t0"
    }
    \$3 == "start_codon" {
        print \$1"\t"(\$4-1)"\t"\$5"\tstart_codon\t0\t"\$7"\t"\$4"\t"\$5"\t0\t1\t"(\$5-\$4+1)"\t0"
    }
    \$3 == "stop_codon" {
        print \$1"\t"(\$4-1)"\t"\$5"\tstop_codon\t0\t"\$7"\t"\$4"\t"\$5"\t0\t1\t"(\$5-\$4+1)"\t0"
    }
    ' filtered.gtf > filtered.bed

    {
      echo "===== FILTERED GTF Feature Type Distribution ====="
      cut -f3 filtered.gtf | sort | uniq -c | sort -nr
      echo ""
      echo "===== FILTERED transcript_biotype Distribution ====="
      grep -w 'transcript' filtered.gtf | grep -o 'transcript_biotype "[^"]*"' | sort | uniq -c | sort -nr
      echo ""
      echo "===== EXCLUDED GTF Feature Type Distribution ====="
      cut -f3 excluded.gtf | sort | uniq -c | sort -nr
      echo ""
      echo "===== EXCLUDED transcript_biotype Distribution ====="
      grep -w 'transcript' excluded.gtf | grep -o 'transcript_biotype "[^"]*"' | sort | uniq -c | sort -nr
    } > gtf_report.txt
    """
}
