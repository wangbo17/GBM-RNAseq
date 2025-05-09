#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

genes_start <- which(args == "--genes") + 1
isoforms_start <- which(args == "--isoforms") + 1

genes_end <- isoforms_start - 2
isoforms_end <- length(args)

gene_files <- args[genes_start:genes_end]
isoform_files <- args[isoforms_start:isoforms_end]

library(readr)
library(dplyr)

read_matrix <- function(files, id_col) {
  ids <- NULL
  expected_count_list <- list()
  TPM_list <- list()
  FPKM_list <- list()

  for (file in files) {
    sample_name <- gsub("_Aligned\\.(genes|isoforms)\\.results", "", basename(file))
    df <- read_tsv(file, col_types = cols()) %>%
      select(all_of(id_col), expected_count, TPM, FPKM)

    if (is.null(ids)) ids <- df[[id_col]]
    expected_count_list[[sample_name]] <- df$expected_count
    TPM_list[[sample_name]] <- df$TPM
    FPKM_list[[sample_name]] <- df$FPKM
  }

  list(
    expected_count = as.data.frame(cbind(!!id_col := ids, do.call(cbind, expected_count_list))),
    TPM            = as.data.frame(cbind(!!id_col := ids, do.call(cbind, TPM_list))),
    FPKM           = as.data.frame(cbind(!!id_col := ids, do.call(cbind, FPKM_list)))
  )
}

write_matrix <- function(matrix_list, prefix) {
  write_tsv(matrix_list$expected_count, paste0(prefix, "_expected_count_matrix.tsv"))
  write_tsv(matrix_list$TPM,            paste0(prefix, "_TPM_matrix.tsv"))
  write_tsv(matrix_list$FPKM,           paste0(prefix, "_FPKM_matrix.tsv"))
}

write_matrix(read_matrix(gene_files, "gene_id"), "gene")
write_matrix(read_matrix(isoform_files, "transcript_id"), "isoform")
