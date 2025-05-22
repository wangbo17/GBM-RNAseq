#!/usr/bin/env Rscript

library(readr)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

genes_start <- which(args == "--genes") + 1
isoforms_start <- which(args == "--isoforms") + 1

genes_end <- isoforms_start - 2
isoforms_end <- length(args)

gene_files <- args[genes_start:genes_end]
isoform_files <- args[isoforms_start:isoforms_end]

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
    expected_count = {
      mat <- as.data.frame(cbind(ids, do.call(cbind, expected_count_list)))
      names(mat)[1] <- id_col
      mat
    },
    TPM = {
      mat <- as.data.frame(cbind(ids, do.call(cbind, TPM_list)))
      names(mat)[1] <- id_col
      mat
    },
    FPKM = {
      mat <- as.data.frame(cbind(ids, do.call(cbind, FPKM_list)))
      names(mat)[1] <- id_col
      mat
    }
  )
}

write_matrix <- function(matrix_list, prefix) {
  write_csv(matrix_list$expected_count, paste0(prefix, "_expected_count_matrix.csv"))
  write_csv(matrix_list$TPM, paste0(prefix, "_TPM_matrix.csv"))
  write_csv(matrix_list$FPKM, paste0(prefix, "_FPKM_matrix.csv"))
}

write_matrix(read_matrix(gene_files, "gene_id"), "gene")
write_matrix(read_matrix(isoform_files, "transcript_id"), "isoform")
