#!/bin/bash
#$ -V -cwd
#$ -l node_type=40core-192G
#$ -l h_rt=48:00:00
#$ -l h_vmem=8G
#$ -pe smp 1
#$ -m be
#$ -j y
#$ -o rnaseq.out

# ==============================
# Run Pipeline
# ==============================
nextflow run rnaseq_pe.nf --input_csv data/paired-end.csv -resume
