rm(list = ls())
library(dplyr)
library(ggplot2)
library(annotate)
setwd("~/Dropbox/CodeShareHub/EnrichmentPlot/")
source("run_enrichment_v4.2_GOBP.R")

data = readRDS("test/test-dataset/DESeq2_table.RDS")
colnames(data)

# Enrichment test
res = run_all(data, lfc_column = "log2FoldChange", 
              padj_column = "padj", 
              gene_name_column = 'gene_name', 
              lfc_thresh = 0.25, 
              padj_cutoff = 0.05, 
              organism = 'human', 
              minGSSize = 100, 
              maxGSSize = 2000,
              output_file = 'test/test-result/enrichment_result_bp.xlsx')

res_pl = plot_all(res,
                  width = 22.5,
                  height = 20,
                  output_file = 'test/test-result/enrichment_result_bp.pdf')
