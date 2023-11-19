library(dplyr)
library(ggplot2)
source("run_enrichment_v3.2.R")

data = readRDS("test-dataset/DESeq2_table.RDS")

colnames(data)
# runORA for up-DEGs
res = run_all(data, lfc_column = "log2FoldChange", 
              padj_column = "padj", 
              gene_name_column = 'gene_name', 
              lfc_thresh = 0.25, 
              padj_cutoff = 0.05, 
              organism = 'human', 
              minGSSize = 100, 
              maxGSSize = 2000,
              output_file = 'enrichment_result.xlsx')


plot_volcano(res$DEGs)

plot_treeplot(res)
