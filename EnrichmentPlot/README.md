# EnrichmentPlot

```
library(dplyr)
library(ggplot2)
library(annotate)
source("run_enrichment_v7.R")

data = readRDS("test/DESeq2_KD.WT_RNA-Seq.RDS")

# Enrichment test
res = run_all(deg_sample, 
              organism = 'mouse', 
              lfc_column = "avg_log2FC", 
              padj_column = "p_val_adj", 
              gene_name_column = 'gene_name', 
              lfc_thresh = 0.25, 
              padj_cutoff = 0.05, 
              minGSSize = 100, 
              maxGSSize = 1000,
              output_file = 'sample_enrichment_result.xlsx')

# Visualization
res_pl = plot_all(res,
                  width = 20,
                  height = 17.5,
                  lfc_thresh = 0.25, 
                  padj_cutoff = 0.05, 
                  output_file = 'enrichment_result.pdf')
```
