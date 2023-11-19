# EnrichmentPlot

<p align="center">
    <img src="https://github.com/joonan-lab/CodeShareHub/assets/47490862/62fcf4c1-3fd0-4f99-84bc-1e17e921ca45" alt="image" width="800"/>
</p>


```
library(dplyr)
library(ggplot2)
source("run_enrichment_v4.R")

data = readRDS("test-dataset/DESeq2_table.RDS")
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
              output_file = 'test-result/enrichment_result.xlsx')

res_pl = plot_all(res,
                  width = 22.5,
                  height = 20,
                  output_file = 'test-result/enrichment_result.pdf')
```
