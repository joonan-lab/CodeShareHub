library(dplyr)
library(ggplot2)
source("~/Desktop/run_enrichment_v2.2.R")

data = readRDS("~/Downloads/DESeq2_table.RDS")

# runORA for up-DEGs
res = run_ORA(data, 
              df_source = "DESeq2", 
              organism = "hsapiens", 
              change = "pos")

res2 = run_GSEA(data,
                df_source = "DESeq2",
                bg_file = "~/Downloads/h.all.v2022.1.Hs.symbols.gmt")

plot_bar(res = res, 
         isORA = T, n_term = 20, 
         width = 7.5, height = 7.5, filename = "tmp.pdf")

plot_bar(res = res2, 
        isORA = F, n_term = 20, 
        width = 7.5, height = 7.5, filename = "tmp2.pdf")

