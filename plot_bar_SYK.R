library(dplyr)
source("run_enrichment_v3.2.R")

d = readRDS("test-dataset/DESeq2_table.RDS")

res_res = run_all(d, lfc = 'log2FoldChange', padj = 'padj', gene_name = 'gene_name', 
                   lfc_thresh = 0.25, padj_cutoff = 0.05,
                   organism = 'human')

plot_bar <- function(res){
  
  ### ORA
  res_ORA1 = res[["ORA_up"]]
  res_ORA1$direction = "Up"
  res_ORA1 = res_ORA1 %>% top_n(10, -p.adjust)
  
  res_ORA2 = res[["ORA_down"]]
  res_ORA2$direction = "Dn"
  res_ORA2 = res_ORA2 %>% top_n(10, -p.adjust)
  
  res_ORA = rbind.data.frame(res_ORA1, res_ORA2)
  
  res_ORA$padj_Dir = ifelse(res_ORA$direction == "Up", -log10(res_ORA$p.adjust), log10(res_ORA$p.adjust))
  res_ORA$binary_Dir = ifelse(res_ORA$direction == "Up", 1, 0)
  
  res_ORA$Description_ID_Dir = ifelse(res_ORA$direction == "Up", paste0(res_ORA$Description_ID, "   "),
                                       paste0("   ", res_ORA$Description_ID))
  
  ## Visualize
  p1 = ggplot(res_ORA, aes(x = padj_Dir, y = reorder(Description_ID_Dir, -log10(p.adjust)))) +
    geom_bar(stat = "identity", 
             color = 'black', lwd = 0.5,
             aes(fill = direction), alpha= 0.75,
             position = "identity", show.legend = F) +
    scale_fill_manual(values = c("Up" = "#E64B35", "Dn" = "#3182bd")) +
    geom_text(aes(y=Description_ID_Dir, x=0, label= Description_ID_Dir), hjust=res_ORA$binary_Dir) + 
    xlim(c(-max(-log10(res_ORA$p.adjust)), max(-log10(res_ORA$p.adjust)))) +
    #coord_fixed(ratio = 0.15) +
    labs(x = '-log10(FDR)', y = '') +
    #theme_minimal(base_size = 7) +
    geom_vline(xintercept = 0, colour = "black", size = 0.5) + # Add 0 lines
    theme_classic(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      axis.text.x = element_text(size = 10, colour="black", vjust = 1, hjust = 0.5),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(size = 12, vjust = 0.2),
      axis.line.y = element_blank(),
      plot.title = element_text(size = 12, face = "bold", vjust = 1.2, hjust = 0.5),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
  
  ### GSEA
  res_GSEA = res[["GSEA"]]
  res_GSEA = res_GSEA %>% top_n(20, abs(NES))
  
  res_GSEA$binary_Dir = ifelse(res_GSEA$NES > 0, 1, 0)
  
  res_GSEA$Description_ID_Dir = ifelse(res_GSEA$NES > 0, paste0(res_GSEA$Description_ID, "   "),
                                        paste0("   ", res_GSEA$Description_ID))
  
  ## Visualize
  p2 = ggplot(res_GSEA, aes(x = NES, y = reorder(Description_ID_Dir, NES))) +
    geom_bar(stat = "identity", 
             color = 'black', lwd = 0.5,
             aes(fill = NES), alpha= 1,
             position = "identity", show.legend = F) +
    scale_fill_gradient2(low = '#90C0DF', mid = 'white', high = '#C593C2', midpoint = 0) + 
    geom_text(aes(y=Description_ID_Dir, x=0, label= Description_ID_Dir), hjust=res_GSEA$binary_Dir) + 
    xlim(c(-max(res_GSEA$NES), max(res_GSEA$NES))) +
    #coord_fixed(ratio = 0.15) +
    labs(x = 'NES', y = '') +
    #theme_minimal(base_size = 7) +
    geom_vline(xintercept = 0, colour = "black", size = 0.5) + # Add 0 lines
    theme_classic(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      axis.text.x = element_text(size = 10, colour="black", vjust = 1, hjust = 0.5),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(size = 12, vjust = 0.2),
      axis.line.y = element_blank(),
      plot.title = element_text(size = 12, face = "bold", vjust = 1.2, hjust = 0.5),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
  
  p3 = plot_grid(p1, p2, labels = c("ORA", "GSEA"), label_size = 15)
  
  return(p3)
}

plot_bar(res = test_res)


