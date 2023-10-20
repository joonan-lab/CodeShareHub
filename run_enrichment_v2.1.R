library(tidyverse)
library(gprofiler2)
library(fgsea)
library(writexl)

run_ORA <- function(df, df_source = NA, organism = NA, change = NA, lfc_thresh = 0.25, p.adj_thresh = 0.05, outtag = NA){
  if(is.na(df_source)){
    stop('argument "df_source" is missing, supports Seurat or DESeq')
  }
  if(is.na(change)){
    stop('argument "change" is missing, must be pos or neg')
  }
  if(is.na(organism)){
    stop('argument "organism" is missing, supports hsapiens or mmusculus')
  }
  
  
  if(change == 'pos'){
    cat(paste('Extracting Up-regulated genesets with log2FC ', lfc_thresh, ' & p-value threshold ', p.adj_thresh,  '\n Users can change cutoffs by adjusting the options "lfc_thresh" & "p.adj_thresh". \n\n', sep = ''))
  } else {
    cat(paste('Extracting Down-regulated genesets with log2FC', lfc_thresh, ' & p-value threshold ', p.adj_thresh,  '\n Users can change cutoffs by adjusting the options "lfc_thresh" & "p.adj_thresh". \n\n', sep = ''))
  }
  
  if(df_source == 'DESeq'){
    if(change == 'pos'){
      geneset = df %>% filter(log2FoldChange >= lfc_thresh & padj < p.adj_thresh) %>% pull(gene_name)
    }else if(change == 'neg'){
      geneset = df %>% filter(log2FoldChange <= -lfc_thresh & padj < p.adj_thresh) %>% pull(gene_name)
     }
  }
  
  if(df_source == 'Seurat'){
    df$gene_name = rownames(df)
    if(change == 'pos'){
      geneset = df %>% filter(avg_log2FC >= lfc_thresh & p_val_adj < p.adj_thresh) %>% pull(gene_name)
    }else if(change == 'neg'){
      geneset = df %>% filter(avg_log2FC <= -lfc_thresh & p_val_adj < p.adj_thresh) %>% pull(gene_name)
    }
  }
  
  cat("Calculating Over Representation... \n\n")
  res_go1 = gost(query = geneset, 
                 organism = organism,
                 evcodes = TRUE,
                 #sources = 'ALL',
                 correction_method= "fdr",
                 significant = TRUE,
  )
  res_go1_df = as.data.frame(res_go1$result)
  res_go1_result <- res_go1_df %>%
    dplyr::select(source, term_id, term_name, p_value, term_size:intersection_size, intersection) %>%
    dplyr::rename("p.adjust" = 4) %>%
    arrange(p.adjust)
  
  if(is.na(outtag)){
    outtag = paste('table.ORA.', df_source, '_', organism, '_', change, '_log2FC.', lfc_thresh, '_p.', p.adj_thresh, '.xlsx',sep = '')
    cat(paste('Saving results to ', outtag, '! \n Set the "outtag" option for custom file name. \n\n',sep = ''))
  } else {
    cat(paste('Saving results to ', outtag, '!\n\n',sep = ''))
  }
  
  write_xlsx(res_go1_result, outtag)
  cat("Done! \n\n")
  return(res_go1_result)
}

run_GSEA <- function(df, df_source = NA, bg_file = NA, minSize = 10, maxSize = 2000, nperm = 10000, outtag = NA){
  if(is.na(df_source)){
    stop('argument "df_source" is missing, supports Seurat or DESeq')
  }
  if(is.na(bg_file)){
    stop('argument "bg_file" is missing, please provide the directory for gmt files.')
  }
  
  cat('Extracting ranked genesets from dataframe... \n\n')
  if(df_source == 'DESeq'){
    df1 = df %>% arrange(-log2FoldChange)
    geneset = df1 %>% pull(log2FoldChange)
    names(geneset) = df1$gene_name
    }
  
  if(df_source == 'Seurat'){
    df$gene_name = rownames(df)
    df1 = df %>% arrange(-avg_log2FC)
    geneset = df1 %>% pull(avg_log2FC)
    names(geneset) = df1$gene_name
    }
  
  cat(paste('Calculating Gene Set Enrichment Scores with minSize ', minSize, ', maxSize ', maxSize, ', nperm ', nperm, '... \n', " Adjust the values by providing values in the options. \n\n",sep = ''))
  pp <- fgsea::gmtPathways(bg_file)
  fgseaRes = fgseaMultilevel(pathways = pp,
                   stats = geneset, 
                   nPermSimple = nperm,
                   minSize = minSize,
                   maxSize = maxSize,
                   )
  fgseaRes2 = fgseaRes %>% dplyr::filter(padj <= 0.05)
  fgseaRes2$leadingEdge <- vapply(fgseaRes2$leadingEdge, paste, collapse = ", ", character(1L))
  fgseaRes2 = fgseaRes2 %>% arrange(padj)
  colnames(fgseaRes2) = c("term_name", "p.value", "p.adjust", "log2err", "ES", "NES", "intersection_size", "intersection")
  
  if(is.na(outtag)){
    outtag = paste('table.GSEA.', df_source, "_minSize.", minSize, "_maxSize.", maxSize, "_nPerm.", nperm, '.xlsx', sep = '')
    cat(paste('Saving results to ', outtag, '! \n Set the "outtag" option for custom file name. \n\n',sep = ''))} else {
    cat(paste('Saving results to ', outtag, '!\n\n',sep = ''))
    }
  write_xlsx(fgseaRes2, outtag)
  
  cat('Done! \n')
  return(fgseaRes2)
}


plot_bar <- function(filename, res, isORA, n_term, width, height){
  res_term = res %>% top_n(n_term, -p.adjust)
  if (isORA){
    res_term$term_name = paste(res_term$term_name, " (", res_term$term_id, ")", sep = "")
    res_term$color = res_term$p.adjust
    legend_title = "FDR"
  }
  else{
    res_term$color = res_term$NES
    legend_title = "NES"
  }
  
  plot <- ggplot(res_term, 
                 aes(x = reorder(term_name, -log10(p.adjust)), y = -log10(p.adjust))) + 
    geom_bar(stat = "identity", position = "identity",
             color = 'black', lwd = 0.5, alpha = 1, show.legend = T,
             aes(fill = color)) +
    scale_fill_viridis_c() +
    xlab('') +
    ylab(expression(-log[10]("FDR"))) + # y-axis label
    labs(fill = legend_title) +
    coord_flip() +
    theme_classic(base_size = 12) +
    ylim(0, max(-log10(res_term$p.adjust))) +
    theme(
      aspect.ratio = 1.5,
      panel.grid.major = element_blank(),
      axis.text.x=element_text(size = 10, colour="black", vjust = 1.1, hjust = 1),
      axis.text.y = element_text(size = 10, colour = "black"),
      axis.title.x = element_text(size = 12, vjust = 0.2),
      plot.title = element_text(size = 12, face = "bold", vjust = 1.2, hjust = 0.5),
      legend.position="right",
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
    )
  ggsave(filename = filename, plot = plot, width = width, height = height)
  return(plot)
}


plot_dotplot <- function(filename, res, isORA, n_term, width, height){
  res_term = res %>% top_n(n_term, -p.adjust)
  if (isORA){
    res_term$term_name = paste(res_term$term_name, " (", res_term$term_id, ")", sep = "")
    res_term$color = res_term$p.adjust
    legend_title = "FDR"
  }
  else{
    res_term$color = res_term$NES
    res_term$size = res_term$intersection_size
    legend_title = "NES"
  }
  
  plot <- ggplot(res_term, 
                 aes(x =reorder(term_name, -log10(p.adjust)) , y = -log10(p.adjust))) + 
    geom_point(color = 'black', alpha = 1, show.legend = T, shape = 21, 
               aes(fill = color, size = intersection_size)) +
    scale_fill_viridis_c() +
    scale_size(range = c(1,3.5))+
    ylab(expression(-log[10]("FDR"))) +
    xlab('') + # y-axis label
    labs(fill = legend_title) +
    coord_flip() +
    theme_classic(base_size = 12) +
    ylim(min(-log10(res_term$p.adjust))-0.1, max(-log10(res_term$p.adjust))) +
    theme(
      aspect.ratio = 1.5,
      panel.grid.major = element_blank(),
      axis.text.x=element_text(size = 10, colour="black", vjust = 1.1, hjust = 1),
      axis.text.y = element_text(size = 10, colour = "black"),
      axis.title.x = element_text(size = 12, vjust = 0.2),
      plot.title = element_text(size = 12, face = "bold", vjust = 1.2, hjust = 0.5),
      legend.position="right",
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
    )
  ggsave(filename = filename, plot = plot, width = width, height = height)
  return(plot)
}

