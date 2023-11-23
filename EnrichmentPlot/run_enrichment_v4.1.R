library(tidyverse)
library(clusterProfiler)
library(ReactomePA)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(cowplot)
library(ggrepel)
library(writexl)
library(ggpubr)
library(enrichplot)
library(ggnewscale)
library(patchwork)

prepare_dataset <- function(df, lfc_column = NA, padj_column = NA, gene_name_column = NA, lfc_thresh = NA, padj_cutoff = NA){
  d = data.frame('SYMBOL' = df[,gene_name_column], 
                 'log2FC' = df[,lfc_column],
                 'p.adjust' = df[,padj_column])
  d = d %>% filter(!is.na(p.adjust))
  d$change = ifelse(d$p.adjust > padj_cutoff, 'None', 
                    ifelse(d$log2FC >= lfc_thresh, 'Up',
                           ifelse(d$log2FC <= -lfc_thresh, 'Down', 'None')))
  d$change = factor(d$change, levels = c('Up', 'Down', 'None'))
  return(d)
}

run_ORA <- function(df, organism = NA, minGSSize = 100, maxGSSize = 2000){
  if(organism == 'human'){
    org_db = org.Hs.eg.db
    kegg_org = 'hsa'
    WP_org = 'Homo sapiens'
    PA_org = 'human'
  } 
  else if (organism == 'mouse'){
    org_db = org.Mm.eg.db
    kegg_org = 'mmu'
    WP_org = 'Mus musculus'
    PA_org = 'mouse'
  } else {
    stop('Organism not supported. Possible options: "human", "mouse"')
  }
  
  up_deg = df %>% filter(change == 'Up') %>% pull(SYMBOL)
  down_deg = df %>% filter(change == 'Down') %>% pull(SYMBOL)
  
  cat("\nConverting all DEGs to ENTREZ ID... \n")
  ENTREZID_up <- bitr(up_deg, fromType = "SYMBOL",
                      toType = c("ENTREZID"),
                      OrgDb = org_db, drop = TRUE) %>% pull(ENTREZID)
  ENTREZID_dn <- bitr(down_deg, fromType = "SYMBOL",
                      toType = c("ENTREZID"),
                      OrgDb = org_db, drop = TRUE) %>% pull(ENTREZID)
  # Run ORA
  cat("\nRunning ORA for up-regulated DEGs... \n")
  cat(" - Running enrichGO...\n")
  go_up_res = enrichGO(
    gene = ENTREZID_up,
    OrgDb = org_db,
    keyType = "ENTREZID",
    ont = "all",
    pvalueCutoff = 0.05,
    pAdjustMethod = "fdr",
    qvalueCutoff = 0.2,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    readable = FALSE,
    pool = FALSE
  )
  
  cat(" - Running enrichKEGG... \n")
  kegg_up_res = enrichKEGG(
    gene = ENTREZID_up,
    organism = kegg_org,
    keyType = "ncbi-geneid",
    pvalueCutoff = 0.05,
    pAdjustMethod = "fdr",
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    qvalueCutoff = 0.2,
    use_internal_data = FALSE
  )
  
  cat(" - Running enrichWP... \n")
  wp_up_res = enrichWP(gene = ENTREZID_up, 
                       organism = WP_org,
                       minGSSize = minGSSize,
                       maxGSSize = maxGSSize)
  
  cat(" - Running enrichPathway... \n\n")
  pa_up_res = enrichPathway(gene= ENTREZID_up, 
                            pvalueCutoff = 0.05,
                            organism = PA_org,
                            minGSSize = minGSSize,
                            maxGSSize = maxGSSize)
  
  up_res = merge_result(list(go = go_up_res, 
                             kegg = kegg_up_res, 
                             wp = wp_up_res, 
                             reactome = pa_up_res))
  up_res_df = up_res@compareClusterResult
  up_res_df$p.adjust = p.adjust(up_res_df$pvalue, 'fdr')
  firstup <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
  }
  up_res_df$Description <- tolower(up_res_df$Description)
  up_res_df$Description <- firstup(up_res_df$Description)
  
  up_res_df$Description_ID = paste(up_res_df$Description, " (", up_res_df$ID, ")", sep = "")
  
  cat("Running ORA for down-regulated DEGs... \n")
  cat(" - Running enrichGO... \n")
  go_down_res = enrichGO(
    gene = ENTREZID_dn,
    OrgDb = org_db,
    keyType = "ENTREZID",
    ont = "all",
    pvalueCutoff = 0.05,
    pAdjustMethod = "fdr",
    qvalueCutoff = 0.2,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize
  )
  
  cat(" - Running enrichKEGG.... \n")
  kegg_down_res = enrichKEGG(
    gene = ENTREZID_dn,
    organism = kegg_org,
    keyType = "ncbi-geneid",
    pvalueCutoff = 0.05,
    pAdjustMethod = "fdr",
    minGSSize = minGSSize,
    maxGSSize = maxGSSize
  )
  
  cat(" - Running enrichWP... \n")
  wp_down_res = enrichWP(gene = ENTREZID_dn, 
                         organism = WP_org,
                         minGSSize = minGSSize,
                         maxGSSize = maxGSSize)
  
  cat(" - Running enrichPathway... \n\n")
  pa_down_res = enrichPathway(gene= ENTREZID_dn, 
                              pvalueCutoff = 0.05,
                              organism = PA_org,
                              minGSSize = minGSSize,
                              maxGSSize = maxGSSize)
  
  down_res = merge_result(list(go = go_down_res,
                               kegg = kegg_down_res, 
                               wiki = wp_down_res, 
                               reactome = pa_down_res))
  down_res_df = down_res@compareClusterResult
  down_res_df$p.adjust = p.adjust(down_res_df$pvalue, 'fdr')
  down_res_df$Description <- tolower(down_res_df$Description)
  down_res_df$Description <- firstup(down_res_df$Description)
  
  down_res_df$Description_ID = paste(down_res_df$Description, " (", down_res_df$ID, ")", sep = "")
  res = list(Up = up_res_df, Down = down_res_df, ego_up = go_up_res, ego_down = go_down_res)
  cat("ORA complete :) \n\n")
  
  return(res)
}

run_GSEA <- function(df, organism = NA, minGSSize = 100, maxGSSize = 2000){
  
  ### Convert gene format
  cat("Converting all pre-ranked genes to ENTREZ ID... \n")
  if (organism == 'human'){
    ENTREZID <- bitr(df$SYMBOL, fromType = "SYMBOL",
                   toType = c("ENTREZID"),
                   OrgDb = org.Hs.eg.db, drop = TRUE)
  } else if (organism == 'mouse'){
    ENTREZID <- bitr(df$SYMBOL, fromType = "SYMBOL",
                     toType = c("ENTREZID"),
                     OrgDb = org.Mm.eg.db, drop = TRUE)
  }
  else {
    print('Organism not supported. Possible options: "human", "mouse"')
  }
  df1 = merge(df, ENTREZID, by = "SYMBOL")
  
  ### Make gene list
  df2 = df1 %>% arrange(-log2FC)
  geneList = df2 %>% pull(log2FC)
  names(geneList) = df2$ENTREZID
  
  ########## Pathway Analysis ##########
  if (organism == "human"){
    ## GO
    cat("\nRunning GSEA for pre-ranked genes... \n")
    cat(" - Running gseGO...\n")
    ego <- gseGO(geneList     = geneList,
                 OrgDb        = org.Hs.eg.db,
                 ont          = "ALL",
                 minGSSize    = minGSSize,
                 maxGSSize    = maxGSSize,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "fdr", 
                 verbose      = FALSE)
    ## KEGG
    cat(" - Running gseKEGG...\n")
    kk <- gseKEGG(geneList     = geneList,
                  organism     = 'hsa',
                  minGSSize    = minGSSize,
                  maxGSSize    = maxGSSize,
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "fdr", 
                  verbose      = FALSE)
    ## WikiPathway
    cat(" - Running gseWP...\n")
    wp <- gseWP(geneList = geneList, 
                organism = "Homo sapiens",
                minGSSize    = minGSSize,
                maxGSSize    = maxGSSize,
                pvalueCutoff = 0.05,
                pAdjustMethod = "fdr", 
                verbose      = FALSE)
    ## Reactome
    cat(" - Running gsePathway...\n\n")
    y <- gsePathway(geneList, 
                    organism = "human",
                    minGSSize    = minGSSize,
                    maxGSSize    = maxGSSize,
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "fdr", 
                    verbose = FALSE)
  } else if (organism == "mouse") {
    ## GO
    cat("\nRunning GSEA for pre-ranked genes... \n")
    cat(" - Running gseGO...\n")
    ego <- gseGO(geneList     = geneList,
                 OrgDb        = org.Mm.eg.db,
                 ont          = "ALL",
                 minGSSize    = minGSSize,
                 maxGSSize    = maxGSSize,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "fdr", 
                 verbose      = FALSE)
    ## KEGG
    cat(" - Running gseKEGG...\n")
    kk <- gseKEGG(geneList     = geneList,
                  organism     = 'mmu',
                  minGSSize    = minGSSize,
                  maxGSSize    = maxGSSize,
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "fdr", 
                  verbose      = FALSE)
    ## WikiPathway
    cat(" - Running gseWP...\n")
    wp <- gseWP(geneList = geneList, 
                organism = "Mus musculus",
                minGSSize    = minGSSize,
                maxGSSize    = maxGSSize,
                pvalueCutoff = 0.05,
                pAdjustMethod = "fdr", 
                verbose      = FALSE)
    ## Reactome
    cat(" - Running gsePathway...\n\n")
    y <- gsePathway(geneList, 
                    organism = "mouse",
                    minGSSize    = minGSSize,
                    maxGSSize    = maxGSSize,
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "fdr", 
                    verbose = FALSE)
  } else {
    print('Organism not supported. Possible options: "human", "mouse"')
  }
  
  ## Merge all results
  merged_res = merge_result(list(go = ego, kegg = kk, wiki = wp, reactome = y))
  merged_res2 = as.data.frame(merged_res@compareClusterResult)
  
  merged_res2$p.adjust = p.adjust(merged_res2$pvalue, 'fdr')
  
  firstup <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
  }
  merged_res2$Description <- tolower(merged_res2$Description)
  merged_res2$Description <- firstup(merged_res2$Description)
  
  merged_res2$Description_ID = paste(merged_res2$Description, " (", merged_res2$ID, ")", sep = "")
  
  res_gsea = list(GSEA = merged_res2, ego_gsea = ego)
  
  cat("GSEA complete :] \n\n")
  return(res_gsea)
}

run_all <- function(df, 
                    lfc_column = NA, 
                    padj_column = NA, 
                    gene_name_column = NA, 
                    lfc_thresh = NA, 
                    padj_cutoff = NA,
                    organism = NA,
                    minGSSize = 100, 
                    maxGSSize = 2000,
                    output_file = 'results.xlsx'){
  prep_df = prepare_dataset(df, lfc_column = lfc_column, padj_column = padj_column, gene_name_column = gene_name_column, lfc_thresh = lfc_thresh, padj_cutoff = padj_cutoff)
  ORA_res = run_ORA(prep_df, organism = organism, minGSSize = minGSSize, maxGSSize = maxGSSize)
  GSEA_res = run_GSEA(prep_df, organism = organism, minGSSize = minGSSize, maxGSSize = maxGSSize)
  res_list = list(DEGs = prep_df, 
                  ORA_up = ORA_res[['Up']], 
                  ORA_down = ORA_res[['Down']], 
                  GSEA = GSEA_res[['GSEA']])
  
  wd_path = getwd()
  cat(paste0('Saving results to ', '"', wd_path, '/', output_file, '"\n'))
  write_xlsx(x = res_list, path = output_file)
  
  res_list[['ORA_ego_up']] = ORA_res[['ego_up']]
  res_list[['ORA_ego_down']] = ORA_res[['ego_down']]
  res_list[['GSEA_ego']] = GSEA_res[['ego_gsea']]
  return(res_list)
}

##### Visualization #####
plot_volcano <- function(df){
  df1 = df$DEG
  #set x-axis limit for visualization
  x_range = max(abs(df1$log2FC)) * 1.05
  
  #set gene symbol to print
  gene_up = df1 %>% filter(change == 'Up') %>% slice_min(., p.adjust,n = 10) %>% pull(SYMBOL)
  gene_down = df1 %>% filter(change == 'Down') %>% slice_min(., p.adjust, n = 10) %>% pull(SYMBOL)
  res.gene = c(gene_up, gene_down)
  
  p = ggplot(df1, aes(x = log2FC, y = -log10(p.adjust))) + 
    geom_point(aes(fill = change), 
               shape = 21, alpha = 0.75,
               na.rm = F, stroke = 0, size=2.5) + # Make dots bigger
    theme_classic(base_size = 15) + # change theme
    labs(title = paste('DEG volcano'),) + # Add a title + condition information
    xlab(expression(log[2]("Expr A" / "B"))) + # x-axis label
    ylab(expression(-log[10]("FDR"))) + # y-axis label
    geom_hline(yintercept = 1.3, colour = "grey40", linetype='dashed', size = 0.7) + # Add cutoffs
    geom_vline(xintercept = 0.5, linetype='dashed', colour = "grey40", size = 0.7) + # Add 0 lines
    geom_vline(xintercept = -0.5, linetype='dashed', colour = "grey40", size = 0.7) + # Add 0 lines
    xlim(-x_range, x_range) + 
    geom_label_repel(
      data          = subset(df1, df1$SYMBOL %in% res.gene & df1$change == 'Down'),
      nudge_x       = -x_range - subset(df1, df1$SYMBOL %in% res.gene & df1$change == 'Down')$log2FC,
      segment.size  = 0.5,
      segment.color = "grey50",
      direction     = "y",
      aes(label = SYMBOL), 
      size = 4, box.padding = 0.5, max.overlaps = Inf, fill = "white") +
    geom_label_repel(
      data          = subset(df1, df1$SYMBOL %in% res.gene & df1$change == 'Up'),
      nudge_x       = x_range - subset(df1, df1$SYMBOL %in% res.gene & df1$change == 'Up')$log2FC,
      segment.size  = 0.5,
      segment.color = "grey50",
      direction     = "y",
      aes(label = SYMBOL), 
      size = 4, box.padding = 0.5, max.overlaps = Inf, fill = "white") +
    theme(
      panel.grid.major = element_line(colour = "grey80", size = 0.5),
      panel.grid.minor = element_line(colour = "grey90", size = 0.25),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      legend.title = element_blank(),
      title = element_text(size = 12, face = 'bold'),
      # legend.position = c(0.5, 0.9),
      # legend.direction = "horizontal",
      legend.background = element_rect(fill = "white"),
      legend.margin=margin(0.1,0.1,0.1,0.1)) +
    scale_fill_manual(values = c("Up" = "#E64B35", 
                                 "Down" = "#3182bd", 
                                 "None" = "grey"))
  
  return(p)
}

plot_treeplot = function(df){
  offspring.tbl_tree_item <- getFromNamespace("offspring", "tidytree")
  assign("offspring.tbl_tree_item", offspring.tbl_tree_item, envir = .GlobalEnv)
  
  termsim_ora_up <- pairwise_termsim(df$ORA_ego_up, method = "JC")
  termsim_ora_down <- pairwise_termsim(df$ORA_ego_down, method = "JC")
  termsim_gsea <- pairwise_termsim(df$GSEA_ego, method = "JC")
  
  p1 = treeplot(termsim_ora_up, showCategory = 12, fontsize = 3.5,
                cluster.params = list(n = 4, label_words_n = 4),
                offset.params = list(tiplab = 0.2), 
                hilight.params = list(hilight = F)) 
  #+  ggtitle('ORA Up')
  
  p2 = treeplot(termsim_ora_down, showCategory = 12, fontsize = 3.5,
                cluster.params = list(n = 4, label_words_n = 4),
                offset.params = list(tiplab = 0.2), 
                hilight.params = list(hilight = F)) 
  #+  ggtitle('ORA Down')
  
  p3 = treeplot(termsim_gsea, showCategory = 12, fontsize = 3.5,
                cluster.params = list(n = 4, label_words_n = 4),
                offset.params = list(tiplab = 0.2), 
                hilight.params = list(hilight = F)) 
  #+  ggtitle('GSEA')
  
  g1 = ggarrange(p1, p2, p3, ncol=1, common.legend = TRUE, legend="right", 
                 labels = c("ORA Up", "ORA Down", "GSEA"))
  return(g1)
}

plot_bar <- function(df){
  
  ### ORA
  df_ORA1 = df[["ORA_up"]]
  df_ORA1$direction = "Up"
  df_ORA1 = df_ORA1 %>% top_n(10, -p.adjust)
  
  df_ORA2 = df[["ORA_down"]]
  df_ORA2$direction = "Dn"
  df_ORA2 = df_ORA2 %>% top_n(10, -p.adjust)
  
  df_ORA = rbind.data.frame(df_ORA1, df_ORA2)
  
  df_ORA$padj_Dir = ifelse(df_ORA$direction == "Up", -log10(df_ORA$p.adjust), log10(df_ORA$p.adjust))
  df_ORA$binary_Dir = ifelse(df_ORA$direction == "Up", 1, 0)
  
  df_ORA$Description_ID_Dir = ifelse(df_ORA$direction == "Up", paste0(df_ORA$Description_ID, "   "),
                                     paste0("   ", df_ORA$Description_ID))
  
  ## Visualize
  p1 = ggplot(df_ORA, aes(x = padj_Dir, y = reorder(Description_ID_Dir, log10(p.adjust)))) +
    geom_bar(stat = "identity", 
             color = 'black', lwd = 0.5,
             aes(fill = direction), alpha= 0.75,
             position = "identity", show.legend = F) +
    scale_fill_manual(values = c("Up" = "#E64B35", "Dn" = "#3182bd")) +
    geom_text(aes(y=Description_ID_Dir, x=0, label= Description_ID_Dir), hjust=df_ORA$binary_Dir) + 
    xlim(c(-max(-log10(df_ORA$p.adjust)), max(-log10(df_ORA$p.adjust)))) +
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
      #plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      plot.margin = margin(t = 10, b = 10))
  
  ### GSEA
  df_GSEA = df[["GSEA"]]
  df_GSEA = df_GSEA %>% top_n(20, abs(NES))
  
  df_GSEA$binary_Dir = ifelse(df_GSEA$NES > 0, 1, 0)
  
  df_GSEA$Description_ID_Dir = ifelse(df_GSEA$NES > 0, paste0(df_GSEA$Description_ID, "   "),
                                      paste0("   ", df_GSEA$Description_ID))
  
  ## Visualize
  p2 = ggplot(df_GSEA, aes(x = NES, y = reorder(Description_ID_Dir, NES))) +
    geom_bar(stat = "identity", 
             color = 'black', lwd = 0.5,
             aes(fill = NES), alpha= 1,
             position = "identity", show.legend = F) +
    scale_fill_gradient2(low = '#169194', mid = 'white', high = '#C593C2', midpoint = 0) + 
    geom_text(aes(y=Description_ID_Dir, x=0, label= Description_ID_Dir), hjust=df_GSEA$binary_Dir) + 
    xlim(c(-max(df_GSEA$NES), max(df_GSEA$NES))) +
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
      #plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      plot.margin = margin(t = 10, b = 10))
  
  g1 = plot_grid(p1,p2,ncol = 1, labels = c("ORA", "GSEA"), label_size = 15)
  
  return(g1)
}

plot_dotplot <- function(df){
  df_ORA <- rbind(df[['ORA_up']] %>% mutate(direction = 'Up'), 
                  df[['ORA_down']] %>% mutate(direction = 'Down'))
  df_GSEA <- df[['GSEA']] %>% mutate(direction = ifelse(NES > 0, 'Up', 'Down'))
  
  df_ORA$GeneRatio = sapply(df_ORA$GeneRatio, function(x) eval(parse(text=x)))
  df_ORA <- df_ORA %>% mutate(log_padj = -log10(p.adjust)) %>% 
    mutate(log_padj = ifelse(direction == 'Up',log_padj, -log_padj)) %>% 
    mutate(Description_ID_Dir = ifelse(df_ORA$direction == "Up", paste0(df_ORA$Description_ID, "   "),
                                       paste0("   ", df_ORA$Description_ID))) %>%
    mutate(binary_Dir = ifelse(direction == 'Up', 1, 0)) %>% 
    group_by(direction) %>% 
    arrange(desc(abs(log_padj))) %>%
    slice_head(n = 10) %>% ungroup() %>% 
    mutate(Description_ID = reorder(Description_ID_Dir, log_padj))
  
  df_GSEA <- df_GSEA %>% mutate(log_padj = -log10(p.adjust)) %>%
    mutate(Description_ID_Dir = ifelse(df_GSEA$direction == "Up", paste0(df_GSEA$Description_ID, "   "),
                                       paste0("   ", df_GSEA$Description_ID))) %>%
    mutate(binary_Dir = ifelse(direction == 'Up', 1, 0)) %>% 
    arrange(desc(abs(NES))) %>%
    slice_head(n = 20) %>% 
    mutate(Description_ID = reorder(Description_ID_Dir, NES))
  
  # plot ORA 
  p0_ORA <- df_ORA %>% ggplot(aes(log_padj, Description_ID_Dir)) + 
    geom_segment( aes(x=0, xend=log_padj, y=Description_ID, yend=Description_ID), color = 'black', size = 0.6) +
    geom_point(data = subset(df_ORA, direction == "Up"), aes(color = GeneRatio, size = Count)) +
    scale_color_gradient(low = '#f0c4bd', high = "#E64B35") + 
    new_scale_color()+
    geom_point(data = subset(df_ORA, direction == "Down"), aes(color = GeneRatio, size = Count)) +
    scale_color_gradient(low = '#dcf0fc', high = "#3182bd") + 
    geom_text(aes(y=Description_ID_Dir, x=0, label= Description_ID_Dir), hjust= df_ORA$binary_Dir)+
    geom_vline(xintercept = 0)+
    scale_x_continuous(limits = limits = c(-max(abs(df_ORA$log_padj)), max(abs(df_ORA$log_padj))))+
    labs(x = '-log10(FDR)', y = '')+
    theme_classic(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      axis.text.x = element_text(size = 10, colour="black", vjust = 1, hjust = 0.5),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(size = 12, vjust = 0.2),
      axis.line.y = element_blank(),
      legend.position = 'none',
      plot.title = element_text(size = 12, face = "bold", vjust = 1.2, hjust = 0.5),
      #plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      plot.margin = margin(t = 10, b = 10))
  
  clr_plt_1 <- df_ORA %>% ggplot(aes(log_padj, Description_ID_Dir)) + 
    geom_segment( aes(x=0, xend=log_padj, y=Description_ID, yend=Description_ID), color = 'black', size = 0.6) +
    geom_point(data = subset(df_ORA, direction == "Up"), aes(color = GeneRatio, size = Count)) +
    scale_color_gradient(name = 'Up \n Gene Ratio', low = '#f0c4bd', high = "#E64B35") +
    scale_size_continuous(guide = 'none')
  clr_plt_2 <- df_ORA %>% ggplot(aes(log_padj, Description_ID_Dir)) + 
    geom_segment( aes(x=0, xend=log_padj, y=Description_ID, yend=Description_ID), color = 'black', size = 0.6) +
    geom_point(data = subset(df_ORA, direction == "Down"), aes(color = GeneRatio, size = Count)) +
    scale_color_gradient(name = 'Down \n Gene Ratio', low = '#dcf0fc', high = "#3182bd") +
    scale_size_continuous(guide = 'none')
  size_plt <- df_ORA %>% ggplot(aes(log_padj, Description_ID_Dir)) + 
    geom_segment( aes(x=0, xend=log_padj, y=Description_ID, yend=Description_ID), color = 'black', size = 0.6) +
    geom_point(aes(size = Count)) +
    scale_size_continuous(name = 'Gene Count')
  clr_lgd_1 <- get_legend(clr_plt_1)
  clr_lgd_2 <- get_legend(clr_plt_2)
  size_lgd <- get_legend(size_plt)
  
  blank_p <- plot_spacer() + theme_void()
  leg12 <- plot_grid(clr_lgd_1, 
                     size_lgd,
                     blank_p,
                     nrow = 3
  )
  leg30 <- plot_grid(clr_lgd_2, blank_p,
                     blank_p, 
                     nrow = 3
  )
  leg123_ORA <- plot_grid(leg12, leg30,
                          ncol = 2
  )
  p_ORA <- plot_grid(p0_ORA + theme(plot.margin = margin(r = 15, t = 10, l = 0, b = 10)),
                     leg123_ORA + theme(plot.margin = margin(l = 15)),
                     nrow = 1,
                     align = "h",
                     axis = "t",
                     rel_widths = c(1, 0.3))
  
  # plot GSEA
  p0_GSEA <- df_GSEA %>% ggplot(aes(NES, Description_ID_Dir)) + 
    geom_segment( aes(x=0, xend=NES, y=Description_ID, yend=Description_ID), color = 'black', size = 0.6) +
    geom_point(data = subset(df_GSEA, direction == "Up"), aes(color = NES, size = log_padj)) +
    scale_color_gradient(low = '#d9c1d7', high = "#c449bd") + 
    new_scale_color()+
    geom_point(data = subset(df_GSEA, direction == "Down"), aes(color = NES, size = log_padj)) +
    scale_color_gradient(low = '#c8e7e8', high = "#169194") + 
    geom_text(aes(y=Description_ID_Dir, x=0, label= Description_ID_Dir), hjust= df_GSEA$binary_Dir)+
    geom_vline(xintercept = 0)+
    scale_x_continuous(limits = c(-max(df_GSEA$NES), max(df_GSEA$NES)))+
    labs(x = 'NES', y = '')+
    theme_classic(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      axis.text.x = element_text(size = 10, colour="black", vjust = 1, hjust = 0.5),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(size = 12, vjust = 0.2),
      axis.line.y = element_blank(),
      legend.position = 'none',
      plot.title = element_text(size = 12, face = "bold", vjust = 1.2, hjust = 0.5),
      #plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      plot.margin = margin(t = 10, b = 10))
  
  clr_plt_1 <- df_GSEA %>% ggplot(aes(NES, Description_ID_Dir)) + 
    geom_segment( aes(x=0, xend=log_padj, y=Description_ID, yend=Description_ID), color = 'black', size = 0.6) +
    geom_point(data = subset(df_GSEA, direction == "Up"), aes(color = NES, size = log_padj)) +
    scale_color_gradient(name = 'Up NES', low = '#d9c1d7', high = "#c449bd") +
    scale_size_continuous(guide = 'none')
  clr_plt_2 <- df_GSEA %>% ggplot(aes(NES, Description_ID_Dir)) + 
    geom_segment( aes(x=0, xend=log_padj, y=Description_ID, yend=Description_ID), color = 'black', size = 0.6) +
    geom_point(data = subset(df_GSEA, direction == "Down"), aes(color = NES, size = log_padj)) +
    scale_color_gradient(name = 'Down NES', low = '#c8e7e8', high = "#169194") +
    scale_size_continuous(guide = 'none')
  size_plt <- df_GSEA %>% ggplot(aes(NES, Description_ID_Dir)) + 
    geom_segment( aes(x=0, xend=NES, y=Description_ID, yend=Description_ID), color = 'black', size = 0.6) +
    geom_point(aes(size = log_padj)) +
    scale_size_continuous(name = '-log10(FDR)')
  
  clr_lgd_1 <- get_legend(clr_plt_1)
  clr_lgd_2 <- get_legend(clr_plt_2)
  size_lgd <- get_legend(size_plt)
  
  blank_p <- plot_spacer() + theme_void()
  leg12 <- plot_grid(clr_lgd_1, 
                     size_lgd,
                     blank_p,
                     nrow = 3
  )
  leg30 <- plot_grid(clr_lgd_2, blank_p,
                     blank_p, 
                     nrow = 3
  )
  leg123_GSEA <- plot_grid(leg12, leg30,
                           ncol = 2
  )
  p_GSEA <- plot_grid(p0_GSEA + theme(plot.margin = margin(r = 15, t = 10, l = 0, b = 10)),
                      leg123_GSEA + theme(plot.margin = margin(l = 15)),
                      nrow = 1,
                      align = "h",
                      axis = "t",
                      rel_widths = c(1, 0.3))
  p <- plot_grid(p_ORA, p_GSEA, ncol = 1, labels = c('ORA', 'GSEA'))
  return(p)
}

plot_all <- function(df, output_file = 'results.pdf', width = 22.5, height = 20){
  p1 = plot_volcano(df)
  p2 = plot_treeplot(df)
  p3 = plot_bar(df)
  p4 = plot_dotplot(df)
  
  g0 = plot_grid(p1, p2, nrow = 1, rel_widths = c(1, 1.5))
  g1 = plot_grid(p3, p4, nrow = 1, rel_widths = c(1,1.25))
  g2 = plot_grid(g0, g1, nrow = 2, rel_heights = c(1, 1.5))
  
  ggsave(g2, file= output_file, width = width, height = height)
  return(g2)
}

