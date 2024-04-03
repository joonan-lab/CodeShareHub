library(tidyverse)
library(clusterProfiler)
library(ReactomePA)
library(annotate)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(cowplot)
library(ggrepel)
library(writexl)
library(ggpubr)
library(enrichplot)
library(ggnewscale)
library(patchwork)
library(ggrastr)

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

######## Run Enrich Test ########
## Prepare Dataset 
prepare_dataset <- function(df, 
                            lfc_column = NA, 
                            padj_column = NA, 
                            gene_name_column = NA, 
                            lfc_thresh = 0.1, 
                            padj_cutoff = 0.05){
  print("Prepare Dataset")
  d = data.frame('SYMBOL' = df[,gene_name_column], 
                 'log2FC' = df[,lfc_column],
                 'p.adjust' = df[,padj_column])
  colnames(d) = c("SYMBOL", "log2FC", "p.adjust")
  d = d %>% filter(!is.na(p.adjust))
  d$change = ifelse(d$p.adjust > padj_cutoff, 'None', 
                    ifelse(d$log2FC >= lfc_thresh, 'Up',
                           ifelse(d$log2FC <= -lfc_thresh, 'Down', 'None')))
  d$change = factor(d$change, levels = c('Up', 'Down', 'None'))
  return(d)
}

## Run ORA
run_ORA <- function(df, 
                    organism = NA, 
                    minGSSize = 100, 
                    maxGSSize = 1000){
  print("Run ORA")
  
  ## Setting & Gene conversion.
  if(organism == 'human'){
    org_db = org.Hs.eg.db
    #kegg_org = 'hsa'
    #WP_org = 'Homo sapiens'
    PA_org = 'human'
  } 
  else if (organism == 'mouse'){
    org_db = org.Mm.eg.db
    #kegg_org = 'mmu'
    #WP_org = 'Mus musculus'
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
  # go_up_res_all = enrichGO(
  #   gene = ENTREZID_up,
  #   OrgDb = org_db,
  #   keyType = "ENTREZID",
  #   ont = "ALL",
  #   pvalueCutoff = 0.05,
  #   pAdjustMethod = "fdr",
  #   minGSSize = minGSSize,
  #   maxGSSize = maxGSSize,
  #   readable = TRUE,
  #   pool = FALSE
  # )
  
  go_up_res_bp = enrichGO(
    gene = ENTREZID_up,
    OrgDb = org_db,
    keyType = "ENTREZID",
    ont = "BP",
    pvalueCutoff = 0.05,
    pAdjustMethod = "fdr",
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    readable = TRUE
    )
  
  # go_up_res_mf = enrichGO(
  #   gene = ENTREZID_up,
  #   OrgDb = org_db,
  #   keyType = "ENTREZID",
  #   ont = "MF",
  #   pvalueCutoff = 0.05,
  #   pAdjustMethod = "fdr",
  #   minGSSize = minGSSize,
  #   maxGSSize = maxGSSize,
  #   readable = TRUE
  # )
  
  # go_up_res_cc = enrichGO(
  #   gene = ENTREZID_up,
  #   OrgDb = org_db,
  #   keyType = "ENTREZID",
  #   ont = "CC",
  #   pvalueCutoff = 0.05,
  #   pAdjustMethod = "fdr",
  #   minGSSize = minGSSize,
  #   maxGSSize = maxGSSize,
  #   readable = TRUE
  # )
  
  # cat(" - Running enrichPathway... \n\n")
  # pa_up_res = enrichPathway(
  #   gene= ENTREZID_up,
  #   organism = PA_org,
  #   pvalueCutoff = 0.05,
  #   pAdjustMethod = "fdr",
  #   minGSSize = minGSSize,
  #   maxGSSize = maxGSSize,
  #   readable = TRUE
  #   )
  
  # cat(" - Running enrichKEGG... \n")
  # kegg_up_res = enrichKEGG(
  #   gene = ENTREZID_up,
  #   organism = kegg_org,
  #   keyType = "ncbi-geneid",
  #   pvalueCutoff = 0.05,
  #   pAdjustMethod = "fdr",
  #   minGSSize = minGSSize,
  #   maxGSSize = maxGSSize,
  #   qvalueCutoff = 0.2,
  #   use_internal_data = FALSE
  # )

  # cat(" - Running enrichWP... \n")
  # wp_up_res = enrichWP(gene = ENTREZID_up,
  #                      organism = WP_org,
  #                      minGSSize = minGSSize,
  #                      maxGSSize = maxGSSize)

  # up_res = merge_result(list(
  #   go_all = go_up_res_all,
  #   go_bp = go_up_res_bp,
  #   go_cc = go_up_res_cc,
  #   go_mf = go_up_res_mf,
  #   # kegg = kegg_up_res,
  #   # wp = wp_up_res,
  #   reactome = pa_up_res
  #   ))
  # up_res_df = up_res@compareClusterResult
  
  up_res_list = list(
    # up_go_all = go_up_res_all,
    up_go_bp = go_up_res_bp
    # up_go_cc = go_up_res_cc,
    # # up_go_mf = go_up_res_mf,
    # up_reactome = pa_up_res
    )

  cat("Running ORA for down-regulated DEGs... \n")
  cat(" - Running enrichGO... \n")
  # go_down_res_all = enrichGO(
  #   gene = ENTREZID_dn,
  #   OrgDb = org_db,
  #   keyType = "ENTREZID",
  #   ont = "ALL",
  #   pvalueCutoff = 0.05,
  #   pAdjustMethod = "fdr",
  #   minGSSize = minGSSize,
  #   maxGSSize = maxGSSize,
  #   readable = TRUE,
  #   pool = FALSE
  # )
  
  go_down_res_bp = enrichGO(
    gene = ENTREZID_dn,
    OrgDb = org_db,
    keyType = "ENTREZID",
    ont = "BP",
    pvalueCutoff = 0.05,
    pAdjustMethod = "fdr",
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    readable = TRUE
  )
  
  # go_down_res_mf = enrichGO(
  #   gene = ENTREZID_dn,
  #   OrgDb = org_db,
  #   keyType = "ENTREZID",
  #   ont = "MF",
  #   pvalueCutoff = 0.05,
  #   pAdjustMethod = "fdr",
  #   minGSSize = minGSSize,
  #   maxGSSize = maxGSSize,
  #   readable = TRUE
  # )
  
  # go_down_res_cc = enrichGO(
  #   gene = ENTREZID_dn,
  #   OrgDb = org_db,
  #   keyType = "ENTREZID",
  #   ont = "CC",
  #   pvalueCutoff = 0.05,
  #   pAdjustMethod = "fdr",
  #   minGSSize = minGSSize,
  #   maxGSSize = maxGSSize,
  #   readable = TRUE
  # )
  
  # cat(" - Running enrichKEGG.... \n")
  # kegg_down_res = enrichKEGG(
  #   gene = ENTREZID_dn,
  #   organism = kegg_org,
  #   keyType = "ncbi-geneid",
  #   pvalueCutoff = 0.05,
  #   pAdjustMethod = "fdr",
  #   minGSSize = minGSSize,
  #   maxGSSize = maxGSSize
  # )

  # cat(" - Running enrichWP... \n")
  # wp_down_res = enrichWP(gene = ENTREZID_dn,
  #                        organism = WP_org,
  #                        minGSSize = minGSSize,
  #                        maxGSSize = maxGSSize)

  # cat(" - Running enrichPathway... \n\n")
  # pa_down_res = enrichPathway(gene= ENTREZID_dn,
  #                             pvalueCutoff = 0.05,
  #                             organism = PA_org,
  #                             pAdjustMethod = "fdr",
  #                             minGSSize = minGSSize,
  #                             maxGSSize = maxGSSize,
  #                             readable = TRUE)
  
  down_res_list = list(
    # down_go_all = go_down_res_all,
    down_go_bp = go_down_res_bp
    # down_go_cc = go_down_res_cc,
    # down_go_mf = go_down_res_mf,
    # down_reactome = pa_down_res
    )
  
  # # Function to convert gene IDs to gene symbols
  # org = ifelse(organism == 'human', "org.Hs.eg", "org.Mm.eg")
  # 
  # convert_gene_ids <- function(gene_ids) {
  #   gene_id_list <- strsplit(gene_ids, "/")
  #   entrez_ids <- unlist(gene_id_list)
  #   gene_symbols <- getSYMBOL(entrez_ids, data = org)
  #   result <- paste(gene_symbols, collapse = ", ")
  #   return(result)
  # }
  # 
  # up_res_df$gene_symbols <- sapply(up_res_df$geneID, convert_gene_ids)  
  # down_res_df$gene_symbols <- sapply(down_res_df$geneID, convert_gene_ids)
  
  # up_res_df_bp = up_res_df[up_res_df$Cluster == "go_bp",]
  # down_res_df_bp = down_res_df[down_res_df$Cluster == "go_bp",]
  
  res_ORA_list = append(
    up_res_list, 
    down_res_list,
    )
  
  cat("ORA complete :) \n\n")
  
  return(res_ORA_list)
}

## Run GSEA
run_GSEA <- function(df, 
                     organism = NA, 
                     minGSSize = 100, 
                     maxGSSize = 1000){
  print("Run GSEA")
  
  ## Setting & Gene conversion.
  if(organism == 'human'){
    org_db = org.Hs.eg.db
    # kegg_org = 'hsa'
    # WP_org = 'Homo sapiens'
    PA_org = 'human'
  } 
  else if (organism == 'mouse'){
    org_db = org.Mm.eg.db
    # kegg_org = 'mmu'
    # WP_org = 'Mus musculus'
    PA_org = 'mouse'
  } else {
    stop('Organism not supported. Possible options: "human", "mouse"')
  }
  
  ### Convert gene format
  cat("Converting all pre-ranked genes to ENTREZ ID... \n")
  
  ENTREZID <- bitr(df$SYMBOL, fromType = "SYMBOL",
                    toType = c("ENTREZID"),
                    OrgDb = org_db, drop = TRUE) 
  
  df1 = merge(df, ENTREZID, by = "SYMBOL")
  
  ### Make gene list
  df2 = df1 %>% arrange(-log2FC)
  geneList = df2 %>% pull(log2FC)
  names(geneList) = df2$ENTREZID
  
  ########## Pathway Analysis ##########
  ## GO
  cat("\nRunning GSEA for pre-ranked genes... \n")
  
  cat(" - Running gseGO...\n")
  # ego_all <- gseGO(geneList     = geneList,
  #                 OrgDb        = org_db,
  #                 ont          = "ALL",
  #                 minGSSize    = minGSSize,
  #                 maxGSSize    = maxGSSize,
  #                 pvalueCutoff = 0.05,
  #                 pAdjustMethod = "fdr", 
  #                 verbose      = FALSE,
  #                 eps = 0)
  
  ego_bp <- gseGO(geneList     = geneList,
               OrgDb        = org_db,
               ont          = "BP",
               minGSSize    = minGSSize,
               maxGSSize    = maxGSSize,
               pvalueCutoff = 0.05,
               pAdjustMethod = "fdr", 
               verbose      = FALSE,
               eps = 0)
  
  # ego_mf <- gseGO(geneList     = geneList,
  #                 OrgDb        = org_db,
  #                 ont          = "MF",
  #                 minGSSize    = minGSSize,
  #                 maxGSSize    = maxGSSize,
  #                 pvalueCutoff = 0.05,
  #                 pAdjustMethod = "fdr", 
  #                 verbose      = FALSE,
  #                 eps = 0)
  # 
  # ego_cc <- gseGO(geneList     = geneList,
  #                 OrgDb        = org_db,
  #                 ont          = "CC",
  #                 minGSSize    = minGSSize,
  #                 maxGSSize    = maxGSSize,
  #                 pvalueCutoff = 0.05,
  #                 pAdjustMethod = "fdr", 
  #                 verbose      = FALSE,
  #                 eps = 0)
  
  ## Reactome
  # cat(" - Running gsePathway...\n\n")
  # y <- gsePathway(geneList,
  #                 organism = PA_org,
  #                 minGSSize    = minGSSize,
  #                 maxGSSize    = maxGSSize,
  #                 pvalueCutoff = 0.05,
  #                 pAdjustMethod = "fdr",
  #                 verbose = FALSE,
  #                 eps = 0)
  
  ## Merge all results
  gsea_res_list = list(
    # gsea_go_all = ego_all,
                       gsea_go_bp = ego_bp
                       # gsea_go_cc = ego_cc,
                       # gsea_go_mf = ego_mf,
                       # gsea_reactome = y
                      )
  
  cat("GSEA complete :] \n\n")
  return(gsea_res_list)
}

## Run above codes at one and return list of results (res_list_df, res_list)
run_all <- function(df, 
                    lfc_column = NA, 
                    padj_column = NA, 
                    gene_name_column = NA, 
                    lfc_thresh = 0.1, 
                    padj_cutoff = 0.05,
                    organism = NA,
                    minGSSize = 100, 
                    maxGSSize = 1000,
                    output_file = 'results.xlsx'){
  
  ## Prepare Dataset
  prep_df = prepare_dataset(df, lfc_column = lfc_column, padj_column = padj_column, gene_name_column = gene_name_column, lfc_thresh = lfc_thresh, padj_cutoff = padj_cutoff)
  
  ## Run ORA
  ORA_res = run_ORA(prep_df, organism = organism, minGSSize = minGSSize, maxGSSize = maxGSSize)
  
  ## Run GSEA
  GSEA_res = run_GSEA(prep_df, organism = organism, minGSSize = minGSSize, maxGSSize = maxGSSize)
  
  res_list = append(ORA_res, GSEA_res)
  
  res_df_list = lapply(res_list, function(x){x@result})  
  names(res_df_list) = toupper(names(res_list))
  res_df_list[["DEGs"]] = prep_df
  
  org = ifelse(organism == 'human', "org.Hs.eg", "org.Mm.eg")
  
  convert_gene_ids <- function(gene_ids) {
    gene_id_list <- strsplit(gene_ids, "/")
    entrez_ids <- unlist(gene_id_list)
    gene_symbols <- getSYMBOL(entrez_ids, data = org)
    result <- paste(gene_symbols, collapse = "/")
    return(result)
  }
  
  # res_df_list[["GSEA_GO_ALL"]]$core_enrichment <- sapply(res_df_list[["GSEA_GO_ALL"]]$core_enrichment, convert_gene_ids)
  res_df_list[["GSEA_GO_BP"]]$core_enrichment <- sapply(res_df_list[["GSEA_GO_BP"]]$core_enrichment, convert_gene_ids)  
  # res_df_list[["GSEA_GO_CC"]]$core_enrichment <- sapply(res_df_list[["GSEA_GO_CC"]]$core_enrichment, convert_gene_ids)  
  # res_df_list[["GSEA_GO_MF"]]$core_enrichment <- sapply(res_df_list[["GSEA_GO_MF"]]$core_enrichment, convert_gene_ids)  
  # res_df_list[["GSEA_REACTOME"]]$core_enrichment <- sapply(res_df_list[["GSEA_REACTOME"]]$core_enrichment, convert_gene_ids)  
  
  wd_path = getwd()
  cat(paste0('Saving results to ', '"', wd_path, '/', output_file, '"\n'))
  write_xlsx(x = res_df_list, path = output_file)

  res = list(res_df_list = res_df_list, res_list = res_list)
  
  return(res)
}


######## Visualization ########
## Plot volcano
plot_volcano <- function(res_df_list, 
                         padj_val = padj_cutoff, 
                         lfc_val = lfc_thresh){
  
  df1 = res_df_list$DEGs
  
  # Set x-axis limit for visualization
  x_range = max(abs(df1$log2FC)) * 1.05
  
  # Set gene symbol to print
  gene_up = df1 %>% filter(change == 'Up') %>% slice_min(., p.adjust,n = 15) %>% pull(SYMBOL)
  gene_down = df1 %>% filter(change == 'Down') %>% slice_min(., p.adjust, n = 15) %>% pull(SYMBOL)
  res.gene = c(gene_up, gene_down)
  
  # Plot volcano
  p = ggplot(df1, aes(x = log2FC, y = -log10(p.adjust))) + 
    geom_point_rast(aes(fill = change), 
               shape = 21, alpha = 0.75,
               na.rm = F, stroke = 0, size=2.5) + # Make dots bigger
    theme_classic(base_size = 15) + # change theme
    labs(title = paste('DEG volcano'),) + # Add a title + condition information
    xlab(expression(log[2]("Expr A" / "B"))) + # x-axis label
    ylab(expression(-log[10]("pvalue"))) + # y-axis label
    geom_hline(yintercept = -log10(as.numeric(padj_val)) - 0.01, colour = "grey40", linetype='dashed', size = 0.7) + # Add cutoffs
      geom_vline(xintercept = as.numeric(lfc_val), linetype='dashed', colour = "grey40", size = 0.7) + # Add 0 lines
     geom_vline(xintercept = -as.numeric(lfc_val), linetype='dashed', colour = "grey40", size = 0.7) + # Add 0 lines
    xlim(-x_range, x_range) + 
    geom_label_repel(
      data          = subset(df1, df1$SYMBOL %in% res.gene & df1$change == 'Down') %>% top_n(10, -p.adjust),
      nudge_x       = -x_range - subset(df1, df1$SYMBOL %in% res.gene & df1$change == 'Down')$log2FC,
      segment.size  = 0.5,
      segment.color = "grey50",
      direction     = "y",
      aes(label = SYMBOL), 
      size = 4, box.padding = 0.5, max.overlaps = Inf, fill = "white") +
    geom_label_repel(
      data          = subset(df1, df1$SYMBOL %in% res.gene & df1$change == 'Up')  %>% top_n(10, -p.adjust),
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

## Plot treeplot using GO:BP data
plot_treeplot = function(res_list){
  offspring.tbl_tree_item <- getFromNamespace("offspring", "tidytree")
  assign("offspring.tbl_tree_item", offspring.tbl_tree_item, envir = .GlobalEnv)
  
  termsim_ora_up <- pairwise_termsim(res_list$up_go_all, method = "JC")
  termsim_ora_down <- pairwise_termsim(res_list$down_go_all, method = "JC")
  termsim_gsea <- pairwise_termsim(res_list$gsea_go_all, method = "JC")
  
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
                 labels = c("  ORA Up", "ORA Down", "  GSEA"))
  return(g1)
}

## Plot barplot
plot_bar <- function(res_df, plot_title){
  p = NULL
  
  if (nrow(res_df)==0){
    return (p)
    
  } else{
    ### ORA
    if (str_detect(plot_title, 'GO')){
      p = ggplot(res_df, aes(x = padj_Dir, y = reorder(Description_ID_Dir, padj_Dir))) +
        geom_bar(stat = "identity", 
                 color = 'black', lwd = 0.5,
                 aes(fill = direction), alpha= 0.75,
                 position = "identity", show.legend = F) +
        scale_fill_manual(values = c("Up" = "#E64B35", "Down" = "#3182bd")) +
        geom_text(aes(y=Description_ID_Dir, x=0, label= Description_ID_Dir), hjust=res_df$binary_Dir) + 
        xlim(c(-max(-log10(res_df$p.adjust)), max(-log10(res_df$p.adjust)))) +
        #coord_fixed(ratio = 0.15) +
        labs(x = '-log10(FDR)', y = '', title=plot_title) +
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
      
    }else if (str_detect(plot_title, 'GSEA')){
      ## Visualize
      p = ggplot(res_df, aes(x = NES, y = reorder(Description_ID_Dir, NES))) +
        geom_bar(stat = "identity", 
                 color = 'black', lwd = 0.5,
                 aes(fill = NES), alpha= 1,
                 position = "identity", show.legend = F) +
        scale_fill_gradient2(low = '#169194', mid = 'white', high = '#C593C2', midpoint = 0) + 
        geom_text(aes(y=Description_ID_Dir, x=0, label= Description_ID_Dir), hjust=res_df$binary_Dir) + 
        xlim(c(-max(abs(res_df$NES)), max(abs(res_df$NES)))) +
        #coord_fixed(ratio = 0.15) +
        labs(x = 'NES', y = '', title=plot_title) +
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
    }
    
    
    #g1 = plot_grid(p1,p2,ncol = 1, labels = c("ORA", "GSEA"),
    #               label_x = 0, label_y = 1,
    #               hjust = -0.5, vjust = -0.5)
    
    return(p)
  }
  
}


## Plot all in one
plot_all <- function(res, 
                     output_file = 'results.pres_list', 
                     width = 22.5, 
                     height = 20, 
                     padj_cutoff = 0.05, 
                     lfc_thresh = 0.1){
  res_list_df = res[['res_df_list']]
  res_list = res[['res_list']]
  
  p1 = plot_volcano(res_list_df, padj_val = padj_cutoff, lfc_val = lfc_thresh)
  p2 = plot_treeplot(res_list)
  #bar_df_list <- make_bar_df(res_list_df)
  ### make df
  
  ### ORA
  make_df <- function(res_list_df, db_source, reactome){
    res_df = data.frame()
    
    if (db_source == 'GO'){
      if (reactome == TRUE){
        up_name = "UP_REACTOME"
        down_name = "DOWN_REACTOME"
        name = "GO_REACTOME"
      }else{
        up_name = "UP_GO_ALL"
        down_name = "DOWN_GO_ALL"
        name = "GO_ALL"
      }
      res_go_up = res_list_df[[up_name]]
      res_go_up$direction = 'Up'
      res_go_up = res_go_up %>% top_n(10, -p.adjust)
      
      res_go_down = res_list_df[[down_name]]
      res_go_down$direction = 'Down'
      res_go_down = res_go_down %>% top_n(10, -p.adjust)
      
      res_df = rbind.data.frame(res_go_up, res_go_down)
      
      if (nrow(res_df)!=0){
        res_df$padj_Dir = ifelse(res_df$direction == "Up", -log10(res_df$p.adjust), log10(res_df$p.adjust))
        res_df$binary_Dir = ifelse(res_df$direction == "Up", 1, 0)
        
        res_df$Description_ID = paste0(res_df$Description, " (", res_df$ID, ")")
        res_df$Description_ID_Dir = ifelse(res_df$direction == "Up", paste0(firstup(res_df$Description_ID), "   "),
                                           paste0("   ", firstup(res_df$Description_ID)))
      }
    }else if(db_source=='GSEA'){
      if (reactome == TRUE){
        name = "GSEA_REACTOME"
      }else{
        name = "GSEA_GO_ALL"
      }
      ### GSEA
      res_df = res_list_df[[name]]
      res_df = res_df %>% top_n(20, -log10(p.adjust))
      
      if (nrow(res_df)!=0){
        res_df$direction = ifelse(res_df$NES > 0, 'Up', 'Down')
        res_df$padj_Dir = ifelse(res_df$direction == "Up", -log10(res_df$p.adjust), log10(res_df$p.adjust))
        res_df$binary_Dir = ifelse(res_df$direction == "Up", 1, 0)
        
        res_df$Description_ID = paste0(res_df$Description, " (", res_df$ID, ")")
        res_df$Description_ID_Dir = ifelse(res_df$direction == "Up", paste0(firstup(res_df$Description_ID), "   "),
                                           paste0("   ", firstup(res_df$Description_ID)))
      }
    }
    
    return (res_df)
  }
  
  ####
  bar_df_list = list()
  bar_df_list[['GO_ALL']] = make_df(res_list_df, db_source='GO', reactome=FALSE)
  bar_df_list[['GO_REACTOME']] = make_df(res_list_df, db_source='GO', reactome=TRUE)
  bar_df_list[['GSEA']] = make_df(res_list_df, db_source='GSEA', reactome=FALSE)
  bar_df_list[['GSEA_REACTOME']] = make_df(res_list_df, db_source='GSEA', reactome=TRUE)
  
  p_bar_list = lapply(names(bar_df_list), function(x){plot_bar(bar_df_list[[x]], plot_title=x)})
  #names(p_bar_list) = names(bar_df_list)
  #p4 = plot_dotplot(res_list)
  
  g0 = plot_grid(plotlist=list(p1, p2), nrow = 1, rel_widths = c(1, 1.5))
  g1 = plot_grid(plotlist=p_bar_list, nrow = 2, rel_widths = c(1,1))
  g2 = plot_grid(plotlist=list(g0, g1), nrow = 2, rel_heights = c(1, 1.5))
  
  #ggsave(g2, file= output_file, width = width, height = height)
  print("Save result plots...")
  pdf(paste0(output_file), width = width, height = height)
  print(g2)
  dev.off()
  print("Done!")
  
  return(g2)
}


