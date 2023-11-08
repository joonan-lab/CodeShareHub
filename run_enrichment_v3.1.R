library(tidyverse)
library(clusterProfiler)
library(ReactomePA)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(cowplot)
library(ggrepel)
library(writexl)

prepare_dataset <- function(df, lfc = NA, padj = NA, gene_name = NA, lfc_thresh = NA, padj_cutoff = NA){
  d = data.frame('SYMBOL' = df[,gene_name], 
                 'log2FC' = df[,lfc],
                 'p.adjust' = df[,padj])
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
  ENTREZID <- bitr(df$SYMBOL, fromType = "SYMBOL",
                   toType = c("ENTREZID"),
                   OrgDb = org.Hs.eg.db, drop = TRUE)
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
                    lfc = NA, 
                    padj = NA, 
                    gene_name = NA, 
                    lfc_thresh = NA, 
                    padj_cutoff = NA,
                    organism = NA,
                    minGSSize = 100, 
                    maxGSSize = 2000,
                    out_dir = 'results.xlsx'){
  prep_df = prepare_dataset(df, lfc = lfc, padj = padj, gene_name = gene_name, lfc_thresh = lfc_thresh, padj_cutoff = padj_cutoff)
  ORA_res = run_ORA(prep_df, organism = organism, minGSSize = minGSSize, maxGSSize = maxGSSize)
  GSEA_res = run_GSEA(prep_df, organism = organism, minGSSize = minGSSize, maxGSSize = maxGSSize)
  res_list = list(DEGs = prep_df, 
                  ORA_up = ORA_res[['Up']], 
                  ORA_down = ORA_res[['Down']], 
                  GSEA = GSEA_res[['GSEA']])
  
  wd_path = getwd()
  cat(paste0('Saving results to ', '"', wd_path, '/', out_dir, '"\n'))
  write_xlsx(x = res_list, path = out_dir)
  
  res_list[['ORA_ego_up']] = ORA_res[['ego_up']]
  res_list[['ORA_ego_down']] = ORA_res[['ego_down']]
  res_list[['GSEA_ego']] = GSEA_res[['ego_gsea']]
  return(res_list)
}

##### Visualization #####
plot_volcano <- function(df){
  
  #set x-axis limit for visualization
  x_range = max(abs(df$log2FC)) * 1.05
  
  #set gene symbol to print
  gene_up = df %>% filter(change == 'Up') %>% slice_min(., p.adjust,n = 10) %>% pull(SYMBOL)
  gene_down = df %>% filter(change == 'Down') %>% slice_min(., p.adjust, n = 10) %>% pull(SYMBOL)
  res.gene = c(gene_up, gene_down)
  
  p = ggplot(df, aes(x = log2FC, y = -log10(p.adjust))) + 
    geom_point(aes(fill = change), 
               shape = 21, alpha = 0.75,
               na.rm = F, stroke = 0, size=2.5) + # Make dots bigger
    theme_classic(base_size = 15) + # change theme
    #labs(title = paste('DEG analysis:', condition),) + # Add a title + condition information
    xlab(expression(log[2]("Expr A" / "B"))) + # x-axis label
    ylab(expression(-log[10]("P"))) + # y-axis label
    geom_hline(yintercept = 1.3, colour = "grey40", linetype='dashed', size = 0.7) + # Add cutoffs
    geom_vline(xintercept = 0.5, linetype='dashed', colour = "grey40", size = 0.7) + # Add 0 lines
    geom_vline(xintercept = -0.5, linetype='dashed', colour = "grey40", size = 0.7) + # Add 0 lines
    xlim(-x_range, x_range) + 
    geom_label_repel(
      data          = subset(df, df$SYMBOL %in% res.gene & df$change == 'Down'),
      nudge_x       = -x_range - subset(df, df$SYMBOL %in% res.gene & df$change == 'Down')$log2FC,
      segment.size  = 0.5,
      segment.color = "grey50",
      direction     = "y",
      aes(label = SYMBOL), 
      size = 4, box.padding = 0.5, max.overlaps = Inf, fill = "white") +
    geom_label_repel(
      data          = subset(df, df$SYMBOL %in% res.gene & df$change == 'Up'),
      nudge_x       = x_range - subset(df, df$SYMBOL %in% res.gene & df$change == 'Up')$log2FC,
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
      title = element_text(size = 12),
      legend.position = c(0.5, 0.9),
      legend.direction = "horizontal",
      legend.background = element_rect(fill = "white"),
      legend.margin=margin(0.1,0.1,0.1,0.1)) +
    scale_fill_manual(values = c("Up" = "#E64B35", 
                                 "Down" = "#3182bd", 
                                 "None" = "grey"))
  
  return(p)
}



