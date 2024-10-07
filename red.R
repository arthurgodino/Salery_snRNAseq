# Analysis of red cells only

library(Matrix)
library(Seurat)
library(data.table)
library(dplyr)
library(patchwork)
library(tidyverse)
library(ggpubr)
library(xlsx)
library(cowplot)
library(gplots)

integrated <- readRDS(file = "/Users/arthur/Desktop/Marine_integrated.rds")

#Distrib of Arc expression
integrated$ArcRNA <- GetAssayData(object = integrated, assay = "RNA", slot = "counts")["Arc",]
integrated$ArcSCT <- GetAssayData(object = integrated, assay = "SCT", slot = "counts")["Arc",]
p.ArcRNA <- integrated@meta.data %>%
      group_by(cpp, ArcRNA) %>%
      summarize(n=n()) %>%
ggplot(aes(x = ArcRNA, y = n, fill = cpp)) + 
  geom_col(position = position_dodge(), width = .25) +
  scale_fill_manual(values = c("grey30", "dodgerblue")) +
  scale_x_continuous(breaks = 0:10, limits = c(NA,10)) +
  #scale_y_log10(expand = c(0,0), limits = c(1,10000), breaks = c(1,5,10,50,100,500,1000,5000,10000)) +
  theme_classic()

p.ArcSCT <- integrated@meta.data %>%
  group_by(cpp, ArcSCT) %>%
  summarize(n=n()) %>%
  ggplot(aes(x = ArcSCT, y = n, fill = cpp)) + 
  geom_col(position = position_dodge(), width = .25) +
  scale_fill_manual(values = c("grey30", "dodgerblue")) +
  scale_x_continuous(breaks = 0:10, limits = c(NA,10)) +
  scale_y_log10(expand = c(0,0)) +
  theme_classic()

ggsave(cowplot::plot_grid(p.ArcRNA, p.ArcSCT), file = "/Users/agodino/Desktop/Arc.distrib.pdf", device = pdf, width = 18, height = 8, units = "in")


############################################################################################################################################################
############################################################################################################################################################

# Table of number of nuclei by cluster by conditions by activation status
n.cells <- table(integrated$cluster_pairings_cpp_sorted_Arc) %>% as.data.frame() %>% setNames(c("cluster_pairings_cpp_sorted_Arc", "n")) %>%
  separate(col = "cluster_pairings_cpp_sorted_Arc", into = c("cluster", "pairings", "cpp", "sorted", "Arc.act"), sep = "_") %>%
  pivot_wider(id_cols = c("pairings","cpp", "sorted", "Arc.act"), names_from = cluster, values_from = n) %>%
  arrange(pairings , cpp, sorted, Arc.act) %>%
  replace(is.na(.), 0)
write.csv(n.cells, file = "/Users/arthur/Desktop/n_nuclei_per_cluster_and_condition_and_activation.csv", row.names = F)

# UMAP of red cells, sep by Test or HC
integrated$pairings_Arc <- paste(integrated$pairings, integrated$Arc.expr, sep = "_")
Idents(integrated) <- "cpp"
integrated.HC <- subset(integrated, idents = "HC")
integrated.Test <- subset(integrated, idents = "Test")

Idents(integrated.HC) <- "pairings_Arc"
umap.HC.red <- UMAPPlot(integrated.HC, cols = c("grey85", "grey75", "red2", "red4"), order = c("5P_YES", "1P_YES", "5P_NO", "1P_NO"), shuffle = F) +
  theme_classic(base_size = 5, base_family = "Arial") +
  scale_x_continuous(limits = c(-14,21), expand = c(0,0)) +
  scale_y_continuous(limits = c(-21,14), expand = c(0,0)) +
  theme(plot.background = element_rect(fill="transparent", colour=NA),
        panel.background = element_rect(fill="transparent", colour=NA),
        plot.title = element_blank(),
        plot.margin = margin(.5,.5,.5,.5, "char"),
        aspect.ratio = 1, axis.text = element_text(color = "black"),
        legend.position = c(.9,.9))
ggsave(umap.HC.red, filename = "/Users/arthur/Desktop/red_HC_UMAP.pdf", device = cairo_pdf, width = 15, height = 15, units = "cm")

Idents(integrated.Test) <- "pairings_Arc"
umap.Test.red <- UMAPPlot(integrated.Test, cols = c("grey85", "grey75", "red2", "red4"), order = c("5P_YES", "1P_YES", "5P_NO", "1P_NO"), shuffle = F) +
  theme_classic(base_size = 5, base_family = "Arial") +
  scale_x_continuous(limits = c(-14,21), expand = c(0,0)) +
  scale_y_continuous(limits = c(-21,14), expand = c(0,0)) +
  theme(plot.background = element_rect(fill="transparent", colour=NA),
        panel.background = element_rect(fill="transparent", colour=NA),
        plot.title = element_blank(),
        plot.margin = margin(.5,.5,.5,.5, "char"),
        aspect.ratio = 1, axis.text = element_text(color = "black"),
        legend.position = c(.9,.9))
ggsave(umap.Test.red, filename = "/Users/arthur/Desktop/red_Test_UMAP.pdf", device = cairo_pdf, width = 15, height = 15, units = "cm")


# Prop.graph red by HC/Test
library(ggradar)

temp2 <- integrated.HC@meta.data %>% group_by(pairings_Arc, comb.clusters) %>% summarize(n = n()) %>%
  pivot_wider(names_from = pairings_Arc, values_from = n) %>% replace(is.na(.), 0)
clu2 <- temp2$comb.clusters ; temp2$comb.clusters = NULL 
temp2 <- temp2 %>% as.data.frame() %>% apply(., 2, function(x){100*x/sum(x)}) %>% as.data.frame() %>% mutate(Cell_type = clu2)
radar.red.HC <- temp2 %>%
  pivot_longer(cols=!Cell_type, names_to = "group", values_to = "Percentage") %>%
  #separate(Group, into = c("sorted", "pairings"), sep = "_") %>%
  mutate(Cell_type = factor(Cell_type, levels = c("D1-MSNs", "D2-MSNs", "Interneurons", "Astrocytes", "Oligodendrocytes"))) %>%
  pivot_wider(names_from = Cell_type, values_from = Percentage) %>%
  dplyr::select(c("group", rev(names(.)))) %>%
  ggradar(grid.min = 0, grid.mid = 35, grid.max = 70,
          base.size = 5, font.radar = "Arial",
          legend.text.size = 5,
          values.radar = c("0%", "35%", "70%"),
          group.colours = c("grey85","red2", "grey75", "red4"),
          group.line.width = 1, group.point.size = 2,
          background.circle.colour = "transparent",
          gridline.min.colour = "grey", gridline.mid.colour = "grey", gridline.max.colour = "grey")
ggsave(radar.red.HC, filename = "/Users/arthur/Desktop/red_HC_radar.pdf", device = cairo_pdf, width = 15, height = 15, units = "cm") 

temp2 <- integrated.Test@meta.data %>% group_by(pairings_Arc, comb.clusters) %>% summarize(n = n()) %>%
  pivot_wider(names_from = pairings_Arc, values_from = n) %>% replace(is.na(.), 0)
clu2 <- temp2$comb.clusters ; temp2$comb.clusters = NULL 
temp2 <- temp2 %>% as.data.frame() %>% apply(., 2, function(x){100*x/sum(x)}) %>% as.data.frame() %>% mutate(Cell_type = clu2)
radar.red.Test <- temp2 %>%
  pivot_longer(cols=!Cell_type, names_to = "group", values_to = "Percentage") %>%
  #separate(Group, into = c("sorted", "pairings"), sep = "_") %>%
  mutate(Cell_type = factor(Cell_type, levels = c("D1-MSNs", "D2-MSNs", "Interneurons", "Astrocytes", "Oligodendrocytes"))) %>%
  pivot_wider(names_from = Cell_type, values_from = Percentage) %>%
  dplyr::select(c("group", rev(names(.)))) %>%
  ggradar(grid.min = 0, grid.mid = 35, grid.max = 70,
          base.size = 5, font.radar = "Arial",
          legend.text.size = 5,
          values.radar = c("0%", "35%", "70%"),
          group.colours = c("grey85","red2", "grey75", "red4"),
          group.line.width = 1, group.point.size = 2,
          background.circle.colour = "transparent",
          gridline.min.colour = "grey", gridline.mid.colour = "grey", gridline.max.colour = "grey")
ggsave(radar.red.Test, filename = "/Users/arthur/Desktop/red_Test_radar.pdf", device = cairo_pdf, width = 15, height = 15, units = "cm") 


# DEGs using all reds cells, correcting for pairings & cpp
allDEGs_YESvNO_split_HC <- function(integrated, cluster = "D1-MSNs", lfcT = .25, minpct = 0.25, test = "LR", p.th = .05, outwd = "/Users/arthur/Desktop/") {
  print(paste0("Working on ", cluster))
  comp <- list()
  Idents(integrated) <- "cluster_pairings_cpp_Arc"
  comp[["1P_HC_YESvNO"]] <- FindMarkers(integrated, ident.1 = paste0(cluster, "_1P_HC_YES"), ident.2 = paste0(cluster, "_1P_HC_NO"), assay = "SCT",  min.pct = minpct, logfc.threshold = lfcT, latent.vars = c("sorted"), test.use = test, verbose = FALSE) %>% as.data.frame() %>% mutate(gene = rownames(.)) %>%
    mutate(p_val_adj = p.adjust(.$p_val, "fdr"))
  comp[["5P_HC_YESvNO"]] <- FindMarkers(integrated, ident.1 = paste0(cluster, "_5P_HC_YES"), ident.2 = paste0(cluster, "_5P_HC_NO"), assay = "SCT",  min.pct = minpct, logfc.threshold = lfcT, latent.vars = c("sorted"), test.use = test, verbose = FALSE) %>% as.data.frame() %>% mutate(gene = rownames(.)) %>%
    mutate(p_val_adj = p.adjust(.$p_val, "fdr"))
  comp[["5Pv1P_HC_NO"]] <- FindMarkers(integrated, ident.1 = paste0(cluster, "_5P_HC_NO"), ident.2 = paste0(cluster, "_1P_HC_NO"), assay = "SCT",  min.pct = minpct, logfc.threshold = lfcT, latent.vars = c("sorted"), test.use = test, verbose = FALSE) %>% as.data.frame() %>% mutate(gene = rownames(.)) %>%
    mutate(p_val_adj = p.adjust(.$p_val, "fdr"))
  comp[["5Pv1P_HC_YES"]] <- FindMarkers(integrated, ident.1 = paste0(cluster, "_5P_HC_YES"), ident.2 = paste0(cluster, "_1P_HC_YES"), assay = "SCT",  min.pct = minpct, logfc.threshold = lfcT, latent.vars = c("sorted"), test.use = test, verbose = FALSE) %>% as.data.frame() %>% mutate(gene = rownames(.)) %>%
    mutate(p_val_adj = p.adjust(.$p_val, "fdr"))
  comp[["5Pv1P_HC_YESmNO"]] <- full_join(comp[["5Pv1P_HC_NO"]], comp[["5Pv1P_HC_YES"]], by = "gene", suffix = c(".no", "")) %>%
    dplyr::filter((p_val_adj.no > p.th | is.na(p_val_adj.no)==T) & p_val_adj < p.th) %>%
    dplyr::select(!matches("no"))
  
  wb <- createWorkbook()                                                    
  sheetnames <- names(comp)
  sheets <- lapply(sheetnames, createSheet, wb = wb)
  void <- Map(addDataFrame, comp, sheets)
  saveWorkbook(wb, file = paste0(outwd, cluster, "_YESvNO_HC_lfcT", lfcT, "_minpct", minpct, "_", test, ".xlsx"))
  
  d1 <- full_join(comp[["1P_HC_YESvNO"]], comp[["5P_HC_YESvNO"]], by = "gene", suffix = c(".yvn.HC.1P", ".yvn.HC.5P")) %>%
    full_join(setNames(comp[["5Pv1P_HC_YESmNO"]], paste0(names(comp[["5Pv1P_HC_YESmNO"]]), ".y.HC.5Pv1P")), by = join_by("gene"=="gene.y.HC.5Pv1P")) %>%
    dplyr::select(gene, avg_log2FC.yvn.HC.1P, p_val_adj.yvn.HC.1P, avg_log2FC.yvn.HC.5P, p_val_adj.yvn.HC.5P, avg_log2FC.y.HC.5Pv1P, p_val_adj.y.HC.5Pv1P 
                  
                  ) %>%
    setNames(c("gene", "lfc_yvn.HC.1P", "p_yvn.HC.1P", "lfc_yvn.HC.5P", "p_yvn.HC.5P", "lfc_y.HC.5Pv1P", "p_y.HC.5Pv1P"
               
               )) %>%
    mutate(cluster = cluster)
  
  if (nrow(dplyr::filter(comp[["1P_HC_YESvNO"]], p_val_adj < p.th))>0 ) { 
    DEGcount.yvn.HC.1P <- mutate(dplyr::filter(comp[["1P_HC_YESvNO"]], p_val_adj < p.th), updown = as.character(ifelse(avg_log2FC > 0, "up", "down"))) %>% group_by(updown) %>% summarise(count = n()) %>% mutate(cluster = cluster, pairwise = "1P_HC_YESvNO") } else {
    DEGcount.yvn.HC.1P <- data.frame(updown = c("up", "down"), count = c(0,0)) %>% mutate(cluster = cluster, pairwise = "1P_HC_YESvNO")
    }
  if (nrow(dplyr::filter(comp[["5P_HC_YESvNO"]], p_val_adj < p.th))>0 ) { 
    DEGcount.yvn.HC.5P <- mutate(dplyr::filter(comp[["5P_HC_YESvNO"]], p_val_adj < p.th), updown = as.character(ifelse(avg_log2FC > 0, "up", "down"))) %>% group_by(updown) %>% summarise(count = n()) %>% mutate(cluster = cluster, pairwise = "5P_HC_YESvNO") } else {
      DEGcount.yvn.HC.5P <- data.frame(updown = c("up", "down"), count = c(0,0)) %>% mutate(cluster = cluster, pairwise = "5P_HC_YESvNO")
    }
  if (nrow(dplyr::filter(comp[["5Pv1P_HC_YESmNO"]], p_val_adj < p.th))>0 ) { 
    DEGcount.y.HC.5Pv1P <- mutate(dplyr::filter(comp[["5Pv1P_HC_YESmNO"]], p_val_adj < p.th), updown = as.character(ifelse(avg_log2FC > 0, "up", "down"))) %>% group_by(updown) %>% summarise(count = n()) %>% mutate(cluster = cluster, pairwise = "5Pv1P_HC_YESmNO") } else {
      DEGcount.y.HC.5Pv1P <- data.frame(updown = c("up", "down"), count = c(0,0)) %>% mutate(cluster = cluster, pairwise = "5Pv1P_HC_YESmNO")
    }
  
  DEGcounts <- full_join(DEGcount.yvn.HC.1P, DEGcount.yvn.HC.5P) %>% full_join(DEGcount.y.HC.5Pv1P)
  
  output <- list(exported = comp, DEGs = d1, counts = DEGcounts)
  return(output)
}
allDEGs_YESvNO_split_Test <- function(integrated, cluster = "D1-MSNs", lfcT = .25, minpct = 0.25, test = "LR", p.th = .05, outwd = "/Users/arthur/Desktop/") {
  print(paste0("Working on ", cluster))
  comp <- list()
  Idents(integrated) <- "cluster_pairings_cpp_Arc"
  comp[["1P_Test_YESvNO"]] <- FindMarkers(integrated, ident.1 = paste0(cluster, "_1P_Test_YES"), ident.2 = paste0(cluster, "_1P_Test_NO"), assay = "SCT",  min.pct = minpct, logfc.threshold = lfcT, latent.vars = c("sorted"), test.use = test, verbose = FALSE) %>% as.data.frame() %>% mutate(gene = rownames(.)) %>%
    mutate(p_val_adj = p.adjust(.$p_val, "fdr"))
  comp[["5P_Test_YESvNO"]] <- FindMarkers(integrated, ident.1 = paste0(cluster, "_5P_Test_YES"), ident.2 = paste0(cluster, "_5P_Test_NO"), assay = "SCT",  min.pct = minpct, logfc.threshold = lfcT, latent.vars = c("sorted"), test.use = test, verbose = FALSE) %>% as.data.frame() %>% mutate(gene = rownames(.)) %>%
    mutate(p_val_adj = p.adjust(.$p_val, "fdr"))
  comp[["5Pv1P_Test_NO"]] <- FindMarkers(integrated, ident.1 = paste0(cluster, "_5P_Test_NO"), ident.2 = paste0(cluster, "_1P_Test_NO"), assay = "SCT",  min.pct = minpct, logfc.threshold = lfcT, latent.vars = c("sorted"), test.use = test, verbose = FALSE) %>% as.data.frame() %>% mutate(gene = rownames(.)) %>%
    mutate(p_val_adj = p.adjust(.$p_val, "fdr"))
  comp[["5Pv1P_Test_YES"]] <- FindMarkers(integrated, ident.1 = paste0(cluster, "_5P_Test_YES"), ident.2 = paste0(cluster, "_1P_Test_YES"), assay = "SCT",  min.pct = minpct, logfc.threshold = lfcT, latent.vars = c("sorted"), test.use = test, verbose = FALSE) %>% as.data.frame() %>% mutate(gene = rownames(.)) %>%
    mutate(p_val_adj = p.adjust(.$p_val, "fdr"))
  comp[["5Pv1P_Test_YESmNO"]] <- full_join(comp[["5Pv1P_Test_NO"]], comp[["5Pv1P_Test_YES"]], by = "gene", suffix = c(".no", "")) %>%
    dplyr::filter((p_val_adj.no > p.th | is.na(p_val_adj.no)==T) & p_val_adj < p.th) %>%
    dplyr::select(!matches("no"))
  
  wb <- createWorkbook()                                                    
  sheetnames <- names(comp)
  sheets <- lapply(sheetnames, createSheet, wb = wb)
  void <- Map(addDataFrame, comp, sheets)
  saveWorkbook(wb, file = paste0(outwd, cluster, "_YESvNO_Test_lfcT", lfcT, "_minpct", minpct, "_", test, ".xlsx"))
  
  d1 <- full_join(comp[["1P_Test_YESvNO"]], comp[["5P_Test_YESvNO"]], by = "gene", suffix = c(".yvn.Test.1P", ".yvn.Test.5P")) %>%
    full_join(setNames(comp[["5Pv1P_Test_YESmNO"]], paste0(names(comp[["5Pv1P_Test_YESmNO"]]), ".y.Test.5Pv1P")), by = join_by("gene"=="gene.y.Test.5Pv1P")) %>%
    dplyr::select(gene, 
                  avg_log2FC.yvn.Test.1P, p_val_adj.yvn.Test.1P, avg_log2FC.yvn.Test.5P, p_val_adj.yvn.Test.5P, avg_log2FC.y.Test.5Pv1P, p_val_adj.y.Test.5Pv1P 
    ) %>%
    setNames(c("gene", 
               "lfc_yvn.Test.1P", "p_yvn.Test.1P", "lfc_yvn.Test.5P", "p_yvn.Test.5P", "lfc_y.Test.5Pv1P", "p_y.Test.5Pv1P" 
    )) %>%
    mutate(cluster = cluster)
  
  if (nrow(dplyr::filter(comp[["1P_Test_YESvNO"]], p_val_adj < p.th))>0 ) { 
    DEGcount.yvn.Test.1P <- mutate(dplyr::filter(comp[["1P_Test_YESvNO"]], p_val_adj < p.th), updown = as.character(ifelse(avg_log2FC > 0, "up", "down"))) %>% group_by(updown) %>% summarise(count = n()) %>% mutate(cluster = cluster, pairwise = "1P_Test_YESvNO") } else {
      DEGcount.yvn.Test.1P <- data.frame(updown = c("up", "down"), count = c(0,0)) %>% mutate(cluster = cluster, pairwise = "1P_Test_YESvNO")
    }
  if (nrow(dplyr::filter(comp[["5P_Test_YESvNO"]], p_val_adj < p.th))>0 ) { 
    DEGcount.yvn.Test.5P <- mutate(dplyr::filter(comp[["5P_Test_YESvNO"]], p_val_adj < p.th), updown = as.character(ifelse(avg_log2FC > 0, "up", "down"))) %>% group_by(updown) %>% summarise(count = n()) %>% mutate(cluster = cluster, pairwise = "5P_Test_YESvNO") } else {
      DEGcount.yvn.Test.5P <- data.frame(updown = c("up", "down"), count = c(0,0)) %>% mutate(cluster = cluster, pairwise = "5P_Test_YESvNO")
    }
  if (nrow(dplyr::filter(comp[["5Pv1P_Test_YESmNO"]], p_val_adj < p.th))>0 ) { 
    DEGcount.y.Test.5Pv1P <- mutate(dplyr::filter(comp[["5Pv1P_Test_YESmNO"]], p_val_adj < p.th), updown = as.character(ifelse(avg_log2FC > 0, "up", "down"))) %>% group_by(updown) %>% summarise(count = n()) %>% mutate(cluster = cluster, pairwise = "5Pv1P_Test_YESmNO") } else {
      DEGcount.y.Test.5Pv1P <- data.frame(updown = c("up", "down"), count = c(0,0)) %>% mutate(cluster = cluster, pairwise = "5Pv1P_Test_YESmNO")
    }
  
  DEGcounts <- full_join(DEGcount.yvn.Test.1P, DEGcount.yvn.Test.5P) %>% full_join(DEGcount.y.Test.5Pv1P)
  
  output <- list(exported = comp, DEGs = d1, counts = DEGcounts)
  return(output)
}
allDEGs_YESvNO_together <- function(integrated, cluster = "D1-MSNs", lfcT = .25, minpct = 0.25, test = "LR", p.th = .05, outwd = "/Users/arthur/Desktop/") {
  print(paste0("Working on ", cluster))
  comp <- list()
  Idents(integrated) <- "cluster_Arc"
  comp[["all_YESvNO"]] <- FindMarkers(integrated, ident.1 = paste0(cluster, "_YES"), ident.2 = paste0(cluster, "_NO"), assay = "SCT",  min.pct = minpct, logfc.threshold = lfcT, latent.vars = c("sorted", "pairings", "cpp"), test.use = test, verbose = FALSE) %>% as.data.frame() %>% mutate(gene = rownames(.)) %>%
    mutate(p_val_adj = p.adjust(.$p_val, "fdr"))
  
  wb <- createWorkbook()                                                    
  sheetnames <- names(comp)
  sheets <- lapply(sheetnames, createSheet, wb = wb)
  void <- Map(addDataFrame, comp, sheets)
  saveWorkbook(wb, file = paste0(outwd, cluster, "_YESvNO_lfcT", lfcT, "_minpct", minpct, "_", test, ".xlsx"))
  
  d1 <- setNames(comp[["all_YESvNO"]], paste0(names(comp[["all_YESvNO"]]), ".yvn.all")) %>%
    dplyr::select( 
                  avg_log2FC.yvn.all, p_val_adj.yvn.all) %>%
    setNames(c("gene", 
               "lfc_yvn.all", "p_yvn.all")) %>%
    mutate(cluster = cluster)
  
  
  if (nrow(dplyr::filter(comp[["all_YESvNO"]], p_val_adj < p.th))>0 ) { 
    DEGcount.yvn.all <- mutate(dplyr::filter(comp[["all_YESvNO"]], p_val_adj < p.th), updown = as.character(ifelse(avg_log2FC > 0, "up", "down"))) %>% group_by(updown) %>% summarise(count = n()) %>% mutate(cluster = cluster, pairwise = "all_YESvNO") } else {
      DEGcount.yvn.all <- data.frame(updown = c("up", "down"), count = c(0,0)) %>% mutate(cluster = cluster, pairwise = "all_YESvNO")
    }
  DEGcounts <- DEGcount.yvn.all
  
  output <- list(exported = comp, DEGs = d1, counts = DEGcounts)
  return(output)
}

# For all clusters
compsYvN.HC.LR <- list()
compsYvN.Test.LR <- list()
for (i in c("D1-MSNs", "D2-MSNs")) {
  compsYvN.HC.LR[[i]] <- allDEGs_YESvNO_split_HC(integrated, cluster = i, lfcT = log2(1.15), minpct = .3, test = "LR", p.th = .1)
  compsYvN.Test.LR[[i]] <- allDEGs_YESvNO_split_Test(integrated, cluster = i, lfcT = log2(1.15), minpct = .3, test = "LR", p.th = .1)
}

compsYvN.HC.LR.nolfcT <- list()
compsYvN.Test.LR.nolfcT <- list()
compsYvN.HC.LR.nolfcT[["D1-MSNs"]] <- allDEGs_YESvNO_split_HC(integrated, cluster = "D1-MSNs", lfcT = -Inf, minpct = 0.01, test = "LR", p.th = .1)
compsYvN.Test.LR.nolfcT[["D1-MSNs"]] <- allDEGs_YESvNO_split_Test(integrated, cluster = "D1-MSNs", lfcT = -Inf, minpct = 0.01, test = "LR", p.th = .1)
compsYvN.HC.LR.nolfcT[["D2-MSNs"]] <- allDEGs_YESvNO_split_HC(integrated, cluster = "D2-MSNs", lfcT = -Inf, minpct = 0.01, test = "LR", p.th = .1)
compsYvN.Test.LR.nolfcT[["D2-MSNs"]] <- allDEGs_YESvNO_split_Test(integrated, cluster = "D2-MSNs", lfcT = -Inf, minpct = 0.01, test = "LR", p.th = .1)


## to use, TEst only comparisons with LR testing
final <- compsYvN.Test.LR 
counts <- lapply(final, pluck, "counts") %>% data.table::rbindlist(idcol = "Cluster") %>%
  mutate(updown = factor(updown, levels = c("up", "down")), 
         pairwise = factor(pairwise, levels = c("1P_Test_YESvNO", "1P_Test_YESvNO", "5Pv1P_Test_YESmNO")),
         Cluster = factor(Cluster, levels = c("D1-MSNs", "D2-MSNs")))

p.counts <- counts %>% 
  ggplot(aes(x = Cluster, y = count)) + 
  facet_wrap(~pairwise, ncol = 1) +
  coord_flip() +
  geom_col(aes(fill = updown)) + 
  scale_x_discrete(limits = rev) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("yellow", "blue")) +
  guides(fill = guide_legend(title = "")) +
  ylab("Number of DEGs") +
  theme_classic(base_size = 5, base_family = "Arial") +
  theme(axis.text = element_text(color = "black"), axis.ticks = element_line(color = "black"), axis.title.y = element_blank(),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = c(0.9, 0.9), legend.background = element_rect(fill="transparent"),
        strip.background = element_blank(), strip.placement = "inside", strip.text = element_text(hjust = 0))

# final figure
plot_grid(umap.Test.red, NULL, plot_grid(p.counts, NULL, ncol = 1, align = "hv", axis = "tblr"), plot_grid(d1.plots[["Reln"]], d1.plots[["Kcnip4"]], d1.plots[["Prkca"]],d1.plots[["Pde7b"]], nrow=2, align = "hv", axis = "tblr"), nrow = 2, align = "hv", axis = "trbl") %>%
  ggsave(filename = "/Users/arthur/Desktop/green.pdf", device = cairo_pdf, width = 21.6, height = 27.9, units = "cm") 


# Heatmaps

plot.heatmap <- function(df, ttl) {
  df %>%
    dplyr::select(gene, lfc_yvn.Test.1P, lfc_yvn.Test.5P, lfc_y.Test.5Pv1P) %>%
    pivot_longer(cols = !gene, names_to = "comp", values_to = "log2FC") %>%
    mutate(comp = factor(case_when(comp == "lfc_yvn.Test.1P" ~ "1PTest - YESvNO", comp == "lfc_yvn.Test.5P" ~ "5PTest - YESvNO", comp == "lfc_y.Test.5Pv1P" ~ "YES - 5Pv1P"))) %>%
    ggplot(aes(x= gene, y = comp, fill = log2FC)) +
    geom_tile() +
    ggtitle(ttl) +
    scale_y_discrete(limits=rev) +
    scale_fill_gradientn(colours = c("blue","blue","black","yellow", "yellow"), na.value = "grey80",
                         values = scales::rescale(c(-1.52, -.8, 0, .8, 1.52), from = c(-1.52, 1.52)), 
                         limits = c(-1.52, 1.52),
                         breaks = c(-0.8, 0, 0.8),
                         guide = guide_colorbar(title  = "log2FC", title.vjust = 0.9, frame.colour = "black", frame.linewidth = .2*2.8, ticks.linewidth = .2*2.8, ticks.colour = "black", draw.ulim = F, draw.llim = F, title.position = "top", direction = "vertical", barheight = unit(2, "char"), barwidth = unit(.5, "char"))) +
    theme_classic(base_size = 6, base_family = "Arial") +
    theme(plot.background = element_rect(fill="transparent", colour=NA),
          panel.background = element_rect(fill="transparent", colour=NA),
          aspect.ratio = 1/6,
          plot.margin = margin(.5,.5,.5,.5, "line"),
          plot.title = element_text(size = 6, hjust = .5),
          axis.line = element_blank(), 
          axis.title.y = element_blank(),  axis.title.x = element_blank(), 
          axis.ticks = element_blank(),
          axis.text = element_text(size = 5, colour = "black"), axis.text.x = element_blank(),
          strip.text = element_text(size = 5, color = "black"), strip.background = element_blank(),
          legend.position = "right", legend.justification = c("center", "center"), legend.background = element_rect(fill="transparent"), legend.box.background = element_rect(fill="transparent", color="transparent"),
          legend.text = element_text(size = 5), legend.margin = margin(0,0,0,0), legend.box.margin = margin(0,0,0,0), legend.spacing = unit(-.1, "char"), legend.title = element_text(size=5)) }

hm.d1.1P <- data.frame(gene = dplyr::filter(compsYvN.Test.LR$`D1-MSNs`$DEGs, p_yvn.Test.1P < .1)$gene) %>%
  left_join(compsYvN.Test.LR.nolfcT$`D1-MSNs`$DEGs) %>%
  mutate(gene = factor(gene, levels = arrange(., lfc_yvn.Test.1P)$gene)) %>%
  plot.heatmap(ttl = "D1 - YESvNO_1PTest") # + theme(axis.text.x = element_text(size = 5, colour = "black", angle = 90, hjust=0.95,vjust=0.2))
hm.d1.5P <- data.frame(gene = dplyr::filter(compsYvN.Test.LR$`D1-MSNs`$DEGs, p_yvn.Test.5P < .1)$gene) %>%
  left_join(compsYvN.Test.LR.nolfcT$`D1-MSNs`$DEGs) %>%
  mutate(gene = factor(gene, levels = arrange(., lfc_yvn.Test.5P)$gene)) %>%
  plot.heatmap(ttl = "D1 - YESvNO_5PTest")
hm.d1.5Pv1P <- dplyr::filter(compsYvN.Test.LR$`D1-MSNs`$DEGs, p_y.Test.5Pv1P < .1) %>% dplyr::select(gene, lfc_y.Test.5Pv1P, p_y.Test.5Pv1P, cluster) %>%
  left_join(dplyr::select(compsYvN.Test.LR.nolfcT$`D1-MSNs`$DEGs, gene, lfc_yvn.Test.1P, p_yvn.Test.1P, lfc_yvn.Test.5P, p_yvn.Test.5P)) %>%
  mutate(gene = factor(gene, levels = arrange(., lfc_y.Test.5Pv1P)$gene)) %>%
  plot.heatmap(ttl = "D1 - YES_5Pv1PTest")

hm.D2.1P <- data.frame(gene = dplyr::filter(compsYvN.Test.LR$`D2-MSNs`$DEGs, p_yvn.Test.1P < .1)$gene) %>%
  left_join(compsYvN.Test.LR.nolfcT$`D2-MSNs`$DEGs) %>%
  mutate(gene = factor(gene, levels = arrange(., lfc_yvn.Test.1P)$gene)) %>%
  plot.heatmap(ttl = "D2 - YESvNO_1PTest") # + theme(axis.text.x = element_text(size = 5, colour = "black", angle = 90, hjust=0.95,vjust=0.2))
hm.D2.5P <- data.frame(gene = dplyr::filter(compsYvN.Test.LR$`D2-MSNs`$DEGs, p_yvn.Test.5P < .1)$gene) %>%
  left_join(compsYvN.Test.LR.nolfcT$`D2-MSNs`$DEGs) %>%
  mutate(gene = factor(gene, levels = arrange(., lfc_yvn.Test.5P)$gene)) %>%
  plot.heatmap(ttl = "D2 - YESvNO_5PTest")
hm.D2.5Pv1P <- dplyr::filter(compsYvN.Test.LR$`D2-MSNs`$DEGs, p_y.Test.5Pv1P < .1) %>% dplyr::select(gene, lfc_y.Test.5Pv1P, p_y.Test.5Pv1P, cluster) %>%
  left_join(dplyr::select(compsYvN.Test.LR.nolfcT$`D2-MSNs`$DEGs, gene, lfc_yvn.Test.1P, p_yvn.Test.1P, lfc_yvn.Test.5P, p_yvn.Test.5P)) %>%
  mutate(gene = factor(gene, levels = arrange(., lfc_y.Test.5Pv1P)$gene)) %>%
  plot.heatmap(ttl = "D2 - YES_5Pv1PTest")

all.heatmaps <- list(hm.d1.1P, hm.d1.5P, hm.d1.5Pv1P, hm.D2.1P, hm.D2.5P, hm.D2.5Pv1P) ; rm(hm.d1.1P, hm.d1.5P, hm.d1.5Pv1P, hm.D2.1P, hm.D2.5P, hm.D2.5Pv1P)
maxgenes <- lapply(all.heatmaps, function(x) {nrow(x$data)/3}) %>% unlist() %>% max()

addextra <- function(t, maxadd) {
  toadd <- maxadd - (nrow(t$data)/3)
  if (toadd > 0) {
    tobind <- data.frame(gene = sort(rep(paste0("extra", seq(1,toadd)), 3)), comp = rep(c("1PTest - YESvNO", "5PTest - YESvNO", "YES - 5Pv1P"), toadd), log2FC = rep(NA, toadd*3))
    t$data <- t$data %>% rbind(tobind)
  }
  return(t)
}

all.heatmaps <- lapply(all.heatmaps, addextra, maxadd = 383)

plot_grid(plotlist = all.heatmaps, ncol = 2) %>%
  ggsave(filename = "/Users/arthur/Desktop/red_heatmaps.pdf", device = cairo_pdf, width = 21.6, height = 27.9, units = "cm") 

## Correlations

cor.list <- list()
Idents(integrated.Test) <- "cluster_pairings_cpp_Arc"
for (c in c("D1-MSNs", "D2-MSNs")) {
  for (d in c("1P", "5P")) {
    print(c)
    print(d)
    act <- subset(integrated.Test, idents = paste0(c, "_", d,  "_Test_YES"))
    matrix_mod <- as.matrix(act@assays$SCT@data)
    arc <- as.numeric(matrix_mod["Arc",])
    matrix_mod <- matrix_mod[rownames(matrix_mod)!="Arc",]
    matrix_mod <- matrix_mod[rowSums(matrix_mod)>0,]      # filter out non expressed genes for speed (log(SCTdata) = 0 for all cells)
    sorted <- act@meta.data$sorted
    cor2factor <- function(y, x, f1) {
      model <- summary(lm(y ~ x + f1))
      out1 <- model$coefficients %>% as.data.frame() %>% dplyr::select(Estimate, `Pr(>|t|)`) %>% 
        setNames(c("Estimate", "p.value")) %>% mutate(factor = c("Intercept", "x", "f1")) %>%
        pivot_wider(names_from = factor, values_from = c(Estimate, p.value)) %>%
        dplyr::select(-matches("Intercept"))
    }
    if (ncol(matrix_mod) >= 5) {
      correlations <- apply(matrix_mod,1,cor2factor, x = arc, f1 = sorted) %>%
        data.table::rbindlist(idcol = "gene") %>%
        arrange(Estimate_x) %>%
        dplyr::select(gene, Estimate_x, p.value_x)
      totest <- correlations %>% dplyr::filter(abs(Estimate_x) > .25) %>%
        mutate(p_val_adj_x = p.adjust(.$p.value_x, "fdr"))
      correlations <- left_join(correlations, totest)  
    } else {
      correlations = data.frame(gene = NA, Estimate_x = NA, p.value_x =NA, p_val_adj_x = NA) }
    write.csv(correlations, file = paste("/Users/arthur/Desktop/correlations_in_POS_", c, "_", d, ".csv", sep = ""), row.names = F)
    cor.list[[paste(c, d, sep = "_")]] <- correlations
  }
}

plot.heatmap.3 <- function(df, ttl) {
  df %>%
    dplyr::select(gene, Estimate_x) %>%
    pivot_longer(cols = !gene, names_to = "comp", values_to = "Estimate") %>%
    mutate(comp = factor(case_when(comp == "Estimate_x" ~ "vs Arc expression"))) %>%
    ggplot(aes(x= gene, y = comp, fill = Estimate)) +
    geom_tile() +
    ggtitle(ttl) +
    scale_y_discrete(limits=rev) +
    scale_fill_gradientn(colours = c("deepskyblue","deepskyblue","black","orange", "orange"), na.value = "grey80",
                         values = scales::rescale(c(-1, -.5, 0, .5, 1), from = c(-1, 1)), 
                         limits = c(-1, 1),
                         breaks = c(-1, -0.5, 0, 0.5, 1),
                         guide = guide_colorbar(title  = "Cor.estimate", title.vjust = 0.9, frame.colour = "black", frame.linewidth = .2*2.8, ticks.linewidth = .2*2.8, ticks.colour = "black", draw.ulim = F, draw.llim = F, title.position = "top", direction = "vertical", barheight = unit(2, "char"), barwidth = unit(.5, "char"))) +
    theme_classic(base_size = 6, base_family = "Arial") +
    theme(plot.background = element_rect(fill="transparent", colour=NA),
          panel.background = element_rect(fill="transparent", colour=NA),
          aspect.ratio = 1/6,
          plot.margin = margin(.5,.5,.5,.5, "line"),
          plot.title = element_text(size = 6, hjust = .5),
          axis.line = element_blank(), 
          axis.title.y = element_blank(),  axis.title.x = element_blank(), 
          axis.ticks = element_blank(),
          axis.text = element_text(size = 5, colour = "black"), axis.text.x = element_blank(),
          strip.text = element_text(size = 5, color = "black"), strip.background = element_blank(),
          legend.position = "right", legend.justification = c("center", "center"), legend.background = element_rect(fill="transparent"), legend.box.background = element_rect(fill="transparent", color="transparent"),
          legend.text = element_text(size = 5), legend.margin = margin(0,0,0,0), legend.box.margin = margin(0,0,0,0), legend.spacing = unit(-.1, "char"), legend.title = element_text(size=5)) }
hm.d1.1P.corr <- filter(cor.list$`D1-MSNs_1P`, p_val_adj_x < .1) %>%
  mutate(gene = factor(gene, levels = arrange(., Estimate_x)$gene)) %>%
  plot.heatmap.3(ttl = "D1 - 1P - corr")
hm.d1.5P.corr <- filter(cor.list$`D1-MSNs_5P`, p_val_adj_x < .1) %>%
  mutate(gene = factor(gene, levels = arrange(., Estimate_x)$gene)) %>%
  plot.heatmap.3(ttl = "D1 - 5P - corr")
hm.d2.1P.corr <- filter(cor.list$`D2-MSNs_1P`, p_val_adj_x < .1) %>%
  mutate(gene = factor(gene, levels = arrange(., Estimate_x)$gene)) %>%
  plot.heatmap.3(ttl = "D2 - 1P - corr")
hm.d2.5P.corr <- filter(cor.list$`D2-MSNs_5P`, p_val_adj_x < .1) %>%
  mutate(gene = factor(gene, levels = arrange(., Estimate_x)$gene)) %>%
  plot.heatmap.3(ttl = "D2 - 5P - corr")

all.heatmaps <- list(hm.d1.1P.corr, hm.d1.5P.corr, hm.d2.1P.corr, hm.d2.5P.corr)
maxgenes <- lapply(all.heatmaps, function(x) {nrow(x$data)}) %>% unlist() %>% max()

addextra3 <- function(t, maxadd) {
  toadd <- maxadd - (nrow(t$data)/1)
  if (toadd > 0) {
    tobind <- data.frame(gene = sort(rep(paste0("extra", seq(1,toadd)), 1)), comp = rep(c("vs Arc expression"), toadd), Estimate = rep(NA, toadd*1))
    t$data <- t$data %>% rbind(tobind)
  }
  return(t)
}

all.heatmaps <- lapply(all.heatmaps, addextra3, maxadd = 383)

plot_grid(plotlist = all.heatmaps, ncol = 2) %>%
  ggsave(filename = "/Users/arthur/Desktop/red_corr_heatmaps.pdf", device = cairo_pdf, width = 21.6, height = 27.9, units = "cm") 

# Plot select IEG examples

all.degs.1P <- data.table::rbindlist(list(compsYvN.Test.LR[["D1-MSNs"]]$DEGs, compsYvN.Test.LR[["D2-MSNs"]]$DEGs)) %>%
  dplyr::select(cluster, gene, lfc_yvn.Test.1P, p_yvn.Test.1P) %>%
  setNames(c("cluster", "gene", "value", "padj")) %>%
  mutate(cluster = factor(cluster, levels = c("D1-MSNs", "D2-MSNs")))
all.degs.5P <- data.table::rbindlist(list(compsYvN.Test.LR.nolfcT[["D1-MSNs"]]$DEGs, compsYvN.Test.LR.nolfcT[["D2-MSNs"]]$DEGs)) %>%
  dplyr::select(cluster, gene, lfc_yvn.Test.5P) %>%
  left_join(dplyr::select(data.table::rbindlist(list(compsYvN.Test.LR[["D1-MSNs"]]$DEGs, compsYvN.Test.LR[["D2-MSNs"]]$DEGs)), 1,5,8)) %>%
  setNames(c("cluster", "gene", "value", "padj")) %>%
  mutate(cluster = factor(cluster, levels = c("D1-MSNs", "D2-MSNs")))
all.cors <- rbindlist(cor.list, idcol = "cluster") %>%
  dplyr::select(cluster, gene, Estimate_x, p_val_adj_x) %>%
  setNames(c("cluster", "gene", "value", "padj")) %>%
  rbind(data.frame(cluster = c("D1-MSNs_1P", "D1-MSNs_5P", "D2-MSNs_1P", "D2-MSNs_5P"), gene = c("Arc", "Arc", "Arc", "Arc"), value = rep(NA, 4), padj = rep(NA, 4))) %>%
  mutate(cluster = factor(cluster, levels = c("D1-MSNs_1P", "D1-MSNs_5P", "D2-MSNs_1P", "D2-MSNs_5P")))

iegs <- data.frame(gene = c("Arc", "Egr1", "Egr3", "Fos", "Fosl2", "Homer1", "Junb", "Npas4", "Nr4a1", "Nr4a3"))

plot.select <- function(df, sel, sc = c(-1, -.5, 0, .5, 1), sccols = c("deepskyblue","deepskyblue","black","orange", "orange")) {
  sel %>% left_join(df) %>% mutate(label = case_when(padj < .1 & padj > .05 ~ ".",
                                                     padj < .05 & padj > .01 ~ "*",
                                                     padj < .01 & padj > .001 ~ "**",
                                                     padj < .001 & padj > .0001 ~ "***",
                                                     padj < .0001 ~ "****",
                                                     T ~ "")) %>%
    ggplot(aes(y = gene, fill = value, x= cluster) ) +
    geom_tile() +
    geom_text(aes(label = label), colour = "black", size = 5/2.8) +
    scale_y_discrete(limits=rev, drop = F) +
    scale_fill_gradientn(colours = sccols, na.value = "grey80",
                         values = scales::rescale(sc, from = c(min(sc), max(sc))), 
                         limits = c(min(sc), max(sc)),
                         breaks = c(sc[2], 0, sc[4]),
                         guide = guide_colorbar(title  = "value", title.vjust = 0.9, frame.colour = "black", frame.linewidth = .2*2.8, ticks.linewidth = .2*2.8, ticks.colour = "black", draw.ulim = F, draw.llim = F, title.position = "top", direction = "vertical", barheight = unit(2, "char"), barwidth = unit(.5, "char"))) +
    theme_classic(base_size = 5, base_family = "Arial") +
    theme(plot.background = element_rect(fill="transparent", colour=NA),
          panel.background = element_rect(fill="transparent", colour=NA),
          aspect.ratio = 6,
          plot.margin = margin(.5,.5,.5,.5, "line"),
          plot.title = element_text(size = 6, hjust = .5),
          axis.line = element_blank(), 
          axis.title.y = element_blank(),  axis.title.x = element_blank(), 
          axis.ticks = element_blank(),
          axis.text.x = element_text(size = 5, angle = 45, vjust = 0.8, hjust = 1), axis.text.y = element_text(size = 5, face = "italic"),
          strip.text = element_text(size = 5, color = "black"), strip.background = element_blank(),
          legend.position = "right", legend.justification = c("center", "center"), legend.background = element_rect(fill="transparent"), legend.box.background = element_rect(fill="transparent", color="transparent"),
          legend.text = element_text(size = 5), legend.margin = margin(0,0,0,0), legend.box.margin = margin(0,0,0,0), legend.spacing = unit(-.1, "char"), legend.title = element_text(size=5)) 
}

p.ieg.degs.1P <- plot.select(all.degs.1P, iegs, sc = c(-1.52, -.8, 0, .8, 1.52), sccols = c("blue","blue","black","yellow", "yellow"))
p.ieg.degs.5P <- plot.select(all.degs.5P, iegs, sc = c(-1.52, -.8, 0, .8, 1.52), sccols = c("blue","blue","black","yellow", "yellow"))
p.ieg.cors <- plot.select(all.cors, iegs, sc = c(-1.2, -.5, 0, .5, 1.2), sccols = c("deepskyblue","deepskyblue","black","orange", "orange"))

# final figure
plot_grid(p.ieg.degs.1P, p.ieg.degs.5P, p.ieg.cors, ncol=3, align = "hv", axis = "tblr") %>%
  ggsave(filename = "/Users/arthur/Desktop/red_IEGs.pdf", device = cairo_pdf, width = 21.6, height = 27.9, units = "cm") 

plot_gene_Test <- function(integrated, gene = "Gria4", cluster = "D1-MSNs") {
  Idents(integrated) <- "cluster_cpp"
  d1cells <- subset(integrated, idents = paste0(cluster, "_Test"))
  p <- VlnPlot(d1cells, features = c(gene), split.by = "cluster_pairings_cpp_Arc", assay = "SCT", slot = "data", pt.size = 0) +
    scale_fill_manual(values = c("grey80","grey80","red1","red1")) +
    scale_y_continuous(expand = c(.1,.1)) +
    xlab(cluster) +
    ylab("Normalized gene expression") +
    theme_classic(base_size = 5, base_line_size = .2, base_family = "Arial") +
    theme(legend.position = "bottom", legend.key.size = unit(.1, "line"), legend.direction = "vertical",
          aspect.ratio = 1, 
          axis.text = element_text(color = "black"), axis.ticks = element_line(color="black"),
          axis.text.x = element_blank(), axis.ticks.x = element_blank())
  return(p)
}

gria4 <- plot_gene_Test(integrated, "Gria4", "D1-MSNs")
reln <- plot_gene_Test(integrated, "Reln", "D1-MSNs")
mef2c <- plot_gene_Test(integrated, "Mef2c", "D1-MSNs")
arfgef2 <- plot_gene_Test(integrated, "Arfgef2", "D1-MSNs")

# final figure
plot_grid(NULL, NULL, plot_grid(p.counts, NULL, ncol = 1, align = "hv", axis = "tblr"), plot_grid(gria4, reln, mef2c, arfgef2, nrow=2, align = "hv", axis = "tblr"), nrow = 2, align = "hv", axis = "trbl") %>%
  ggsave(filename = "/Users/arthur/Desktop/red2.pdf", device = cairo_pdf, width = 21.6, height = 27.9, units = "cm") 


# Dot plots example genes
Idents(integrated.Test) <- "comb.clusters"
tempD1 <- subset(integrated.Test, idents = "D1-MSNs")
ex.D1 <- DotPlot(tempD1, features = rev(c("Arhgef3", "Mef2c", "Npas2", "Kcnh5" )), 
              #cluster.idents = F,
              assay = "SCT", dot.scale = 6, group.by = "cluster_pairings_cpp_Arc",
              cols = c("lightgrey", "#F564E3"),
              dot.min = .01, col.min = 0, col.max = 1, scale.min = .1, scale.max = 100) + 
  coord_flip() +
  scale_y_discrete(limits=c("D1-MSNs_1P_Test_NO", "D1-MSNs_5P_Test_NO", "D1-MSNs_1P_Test_YES", "D1-MSNs_5P_Test_YES")) +
  theme_classic(base_family = "Arial", base_size = 5) +
  theme(plot.background = element_rect(fill="transparent", colour=NA),
        panel.background = element_rect(fill="transparent", colour=NA),
        #aspect.ratio = 6/4, 
        axis.title = element_blank(), axis.text = element_text(color = "black"), axis.ticks = element_line(color="black"),
        axis.text.x = element_text(face = "italic", angle = 45, vjust = 0.8, hjust = 1),
        legend.key.size = unit(.5, "line"), legend.position = "bottom", legend.box = "vertical")

tempD2 <- subset(integrated.Test, idents = "D2-MSNs")
ex.D2 <- DotPlot(tempD2, features = rev(c("Arhgef3", "Tet2", "Cacng3", "Hcn1")), 
                 #cluster.idents = F,
                 assay = "SCT", dot.scale = 6, group.by = "cluster_pairings_cpp_Arc",
                 cols = c("lightgrey", "#619CFF"),
                 dot.min = .01, col.min = 0, col.max = 1, scale.min = .1, scale.max = 100) + 
  coord_flip() +
  scale_y_discrete(limits=c("D2-MSNs_1P_Test_NO", "D2-MSNs_5P_Test_NO", "D2-MSNs_1P_Test_YES", "D2-MSNs_5P_Test_YES")) +
  theme_classic(base_family = "Arial", base_size = 5) +
  theme(plot.background = element_rect(fill="transparent", colour=NA),
        panel.background = element_rect(fill="transparent", colour=NA),
        #aspect.ratio = 13/4, 
        axis.title = element_blank(), axis.text = element_text(color = "black"), axis.ticks = element_line(color="black"),
        axis.text.x = element_text(face = "italic", angle = 45, vjust = 0.8, hjust = 1),
        legend.key.size = unit(.5, "line"), legend.position = "bottom", legend.box = "vertical")

ex<-plot_grid(ex.D1, ex.D2, align = "hv", axis = "tblr", ncol=2)
ggsave(ex, filename = "/Users/arthur/Desktop/ex_red_DotPlot.pdf", device = cairo_pdf, width = 9, height = 9, units = "cm") 


# IPA upstream
setwd("/Users/arthur/Desktop/IPA_redTest/")
allUpsReg <- lapply(list.files(pattern=".csv"), read_csv)
names(allUpsReg) <- list.files(pattern=".csv")
df.UpsReg <- rbindlist(allUpsReg, idcol = "File", fill=T) %>%
  mutate(comp = gsub("_UpsReg.csv", "", File), File=NULL) %>%
  dplyr::select(comp, 1, 3,5, 6, 8,9, 10) %>%
  setNames(c("comp", "UpsReg", "type", "biascorr.zscore", "zscore", "pval", "pval_adj", "targets")) %>%
  mutate(ngenes = 1+str_count(targets, ",")) %>%
  pivot_wider(id_cols = c("UpsReg", "type"), names_sep = "__", names_from = "comp", values_from = c("biascorr.zscore", "zscore", "pval", "pval_adj", "targets", "ngenes")) %>%
  mutate(n.NA = rowSums(is.na(.[,15:20]))) %>%
  mutate(ngenes__D1_1P_Test_YESvNO = ifelse(is.na(ngenes__D1_1P_Test_YESvNO)==T, 0, ngenes__D1_1P_Test_YESvNO), 
         ngenes__D1_5P_Test_YESvNO = ifelse(is.na(ngenes__D1_5P_Test_YESvNO)==T, 0, ngenes__D1_5P_Test_YESvNO),
         ngenes__D1_5Pv1P_Test_YES = ifelse(is.na(ngenes__D1_5Pv1P_Test_YES)==T, 0, ngenes__D1_5Pv1P_Test_YES),
         ngenes__D2_1P_Test_YESvNO = ifelse(is.na(ngenes__D2_1P_Test_YESvNO)==T, 0, ngenes__D2_1P_Test_YESvNO), 
         ngenes__D2_5P_Test_YESvNO = ifelse(is.na(ngenes__D2_5P_Test_YESvNO)==T, 0, ngenes__D2_5P_Test_YESvNO),
         ngenes__D2_5Pv1P_Test_YES = ifelse(is.na(ngenes__D2_5Pv1P_Test_YES)==T, 0, ngenes__D2_5Pv1P_Test_YES)) %>%
  rowwise() %>%
  mutate(ordersig = sum(pval_adj__D1_1P_Test_YESvNO, pval_adj__D1_5P_Test_YESvNO, pval_adj__D1_5Pv1P_Test_YES, pval_adj__D2_1P_Test_YESvNO, pval_adj__D2_5P_Test_YESvNO, pval_adj__D2_5Pv1P_Test_YES, na.rm =T),
         avg.ngenes = sum(ngenes__D1_1P_Test_YESvNO, ngenes__D1_5P_Test_YESvNO, ngenes__D1_5Pv1P_Test_YES, ngenes__D2_1P_Test_YESvNO, ngenes__D2_5P_Test_YESvNO, ngenes__D2_5Pv1P_Test_YES, na.rm =T)/6) %>%
  arrange(n.NA, ordersig)%>%
  dplyr::filter(avg.ngenes >= 3) %>%
  dplyr::select(1:8, 21:26) %>%
  ungroup()

chems <- df.UpsReg %>%
  dplyr::filter(grepl("chemical|drug", type)==T) %>%
  dplyr::filter(grepl("toxicant|reagent|non", type)==F) %>%
  head(12) %>%
  mutate(UpsReg = factor(UpsReg, levels = rev(.$UpsReg))) %>%
  pivot_longer(cols = !UpsReg&!type, names_to = c("metric", "comp"), names_sep = "__", values_to = "value") %>%
  pivot_wider(names_from = "metric", values_from = "value") %>%
  mutate(logp = -log10(pval_adj))
p.chems <- chems %>%
  ggplot(aes(x=comp, y=UpsReg)) + 
  geom_point(aes(size = logp, color = biascorr.zscore), shape=15) + 
  scale_size_continuous(limits = c(0.01,15), breaks = c(1, 5, 10)) +
  scale_color_gradientn(colours = c("blue","blue","black","yellow", "yellow"), na.value = "grey80",
                        values = scales::rescale(c(-3.5, -1, 0, 1, 3.5), from = c(-3.5, 3.5)), 
                        limits = c(-3.5, 3.5),
                        breaks = c(-1, 0, 1),
                        guide = guide_colorbar(title  = "Activation z-score", title.vjust = 0.9, frame.colour = "black", frame.linewidth = .2*2.8, ticks.linewidth = .2*2.8, ticks.colour = "black", draw.ulim = F, draw.llim = F, title.position = "top", direction = "vertical", barheight = unit(2, "char"), barwidth = unit(.5, "char"))) +
  theme_classic(base_size = 6, base_family = "Arial") +
  theme(plot.background = element_rect(fill="transparent", colour=NA),
        panel.background = element_rect(fill="transparent", colour=NA),
        #aspect.ratio = 1/1,
        plot.margin = margin(.5,.5,.5,.5, "line"),
        axis.title.y = element_blank(),  axis.title.x = element_blank(), 
        axis.text = element_text(size = 5, colour = "black"),
        strip.text = element_text(size = 5, color = "black"), strip.background = element_blank(),
        legend.position = "right", legend.justification = c("center", "center"), legend.background = element_rect(fill="transparent"), legend.box.background = element_rect(fill="transparent", color="transparent"),
        legend.text = element_text(size = 5), legend.margin = margin(0,0,0,0), legend.box.margin = margin(0,0,0,0), legend.spacing = unit(-.1, "char"), legend.title = element_text(size=5))


genes <- df.UpsReg %>%
  dplyr::filter(grepl("chemical|drug|microRNA", type)==F) %>%
  head(12) %>%
  mutate(UpsReg = factor(UpsReg, levels = rev(.$UpsReg))) %>%
  pivot_longer(cols = !UpsReg&!type, names_to = c("metric", "comp"), names_sep = "__", values_to = "value") %>%
  pivot_wider(names_from = "metric", values_from = "value") %>%
  mutate(logp = -log10(pval_adj))
p.genes <- genes %>%
  ggplot(aes(x=comp, y=UpsReg)) + 
  geom_point(aes(size = logp, color = biascorr.zscore), shape=15) + 
  scale_size_continuous(limits = c(0.01,15), breaks = c(1, 5, 10)) +
  scale_color_gradientn(colours = c("blue","blue","black","yellow", "yellow"), na.value = "grey80",
                        values = scales::rescale(c(-3.5, -1, 0, 1, 3.5), from = c(-3.5, 3.5)), 
                        limits = c(-3.5, 3.5),
                        breaks = c(-1, 0, 1),
                        guide = guide_colorbar(title  = "Activation z-score", title.vjust = 0.9, frame.colour = "black", frame.linewidth = .2*2.8, ticks.linewidth = .2*2.8, ticks.colour = "black", draw.ulim = F, draw.llim = F, title.position = "top", direction = "vertical", barheight = unit(2, "char"), barwidth = unit(.5, "char"))) +
  theme_classic(base_size = 6, base_family = "Arial") +
  theme(plot.background = element_rect(fill="transparent", colour=NA),
        panel.background = element_rect(fill="transparent", colour=NA),
        #aspect.ratio = 1/1,
        plot.margin = margin(.5,.5,.5,.5, "line"),
        axis.title.y = element_blank(),  axis.title.x = element_blank(), 
        axis.text = element_text(size = 5, colour = "black"),
        strip.text = element_text(size = 5, color = "black"), strip.background = element_blank(),
        legend.position = "right", legend.justification = c("center", "center"), legend.background = element_rect(fill="transparent"), legend.box.background = element_rect(fill="transparent", color="transparent"),
        legend.text = element_text(size = 5), legend.margin = margin(0,0,0,0), legend.box.margin = margin(0,0,0,0), legend.spacing = unit(-.1, "char"), legend.title = element_text(size=5))

plot_grid(p.chems, p.genes, align="hv", axis="tblr", nrow=2) %>%
  ggsave(filename = "/Users/arthur/Desktop/UpsReg.pdf", device = cairo_pdf, width = 15, height = 15, units = "cm")
