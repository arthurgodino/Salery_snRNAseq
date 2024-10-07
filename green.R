# Analysis of green cells only

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

integrated <- readRDS(file = "/Users/agodino/Desktop/Marine_integrated.rds")

############################################################################################################################################################
############################################################################################################################################################

# n cells per clusters per condition and seq depth
# Table of number of nuclei by cluster by conditions
n.cells <- table(integrated$cluster_sorted_pairings_cpp) %>% as.data.frame() %>% setNames(c("cluster_sorted_pairings_cpp", "n")) %>%
  separate(col = "cluster_sorted_pairings_cpp", into = c("cluster", "sorted", "pairings", "cpp"), sep = "_") %>%
  pivot_wider(id_cols = c("sorted", "pairings","cpp"), names_from = cluster, values_from = n) %>%
  replace(is.na(.), 0)
write.csv(n.cells, file = "/Users/arthur/Desktop/n_nuclei_per_cluster_and_condition.csv", row.names = F)

qc_to_plot <- FetchData(integrated, vars = c("comb.clusters", "sorted_pairings_cpp", "cpp", "nFeature_SCT", "nCount_SCT", "nFeature_RNA", "nCount_RNA"))
plist <- list()
for (n in names(qc_to_plot)[4:7]) {
plist[[n]] <- ggplot(qc_to_plot, aes(x = comb.clusters, y = !!sym(n), fill=sorted_pairings_cpp, color = cpp)) +
  geom_violin(position=position_dodge(width=.8), scale = "width", drop=F) +
  geom_jitter(position = position_jitterdodge(dodge.width = .8, jitter.width = .05, jitter.height = 0), size=.1, alpha=1, shape = 16) +
  scale_fill_manual(values = c("grey80","grey80","grey30","grey30","green2","green2","green4","green4")) +
  scale_color_manual(values = c("black", "dodgerblue")) +
  scale_y_continuous(limits=c(0,NA), expand = c(0,0)) +
  scale_x_discrete(limits=rev) +
  theme_classic(base_size = 6, base_family = "Helvetica")
}
ggsave(cowplot::plot_grid(plotlist = plist), file = "/Users/agodino/Desktop/vln.pdf", device = pdf, width = 18, height = 8, units = "in")
ggsave(cowplot::plot_grid(plotlist = plist), file = "/Users/agodino/Desktop/vln.png", device = png, width = 18, height = 8, units = "in")


UMIcounts <- full_join(mutate(as.data.frame(integrated$cluster_sorted_pairings_cpp), cell = rownames(as.data.frame(integrated$cluster_sorted_pairings_cpp))), 
                       mutate(as.data.frame(integrated$nCount_RNA), cell = rownames(as.data.frame(integrated$nCount_RNA)))) %>%
  setNames(c("cluster_sorted_pairings_cpp", "cell", "nCount_RNA")) %>%
  group_by(cluster_sorted_pairings_cpp) %>%
  summarise(median = median(nCount_RNA)) %>% #, mean = mean(nCount_RNA), mad = mad(nCount_RNA), sd = sd(nCount_RNA)
  separate(col = "cluster_sorted_pairings_cpp", into = c("cluster", "sorted", "pairings", "cpp"), sep = "_") %>%
  pivot_wider(id_cols = c("sorted", "pairings","cpp"), names_from = cluster, values_from = median)
write.csv(UMIcounts, file = "/Users/arthur/Desktop/median_nCount_RNA_per_cluster_and_condition.csv", row.names = F)
## -> problem with seq depth, need to use sct normalized assay for dea

## DEG by cluster between conditions - focus on HC only (activated vs non activated, 1P vs 5P)

# UMAP of green cells, HC+Test, 1Pv5P
Idents(integrated) <- "cpp"
integrated.HC <- subset(integrated, idents = "HC")
Idents(integrated.HC) <- "comb.clusters"
umap.HC.clust <- UMAPPlot(integrated.HC, label = T, pt.size = .1, label.size = 5/2.8, cols = setNames(scales::hue_pal()(6), levels(integrated$comb.clusters))) + 
  theme_classic(base_size = 5, base_family = "Arial") +
  scale_x_continuous(limits = c(-14,21), expand = c(0,0)) +
  scale_y_continuous(limits = c(-21,14), expand = c(0,0)) +
  theme(plot.background = element_rect(fill="transparent", colour=NA),
        panel.background = element_rect(fill="transparent", colour=NA),
        plot.margin = margin(.5,.5,.5,.5, "char"),
        aspect.ratio = 1, axis.text = element_text(color = "black"),
        legend.position = "none")
ggsave(umap.HC.clust, filename = "/Users/arthur/Desktop/clustering_HC_UMAP.pdf", device = cairo_pdf, width = 15, height = 15, units = "cm")

Idents(integrated.HC) <- "sorted_pairings"
umap.HC.green <- UMAPPlot(integrated.HC, cols = c("grey85", "grey75", "green2", "green4"), order = c("POS_5P", "POS_1P", "NEG_5P", "NEG_1P"), shuffle = F) +
  theme_classic(base_size = 5, base_family = "Arial") +
  scale_x_continuous(limits = c(-14,21), expand = c(0,0)) +
  scale_y_continuous(limits = c(-21,14), expand = c(0,0)) +
  theme(plot.background = element_rect(fill="transparent", colour=NA),
        panel.background = element_rect(fill="transparent", colour=NA),
        plot.margin = margin(.5,.5,.5,.5, "char"),
        aspect.ratio = 1, axis.text = element_text(color = "black"),
        legend.position = c(.9,.9))
ggsave(umap.HC.green, filename = "/Users/arthur/Desktop/green_HC_UMAP.pdf", device = cairo_pdf, width = 15, height = 15, units = "cm")

Idents(integrated) <- "sorted_pairings"
umap.all.green <- UMAPPlot(integrated, cols = c("grey85", "grey75", "green2", "green4"), order = c("POS_5P", "POS_1P", "NEG_5P", "NEG_1P"), shuffle = F) +
  theme_classic(base_size = 5, base_family = "Arial") +
  scale_x_continuous(limits = c(-14,21), expand = c(0,0)) +
  scale_y_continuous(limits = c(-21,14), expand = c(0,0)) +
  theme(plot.background = element_rect(fill="transparent", colour=NA),
        panel.background = element_rect(fill="transparent", colour=NA),
        plot.title = element_blank(),
        plot.margin = margin(.5,.5,.5,.5, "char"),
        aspect.ratio = 1, axis.text = element_text(color = "black"),
        legend.position = c(.9,.9))
ggsave(umap.all.green, filename = "/Users/arthur/Desktop/green_UMAP.pdf", device = cairo_pdf, width = 15, height = 15, units = "cm")


# Prop.graph all green
temp2 <- integrated@meta.data %>% group_by(sorted_pairings, comb.clusters) %>% summarize(n = n()) %>%
  pivot_wider(names_from = sorted_pairings, values_from = n) %>% replace(is.na(.), 0)
clu2 <- temp2$comb.clusters ; temp2$comb.clusters = NULL 
temp2 <- temp2 %>% as.data.frame() %>% apply(., 2, function(x){100*x/sum(x)}) %>% as.data.frame() %>% mutate(Cell_type = clu2)

library(ggradar)
radar.green <- temp2 %>%
  pivot_longer(cols=!Cell_type, names_to = "group", values_to = "Percentage") %>%
  #separate(Group, into = c("sorted", "pairings"), sep = "_") %>%
  mutate(Cell_type = factor(Cell_type, levels = c("D1-MSNs", "D2-MSNs", "Interneurons", "Astrocytes", "Oligodendrocytes"))) %>%
  pivot_wider(names_from = Cell_type, values_from = Percentage) %>%
  dplyr::select(c("group", rev(names(.)))) %>%
  ggradar(grid.min = 0, grid.mid = 35, grid.max = 70,
          base.size = 5, font.radar = "Arial",
          legend.text.size = 5,
          values.radar = c("0%", "35%", "70%"),
          group.colours = c("grey85","grey75","green2","green4"),
          group.line.width = 1, group.point.size = 2,
          background.circle.colour = "transparent",
          gridline.min.colour = "grey", gridline.mid.colour = "grey", gridline.max.colour = "grey")
ggsave(radar.green, filename = "/Users/arthur/Desktop/green_radar.pdf", device = cairo_pdf, width = 15, height = 15, units = "cm") 

# DEGs only using HC to avoid confounding effects of test on transcriptome (activity-dep + priming)

allDEGs_POSvNEG_5Pv1P_HC <- function(integrated, cluster = "D1-MSNs", lfcT = .25, minpct = .25, test = "LR", p.th = .05, outwd = "/Users/arthur/Desktop/") {
  print(paste0("Working on ", cluster))
  Idents(integrated) <- "cluster_sorted_pairings_cpp"
  comp <- list()
  comp[["POSvNEG_1P"]] <- FindMarkers(integrated, ident.1 = paste0(cluster, "_POS_1P_HC"), ident.2 = paste0(cluster, "_NEG_1P_HC"), assay = "SCT",  min.pct = minpct, logfc.threshold = lfcT, test.use = test, verbose = FALSE) %>% as.data.frame() %>% mutate(gene = rownames(.)) %>%
    mutate(p_val_adj = p.adjust(.$p_val, "fdr"))
  comp[["POSvNEG_5P"]] <- FindMarkers(integrated, ident.1 = paste0(cluster, "_POS_5P_HC"), ident.2 = paste0(cluster, "_NEG_5P_HC"), assay = "SCT",  min.pct = minpct, logfc.threshold = lfcT, test.use = test, verbose = FALSE) %>% as.data.frame() %>% mutate(gene = rownames(.)) %>%
    mutate(p_val_adj = p.adjust(.$p_val, "fdr"))
  comp[["NEG_5Pv1P"]] <- FindMarkers(integrated, ident.1 = paste0(cluster, "_NEG_5P_HC"), ident.2 = paste0(cluster, "_NEG_1P_HC"), assay = "SCT",  min.pct = minpct, logfc.threshold = lfcT, test.use = test, verbose = FALSE) %>% as.data.frame() %>% mutate(gene = rownames(.)) %>%
    mutate(p_val_adj = p.adjust(.$p_val, "fdr"))
  comp[["POS_5Pv1P"]] <- FindMarkers(integrated, ident.1 = paste0(cluster, "_POS_5P_HC"), ident.2 = paste0(cluster, "_POS_1P_HC"), assay = "SCT",  min.pct = minpct, logfc.threshold = lfcT, test.use = test, verbose = FALSE) %>% as.data.frame() %>% mutate(gene = rownames(.)) %>%
    mutate(p_val_adj = p.adjust(.$p_val, "fdr"))
  comp[["POS-NEG_5Pv1P"]] <- full_join(comp[["NEG_5Pv1P"]], comp[["POS_5Pv1P"]], by = "gene", suffix = c(".neg", "")) %>%
    dplyr::filter((p_val_adj.neg > p.th | is.na(p_val_adj.neg)==T) & p_val_adj < p.th) %>%
    dplyr::select(!matches("neg"))
  
  wb <- createWorkbook()                                                    
  sheetnames <- names(comp)
  sheets <- lapply(sheetnames, createSheet, wb = wb)
  void <- Map(addDataFrame, comp, sheets)
  saveWorkbook(wb, file = paste0(outwd, cluster, "_POSvNEG_5Pv1P_HConly_lfcT", lfcT, "_minpct", minpct, "_", test, ".xlsx"))
  
  d1 <- full_join(comp[["POSvNEG_1P"]], comp[["POSvNEG_5P"]], by = "gene", suffix = c(".posvneg.1P", ".posvneg.5P")) %>%
    full_join(setNames(comp[["NEG_5Pv1P"]], paste0(names(comp[["NEG_5Pv1P"]]), ".5Pv1Pneg")), by = join_by("gene"=="gene.5Pv1Pneg")) %>%
    full_join(setNames(comp[["POS_5Pv1P"]], paste0(names(comp[["POS_5Pv1P"]]), ".5Pv1Ppos")), by = join_by("gene"=="gene.5Pv1Ppos")) %>%
    full_join(setNames(comp[["POS-NEG_5Pv1P"]], paste0(names(comp[["POS-NEG_5Pv1P"]]), ".5Pv1Pposmneg")), by = join_by("gene"=="gene.5Pv1Pposmneg")) %>%
    dplyr::select(gene, avg_log2FC.posvneg.1P, p_val_adj.posvneg.1P, avg_log2FC.posvneg.5P, p_val_adj.posvneg.5P, avg_log2FC.5Pv1Pneg, p_val_adj.5Pv1Pneg, avg_log2FC.5Pv1Ppos, p_val_adj.5Pv1Ppos, avg_log2FC.5Pv1Pposmneg, p_val_adj.5Pv1Pposmneg) %>%
    setNames(c("gene", "lfc_posvneg.1P", "p_posvneg.1P", "lfc_posvneg.5P", "p_posvneg.5P", "lfc_neg.5Pv1P", "p_neg.5Pv1P", "lfc_pos.5Pv1P", "p_pos.5Pv1P","lfc_posmneg.5Pv1P", "p_posmneg.5Pv1P")) %>%
    mutate(cluster = cluster)
  
  DEGcount.posvneg.1P <- mutate(dplyr::filter(comp[["POSvNEG_1P"]], p_val_adj < p.th), updown = as.character(ifelse(avg_log2FC > 0, "up", "down"))) %>% group_by(updown) %>% summarise(count = n()) %>% mutate(cluster = cluster, pairwise = "POSvNEG_1P")
  DEGcount.posvneg.5P <- mutate(dplyr::filter(comp[["POSvNEG_5P"]], p_val_adj < p.th), updown = as.character(ifelse(avg_log2FC > 0, "up", "down"))) %>% group_by(updown) %>% summarise(count = n()) %>% mutate(cluster = cluster, pairwise = "POSvNEG_5P")
  DEGcount.neg.1Pv5P <- mutate(dplyr::filter(comp[["NEG_5Pv1P"]], p_val_adj < p.th), updown = as.character(ifelse(avg_log2FC > 0, "up", "down"))) %>% group_by(updown) %>% summarise(count = n()) %>% mutate(cluster = cluster, pairwise = "NEG_5Pv1P")
  DEGcount.pos.1Pv5P <- mutate(dplyr::filter(comp[["POS_5Pv1P"]], p_val_adj < p.th), updown = as.character(ifelse(avg_log2FC > 0, "up", "down"))) %>% group_by(updown) %>% summarise(count = n()) %>% mutate(cluster = cluster, pairwise = "POS_5Pv1P")
  DEGcount.posmneg.1Pv5P <- mutate(dplyr::filter(comp[["POS-NEG_5Pv1P"]], p_val_adj < p.th), updown = as.character(ifelse(avg_log2FC > 0, "up", "down"))) %>% group_by(updown) %>% summarise(count = n()) %>% mutate(cluster = cluster, pairwise = "POS-NEG_5Pv1P")
  DEGcounts <- full_join(DEGcount.posvneg.1P, DEGcount.posvneg.5P) %>% full_join(DEGcount.neg.1Pv5P) %>% full_join(DEGcount.pos.1Pv5P) %>% full_join(DEGcount.posmneg.1Pv5P)
  
  output <- list(exported = comp, DEGs = d1, counts = DEGcounts)
  return(output)
}

# For all clusters
comps.HC.LR <- list()
comps.HC.LR.nolfcT <- list()
for (i in c("D1-MSNs", "D2-MSNs")) {
  comps.HC.LR[[i]] <- allDEGs_POSvNEG_5Pv1P_HC(integrated, cluster = i, lfcT = log2(1.15), minpct = .3, test = "LR", p.th = .1)
  comps.HC.LR.nolfcT[[i]] <- allDEGs_POSvNEG_5Pv1P_HC(integrated, cluster = i, lfcT = -Inf, minpct = 0.01, test = "LR", p.th = .1)
}

## to use, HC only comparisons with LR testing
final <- comps.HC.LR 
counts <- lapply(final, pluck, "counts") %>% data.table::rbindlist(idcol = "Cluster") %>%
  mutate(updown = factor(updown, levels = c("up", "down")), 
         pairwise = factor(pairwise, levels = c("POSvNEG_1P", "POSvNEG_5P", "NEG_5Pv1P", "POS_5Pv1P", "POS-NEG_5Pv1P")),
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

plot_gene_HC <- function(integrated, gene = "Kcnd2", cluster = "D1-MSNs") {
  Idents(integrated) <- "cluster_cpp"
  d1cells <- subset(integrated, idents = paste0(cluster, "_HC"))
  p <- VlnPlot(d1cells, features = c(gene), split.by = "cluster_sorted_pairings_cpp", assay = "SCT", slot = "data", pt.size = 0) +
    scale_fill_manual(values = c("grey80","grey30","green2","green4")) +
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

d1.plots = list()
for (g in dplyr::filter(final$`D1-MSNs`$DEGs, p_posmneg.5Pv1P <.1)$gene) {
  d1.plots[[g]] <- plot_gene_HC(integrated, gene = g, cluster = "D1-MSNs")
}

plot_grid(plotlist = d1.plots, ncol = 8, nrow = 10, align = "hv", axis = "trbl") %>%
  ggsave(filename = "/Users/arthur/Desktop/D1_DEGs.pdf", device = cairo_pdf, height = 37, width = 21.6, units = "cm")

# final figure
plot_grid(umap.all.green, NULL, plot_grid(p.counts, NULL, ncol = 1, align = "hv", axis = "tblr"), plot_grid(d1.plots[["Reln"]], d1.plots[["Kcnip4"]], d1.plots[["Prkca"]],d1.plots[["Pde7b"]], nrow=2, align = "hv", axis = "tblr"), nrow = 2, align = "hv", axis = "trbl") %>%
  ggsave(filename = "/Users/arthur/Desktop/green.pdf", device = cairo_pdf, width = 21.6, height = 27.9, units = "cm") 

# Heatmaps
plot.heatmap <- function(df, ttl) {
  df %>%
    dplyr::select(gene, lfc_posvneg.1P, lfc_posvneg.5P, lfc_pos.5Pv1P) %>%
    pivot_longer(cols = !gene, names_to = "comp", values_to = "log2FC") %>%
    mutate(comp = factor(case_when(comp == "lfc_posvneg.1P" ~ "1P - POSvNEG", comp == "lfc_posvneg.5P" ~ "5P - POSvNEG", comp == "lfc_pos.5Pv1P" ~ "POS - 5Pv1P"))) %>%
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

hm.d1.1P <- data.frame(gene = dplyr::filter(comps.HC.LR$`D1-MSNs`$DEGs,p_posvneg.1P < .1)$gene) %>%
  left_join(comps.HC.LR.nolfcT$`D1-MSNs`$DEGs) %>%
  mutate(gene = factor(gene, levels = arrange(., lfc_posvneg.1P)$gene)) %>%
  plot.heatmap(ttl = "D1 - POSvNEG_1P") # + theme(axis.text.x = element_text(size = 5, colour = "black", angle = 90, hjust=0.95,vjust=0.2))
hm.d1.5P <- data.frame(gene = dplyr::filter(comps.HC.LR$`D1-MSNs`$DEGs,p_posvneg.5P < .1)$gene) %>%
  left_join(comps.HC.LR.nolfcT$`D1-MSNs`$DEGs) %>%
  mutate(gene = factor(gene, levels = arrange(., lfc_posvneg.5P)$gene)) %>%
  plot.heatmap(ttl = "D1 - POSvNEG_5P") # + theme(axis.text.x = element_text(size = 5, colour = "black", angle = 90, hjust=0.95,vjust=0.2))
hm.d1.POS <- data.frame(gene = dplyr::filter(comps.HC.LR$`D1-MSNs`$DEGs,p_posmneg.5Pv1P < .1)$gene) %>%
  left_join(comps.HC.LR.nolfcT$`D1-MSNs`$DEGs) %>%
  mutate(gene = factor(gene, levels = arrange(., lfc_pos.5Pv1P)$gene)) %>%
  plot.heatmap(ttl = "D1 - POSmNEG_5Pv1P")# + theme(axis.text.x = element_text(size = 5, colour = "black", angle = 90, hjust=0.95,vjust=0.2))


hm.D2.1P <- data.frame(gene = dplyr::filter(comps.HC.LR$`D2-MSNs`$DEGs,p_posvneg.1P < .1)$gene) %>%
  left_join(comps.HC.LR.nolfcT$`D2-MSNs`$DEGs) %>%
  mutate(gene = factor(gene, levels = arrange(., lfc_posvneg.1P)$gene)) %>%
  plot.heatmap(ttl = "D2 - POSvNEG_1P") # + theme(axis.text.x = element_text(size = 5, colour = "black", angle = 90, hjust=0.95,vjust=0.2))
hm.D2.5P <- data.frame(gene = dplyr::filter(comps.HC.LR$`D2-MSNs`$DEGs,p_posvneg.5P < .1)$gene) %>%
  left_join(comps.HC.LR.nolfcT$`D2-MSNs`$DEGs) %>%
  mutate(gene = factor(gene, levels = arrange(., lfc_posvneg.5P)$gene)) %>%
  plot.heatmap(ttl = "D2 - POSvNEG_5P") # + theme(axis.text.x = element_text(size = 5, colour = "black", angle = 90, hjust=0.95,vjust=0.2))
hm.D2.POS <- data.frame(gene = dplyr::filter(comps.HC.LR$`D2-MSNs`$DEGs,p_posmneg.5Pv1P < .1)$gene) %>%
  left_join(comps.HC.LR.nolfcT$`D2-MSNs`$DEGs) %>%
  mutate(gene = factor(gene, levels = arrange(., lfc_pos.5Pv1P)$gene)) %>%
  plot.heatmap(ttl = "D2 - POSmNEG_5Pv1P")# + theme(axis.text.x = element_text(size = 5, colour = "black", angle = 90, hjust=0.95,vjust=0.2))



all.heatmaps <- list(hm.d1.1P, hm.d1.5P, hm.d1.POS, hm.D2.1P, hm.D2.5P, hm.D2.POS) ; rm(hm.d1.1P, hm.d1.5P, hm.d1.POS, hm.D2.1P, hm.D2.5P, hm.D2.POS)
maxgenes <- lapply(all.heatmaps, function(x) {nrow(x$data)/3}) %>% unlist() %>% max()

addextra <- function(t, maxadd) {
  toadd <- maxadd - (nrow(t$data)/3)
  if (toadd > 0) {
    tobind <- data.frame(gene = sort(rep(paste0("extra", seq(1,toadd)), 3)), comp = rep(c("1P - POSvNEG", "5P - POSvNEG", "POS - 5Pv1P"), toadd), log2FC = rep(NA, toadd*3))
    t$data <- t$data %>% rbind(tobind)
  }
  return(t)
}

all.heatmaps <- lapply(all.heatmaps, addextra, maxadd = maxgenes)

plot_grid(plotlist = all.heatmaps, ncol = 2) %>%
  ggsave(filename = "/Users/arthur/Desktop/green_heatmaps.pdf", device = cairo_pdf, width = 21.6, height = 27.9, units = "cm") 



# D1 vs D2 comparisons
#total genes = 19798

D1.1Psig.5Psig <- comps.HC.LR$`D1-MSNs`$DEGs %>%
                  dplyr::filter(p_posvneg.1P < .1 & p_posvneg.5P < .1 & sign(lfc_posvneg.1P)==sign(lfc_posvneg.5P)) %>%
                  dplyr::select(1:5)
write.csv(D1.1Psig.5Psig, file ="/Users/arthur/Desktop/Overlap_D1_1Pn5P.csv", row.names = F)

D2.1Psig.5Psig <- comps.HC.LR$`D2-MSNs`$DEGs %>%
  dplyr::filter(p_posvneg.1P < .1 & p_posvneg.5P < .1 & sign(lfc_posvneg.1P)==sign(lfc_posvneg.5P)) %>%
  dplyr::select(1:5)
write.csv(D2.1Psig.5Psig, file ="/Users/arthur/Desktop/Overlap_D2_1Pn5P.csv", row.names = F)
  
fish.1P <- fisher.test(matrix(c(36, 93, 158, 19798-36-93-158), nrow=2), alternative = "greater")
fish.5P <- fisher.test(matrix(c(113, 136, 270, 19798-113-136-270), nrow=2), alternative = "greater")
fish.POS <- fisher.test(matrix(c(1, 86, 18, 19798-1-86-18), nrow=2), alternative = "greater")

fish.d1 <- fisher.test(matrix(c(53, 76, 196, 19798-56-73-193), nrow=2), alternative = "greater")
fish.d2 <- fisher.test(matrix(c(49, 145, 334, 19798-49-145-334), nrow=2), alternative = "greater")

# Dot plots example genes
Idents(integrated.HC) <- "comb.clusters"
tempD1 <- subset(integrated.HC, idents = "D1-MSNs")

ex <- DotPlot(tempD1, features = rev(c("Gria3", "Grik5", "Kcnip4", "Reln", "Rap1gap", "Pde7b", "Rgs20", "Prkca")), 
                      #cluster.idents = F,
                      assay = "SCT", dot.scale = 6, group.by = "cluster_sorted_pairings_cpp",
                      cols = c("lightgrey", "#F564E3"),
                      dot.min = .01, col.min = 0, col.max = 1, scale.min = .1) + 
  coord_flip() +
  #scale_y_discrete(labels = levels(integrated$comb.clusters)) +
  theme_classic(base_family = "Arial", base_size = 5) +
  theme(plot.background = element_rect(fill="transparent", colour=NA),
        panel.background = element_rect(fill="transparent", colour=NA),
        aspect.ratio = 8/4, axis.title = element_blank(), axis.text = element_text(color = "black"), axis.ticks = element_line(color="black"),
        axis.text.x = element_text(face = "italic", angle = 45, vjust = 0.8, hjust = 1),
        legend.key.size = unit(.5, "line"), legend.position = "bottom", legend.box = "vertical")
ggsave(ex, filename = "/Users/arthur/Desktop/exD1_DotPlot.pdf", device = cairo_pdf, width = 9, height = 9, units = "cm") 

# IPA upstream
setwd("/Users/arthur/Desktop/IPA_snRNAseq/Genes_n_Chem/")
allUpsReg <- lapply(list.files(pattern=".csv"), read_csv)
names(allUpsReg) <- list.files(pattern=".csv")
df.UpsReg <- rbindlist(allUpsReg, idcol = "File", fill=T) %>%
          mutate(comp = gsub("_UpsReg.csv", "", File), File=NULL) %>%
          dplyr::select(comp, 1, 3,5, 6, 8,9, 10) %>%
          setNames(c("comp", "UpsReg", "type", "biascorr.zscore", "zscore", "pval", "pval_adj", "targets")) %>%
          mutate(ngenes = 1+str_count(targets, ",")) %>%
          pivot_wider(id_cols = c("UpsReg", "type"), names_sep = "__", names_from = "comp", values_from = c("biascorr.zscore", "zscore", "pval", "pval_adj", "targets", "ngenes")) %>%
  mutate(n.NA = rowSums(is.na(.[,15:20]))) %>%
  mutate(ngenes_D1_1P_POSvNEG = ifelse(is.na(ngenes__D1_1P_POSvNEG)==T, 0, ngenes__D1_1P_POSvNEG), 
         ngenes_D1_5P_POSvNEG = ifelse(is.na(ngenes__D1_5P_POSvNEG)==T, 0, ngenes__D1_5P_POSvNEG),
         ngenes_D1_5Pv1P_POS = ifelse(is.na(ngenes__D1_5Pv1P_POS)==T, 0, ngenes__D1_5Pv1P_POS),
         ngenes_D2_1P_POSvNEG = ifelse(is.na(ngenes__D2_1P_POSvNEG)==T, 0, ngenes__D2_1P_POSvNEG), 
         ngenes_D2_5P_POSvNEG = ifelse(is.na(ngenes__D2_5P_POSvNEG)==T, 0, ngenes__D2_5P_POSvNEG),
         ngenes_D2_5Pv1P_POS = ifelse(is.na(ngenes__D2_5Pv1P_POS)==T, 0, ngenes__D2_5Pv1P_POS)) %>%
  rowwise() %>%
          mutate(ordersig = sum(pval_adj__D1_1P_POSvNEG, pval_adj__D1_5P_POSvNEG, pval_adj__D1_5Pv1P_POS, pval_adj__D2_1P_POSvNEG, pval_adj__D2_5P_POSvNEG, pval_adj__D2_5Pv1P_POS, na.rm =T),
                 avg.ngenes = sum(ngenes__D1_1P_POSvNEG, ngenes__D1_5P_POSvNEG, ngenes__D1_5Pv1P_POS, ngenes__D2_1P_POSvNEG, ngenes__D2_5P_POSvNEG, ngenes__D2_5Pv1P_POS, na.rm =T)/6) %>%
          arrange(n.NA, ordersig)%>%
  dplyr::filter(avg.ngenes >= 3) %>%
          dplyr::select(1:8, 21:26) %>%
          ungroup()

chems <- df.UpsReg %>%
        dplyr::filter(grepl("chemical|drug", type)==T) %>%
        head(12) %>%
        mutate(UpsReg = factor(UpsReg, levels = rev(.$UpsReg))) %>%
        pivot_longer(cols = !UpsReg&!type, names_to = c("metric", "comp"), names_sep = "__", values_to = "value") %>%
        pivot_wider(names_from = "metric", values_from = "value") %>%
        mutate(logp = -log10(pval_adj))
p.chems <- chems %>%
           ggplot(aes(x=comp, y=UpsReg)) + 
           geom_point(aes(size = logp, color = biascorr.zscore), shape=15) + 
  scale_size_continuous(limits = c(0.01,10), breaks = c(1, 5, 10)) +
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
  scale_size_continuous(limits = c(0.01,10), breaks = c(1, 5, 10)) +
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
