# Analysis of yellow cells only

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
out.wd = "/Users/arthur/Desktop/"
############################################################################################################################################################
############################################################################################################################################################

### Analysis of reactivation

# UMAP of yellow cells

Idents(integrated) <- "cpp"
integrated.HC <- subset(integrated, idents = "HC")
integrated.Test <- subset(integrated, idents = "Test")

Idents(integrated.HC) <- "pairings_sorted_Arc"
umap.HC.yellow <- UMAPPlot(integrated.HC, cols = rev(c("yellow3","yellow1","red4","red2","green4","green2","grey75","grey85")),
                        order = c("5P_POS_YES", "1P_POS_YES", "5P_NEG_YES", "1P_NEG_YES", "5P_POS_NO", "1P_POS_NO", "5P_NEG_NO", "1P_NEG_NO"), shuffle = F) +
  theme_classic(base_size = 5, base_family = "Arial") +
  scale_x_continuous(limits = c(-14,21), expand = c(0,0)) +
  scale_y_continuous(limits = c(-21,14), expand = c(0,0)) +
  theme(plot.background = element_rect(fill="transparent", colour=NA),
        panel.background = element_rect(fill="transparent", colour=NA),
        plot.title = element_blank(),
        plot.margin = margin(.5,.5,.5,.5, "char"),
        aspect.ratio = 1, axis.text = element_text(color = "black"),
        legend.position = c(.9,.9))
ggsave(umap.HC.yellow, filename = "/Users/arthur/Desktop/yellow_HC_UMAP.pdf", device = cairo_pdf, width = 15, height = 15, units = "cm")

Idents(integrated.Test) <- "pairings_sorted_Arc"
umap.Test.yellow <- UMAPPlot(integrated.Test, cols = rev(c("yellow3","yellow1","red4","red2","green4","green2","grey75","grey85")),
                           order = c("5P_POS_YES", "1P_POS_YES", "5P_NEG_YES", "1P_NEG_YES", "5P_POS_NO", "1P_POS_NO", "5P_NEG_NO", "1P_NEG_NO"), shuffle = F) +
  theme_classic(base_size = 5, base_family = "Arial") +
  scale_x_continuous(limits = c(-14,21), expand = c(0,0)) +
  scale_y_continuous(limits = c(-21,14), expand = c(0,0)) +
  theme(plot.background = element_rect(fill="transparent", colour=NA),
        panel.background = element_rect(fill="transparent", colour=NA),
        plot.title = element_blank(),
        plot.margin = margin(.5,.5,.5,.5, "char"),
        aspect.ratio = 1, axis.text = element_text(color = "black"),
        legend.position = c(.9,.9))
ggsave(umap.Test.yellow, filename = "/Users/arthur/Desktop/yellow_Test_UMAP.pdf", device = cairo_pdf, width = 15, height = 15, units = "cm")


# Prop.graph yellow by HC/Test
library(ggradar)

temp2 <- integrated.HC@meta.data %>% group_by(pairings_sorted_Arc, comb.clusters) %>% summarize(n = n()) %>%
  pivot_wider(names_from = pairings_sorted_Arc, values_from = n) %>% replace(is.na(.), 0)
clu2 <- temp2$comb.clusters ; temp2$comb.clusters = NULL 
temp2 <- temp2 %>% as.data.frame() %>% apply(., 2, function(x){100*x/sum(x)}) %>% as.data.frame() %>% mutate(Cell_type = clu2)
radar.yellow.HC <- temp2 %>%
  pivot_longer(cols=!Cell_type, names_to = "group", values_to = "Percentage") %>%
  #separate(Group, into = c("sorted", "pairings"), sep = "_") %>%
  mutate(Cell_type = factor(Cell_type, levels = c("D1-MSNs", "D2-MSNs", "Interneurons", "Astrocytes", "Oligodendrocytes"))) %>%
  pivot_wider(names_from = Cell_type, values_from = Percentage) %>%
  dplyr::select(c("group", rev(names(.)))) %>%
  ggradar(grid.min = 0, grid.mid = 40, grid.max = 80,
          base.size = 5, font.radar = "Arial",
          legend.text.size = 5,
          values.radar = c("0%", "45%", "80%"),
          group.colours = c("grey85", "red2", "green2","yellow1","grey75","red4","green4","yellow3"),
          group.line.width = 1, group.point.size = 2,
          background.circle.colour = "transparent",
          gridline.min.colour = "grey", gridline.mid.colour = "grey", gridline.max.colour = "grey")
ggsave(radar.yellow.HC, filename = "/Users/arthur/Desktop/yellow_HC_radar.pdf", device = cairo_pdf, width = 15, height = 15, units = "cm") 

temp2 <- integrated.Test@meta.data %>% group_by(pairings_sorted_Arc, comb.clusters) %>% summarize(n = n()) %>%
  pivot_wider(names_from = pairings_sorted_Arc, values_from = n) %>% replace(is.na(.), 0)
clu2 <- temp2$comb.clusters ; temp2$comb.clusters = NULL 
temp2 <- temp2 %>% as.data.frame() %>% apply(., 2, function(x){100*x/sum(x)}) %>% as.data.frame() %>% mutate(Cell_type = clu2)
radar.yellow.Test <- temp2 %>%
  pivot_longer(cols=!Cell_type, names_to = "group", values_to = "Percentage") %>%
  #separate(Group, into = c("sorted", "pairings"), sep = "_") %>%
  mutate(Cell_type = factor(Cell_type, levels = c("D1-MSNs", "D2-MSNs", "Interneurons", "Astrocytes", "Oligodendrocytes"))) %>%
  pivot_wider(names_from = Cell_type, values_from = Percentage) %>%
  dplyr::select(c("group", rev(names(.)))) %>%
  ggradar(grid.min = 0, grid.mid = 40, grid.max = 80,
          base.size = 5, font.radar = "Arial",
          legend.text.size = 5,
          values.radar = c("0%", "40%", "80%"),
          group.colours = c("grey85", "red2", "green2","yellow1","grey75","red4","green4","yellow3"),
          group.line.width = 1, group.point.size = 2,
          background.circle.colour = "transparent",
          gridline.min.colour = "grey", gridline.mid.colour = "grey", gridline.max.colour = "grey")
ggsave(radar.yellow.Test, filename = "/Users/arthur/Desktop/red_yellow_radar.pdf", device = cairo_pdf, width = 15, height = 15, units = "cm") 

# Preferential reactivation

temp <- integrated@meta.data %>% 
  group_by(Group, comb.clusters) %>% summarize(n = n()) %>%
  pivot_wider(names_from = Group, values_from = n) %>% replace(is.na(.), 0) #%>%
#mutate(POS_1P = POS_1P_HC + POS_1P_Test, POS_5P = POS_5P_HC + POS_5P_Test, NEG_1P = NEG_1P_HC + NEG_1P_Test, NEG_5P = NEG_5P_HC + NEG_5P_Test) %>%
#select(comb.clusters, NEG_1P, NEG_5P, POS_1P, POS_5P)
clu <- temp$comb.clusters ; temp$comb.clusters = NULL 
PROP.POS <- temp %>% as.data.frame() %>% apply(., 2, function(x){100*x/sum(x)}) %>% as.data.frame() %>% mutate(Cell_type = clu) %>%
  pivot_longer(cols = !Cell_type, names_to = "group", values_to = "perc.total") %>%
  separate(group, into = c("sorted", "pairings", "cpp"))

temp <- integrated@meta.data %>% 
  group_by(pairings_cpp_sorted_Arc, comb.clusters) %>% summarize(n = n()) %>%
  pivot_wider(names_from = pairings_cpp_sorted_Arc, values_from = n) %>% replace(is.na(.), 0) %>%
  pivot_longer(!comb.clusters, names_to = "pairings_cpp_sorted_Arc", values_to = "n") %>%
  separate(pairings_cpp_sorted_Arc, into = c("pairings", "cpp", "sorted", "Arc"), sep = "_") %>%
  mutate(sortedArc = factor(paste0(sorted, Arc), levels = c("POSYES", "POSNO", "NEGYES", "NEGNO")),
         cluster_pairings_cpp = paste(comb.clusters, pairings, cpp, sep = "_")) %>%
  select(cluster_pairings_cpp, sortedArc, n) %>%
  split(.$cluster_pairings_cpp) %>%
  lapply(pivot_wider, names_from = "sortedArc", values_from = "n") %>%
  lapply(mutate, actrate.POS = 100*POSYES/(POSYES+POSNO), actrate.NEG = 100*NEGYES/(NEGYES+NEGNO)) %>%
  rbindlist() %>%
  mutate(n.POS = POSYES+POSNO, n.NEG = NEGYES+NEGNO) %>%
  select(cluster_pairings_cpp, actrate.POS, actrate.NEG, n.POS, n.NEG) %>%
  separate(cluster_pairings_cpp, into = c("Cell_type", "pairings", "cpp"), sep = "_") %>%
  left_join(PROP.POS)

radar.react.Test <- temp %>% filter(cpp == "Test") %>%
  select(-n.POS, -n.NEG) %>%
  mutate(act.rate = case_when(sorted == "NEG" ~ actrate.NEG, sorted == "POS" ~ actrate.POS)) %>%
  select(-matches("POS|NEG")) %>% replace(is.na(.), 0) %>%
  mutate(act.rate.bin = ifelse(sorted == "NEG", -act.rate, act.rate)) %>%
  mutate(Cell_type = factor(Cell_type, levels = c("D1-MSNs", "D2-MSNs", "Interneurons", "Astrocytes", "Oligodendrocytes"))) %>%
  mutate(sorted_pairings = paste(sorted, pairings, sep ="_")) %>%
  ggplot(aes(y = act.rate.bin, x = perc.total)) +
  facet_wrap(~Cell_type, nrow = 2) +
  geom_rect(xmin = 0, xmax = 80, ymin = -50, ymax = 0, fill = "grey85", alpha = .2) +
  geom_rect(xmin = 0, xmax = 80, ymin = 0, ymax = 50, fill = "green2", alpha = .2) +
  geom_hline(yintercept = c(-25, 0, 25, 50), linetype = "solid") +
  geom_vline(xintercept = c(0, 35, 70), linetype = rep(c("solid", "dotted", "dotted"),5)) +
  geom_bar(aes(fill = sorted_pairings), stat = "identity", width = 10) +
  coord_polar(theta = "y", start = pi) +
  scale_y_continuous(limits = c(-50, 50), breaks = c(-25, 0,25, 50), labels = c(25, 0, 25, 50)) +
  scale_x_continuous(limits = c(-10,80)) +
  scale_fill_manual(values = c("red2", "red4", "yellow1", "yellow3")) +
  theme(panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent"), strip.background = element_rect(fill = "transparent"),
        axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(color = "black"),
        panel.grid = element_blank(),
        axis.title = element_blank())
radar.react.HC <- temp %>% filter(cpp == "HC") %>%
  select(-n.POS, -n.NEG) %>%
  mutate(act.rate = case_when(sorted == "NEG" ~ actrate.NEG, sorted == "POS" ~ actrate.POS)) %>%
  select(-matches("POS|NEG")) %>% replace(is.na(.), 0) %>%
  mutate(act.rate.bin = ifelse(sorted == "NEG", -act.rate, act.rate)) %>%
  mutate(Cell_type = factor(Cell_type, levels = c("D1-MSNs", "D2-MSNs", "Interneurons", "Astrocytes", "Oligodendrocytes"))) %>%
  mutate(sorted_pairings = paste(sorted, pairings, sep ="_")) %>%
  ggplot(aes(y = act.rate.bin, x = perc.total)) +
  facet_wrap(~Cell_type, nrow = 2) +
  geom_rect(xmin = 0, xmax = 80, ymin = -50, ymax = 0, fill = "grey85", alpha = .2) +
  geom_rect(xmin = 0, xmax = 80, ymin = 0, ymax = 50, fill = "green2", alpha = .2) +
  geom_hline(yintercept = c(-25, 0, 25, 50), linetype = "solid") +
  geom_vline(xintercept = c(0, 35, 70), linetype = rep(c("solid", "dotted", "dotted"),5)) +
  geom_bar(aes(fill = sorted_pairings), stat = "identity", width = 10) +
  coord_polar(theta = "y", start = pi) +
  scale_y_continuous(limits = c(-50, 50), breaks = c(-25, 0,25, 50), labels = c(25, 0, 25, 50)) +
  scale_x_continuous(limits = c(-10,80)) +
  scale_fill_manual(values = c("red2", "red4", "yellow2", "yellow4")) +
  theme(panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent"), strip.background = element_rect(fill = "transparent"),
        axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(color = "black"),
        panel.grid = element_blank(),
        axis.title = element_blank())

plot_grid(radar.react.Test, radar.react.HC, nrow = 2) %>%
  ggsave(filename = "/Users/arthur/Desktop/radarreact.pdf", device = cairo_pdf, width = 20, height = 20, units = "cm")


# Stats chi square on counts
temp <- integrated@meta.data %>% 
  group_by(pairings_cpp_sorted_Arc, comb.clusters) %>% summarize(n = n()) %>%
  pivot_wider(names_from = pairings_cpp_sorted_Arc, values_from = n) %>% replace(is.na(.), 0) %>%
  pivot_longer(!comb.clusters, names_to = "pairings_cpp_sorted_Arc", values_to = "n") %>%
  separate(pairings_cpp_sorted_Arc, into = c("pairings", "cpp", "sorted", "Arc"), sep = "_") %>%
  mutate(sortedArc = factor(paste0(sorted, Arc), levels = c("POSYES", "POSNO", "NEGYES", "NEGNO")),
         cluster_pairings_cpp = paste(comb.clusters, pairings, cpp, sep = "_")) %>%
  select(cluster_pairings_cpp, sortedArc, n) %>%
  split(.$cluster_pairings_cpp) %>%
  lapply(pivot_wider, names_from = "sortedArc", values_from = "n") %>%
  lapply(mutate, actrate.POS = 100*POSYES/(POSYES+POSNO), actrate.NEG = 100*NEGYES/(NEGYES+NEGNO)) %>%
  rbindlist() %>%
  mutate(n.POS = POSYES+POSNO, n.NEG = NEGYES+NEGNO)

a <- temp %>% filter(cluster_pairings_cpp == "D1-MSNs_1P_Test") %>%
              select(2:5)
a <- rbind(c(a$NEGNO[1], a$POSNO[1]), c(a$NEGYES[1], a$POSYES[1])) %>% as.table()
capture.output(chisq.test(a), file = "/Users/arthur/Desktop/X2_D1_1P_Test.txt")
a <- temp %>% filter(cluster_pairings_cpp == "D1-MSNs_5P_Test") %>%
  select(2:5)
a <- rbind(c(a$NEGNO[1], a$POSNO[1]), c(a$NEGYES[1], a$POSYES[1])) %>% as.table()
capture.output(chisq.test(a), file = "/Users/arthur/Desktop/X2_D1_5P_Test.txt")
a <- temp %>% filter(cluster_pairings_cpp == "D2-MSNs_1P_Test") %>%
  select(2:5)
a <- rbind(c(a$NEGNO[1], a$POSNO[1]), c(a$NEGYES[1], a$POSYES[1])) %>% as.table()
capture.output(chisq.test(a), file = "/Users/arthur/Desktop/X2_D2_1P_Test.txt")
a <- temp %>% filter(cluster_pairings_cpp == "D2-MSNs_5P_Test") %>%
  select(2:5)
a <- rbind(c(a$NEGNO[1], a$POSNO[1]), c(a$NEGYES[1], a$POSYES[1])) %>% as.table()
capture.output(chisq.test(a), file = "/Users/arthur/Desktop/X2_D2_5P_Test.txt")



## DEG reactivation ######

# DEGs
allDEGs_POSYESvPOSNO_5Pv1P_Test <- function(integrated, cluster = "D1-MSNs", lfcT = .25, minpct = .25, test = "LR", p.th = .05, outwd = "/Users/arthur/Desktop/") {
  print(paste0("Working on ", cluster))
  Idents(integrated) <- "cluster_pairings_cpp_sorted_Arc"
  comp <- list()
  comp[["POSYESvPOSNO_1P"]] <- FindMarkers(integrated, ident.1 = paste0(cluster, "_1P_Test_POS_YES"), ident.2 = paste0(cluster, "_1P_Test_POS_NO"), assay = "SCT",  min.pct = minpct, logfc.threshold = lfcT, test.use = test, verbose = FALSE) %>% as.data.frame() %>% mutate(gene = rownames(.)) %>%
    mutate(p_val_adj = p.adjust(.$p_val, "fdr"))
  comp[["POSYESvPOSNO_5P"]] <- FindMarkers(integrated, ident.1 = paste0(cluster, "_5P_Test_POS_YES"), ident.2 = paste0(cluster, "_5P_Test_POS_NO"), assay = "SCT",  min.pct = minpct, logfc.threshold = lfcT, test.use = test, verbose = FALSE) %>% as.data.frame() %>% mutate(gene = rownames(.)) %>%
    mutate(p_val_adj = p.adjust(.$p_val, "fdr"))
  comp[["POSNO_5Pv1P"]] <- FindMarkers(integrated, ident.1 = paste0(cluster, "_5P_Test_POS_NO"), ident.2 = paste0(cluster, "_1P_Test_POS_NO"), assay = "SCT",  min.pct = minpct, logfc.threshold = lfcT, test.use = test, verbose = FALSE) %>% as.data.frame() %>% mutate(gene = rownames(.)) %>%
    mutate(p_val_adj = p.adjust(.$p_val, "fdr"))
  comp[["POSYES_5Pv1P"]] <- FindMarkers(integrated, ident.1 = paste0(cluster, "_5P_Test_POS_YES"), ident.2 = paste0(cluster, "_1P_Test_POS_YES"), assay = "SCT",  min.pct = minpct, logfc.threshold = lfcT, test.use = test, verbose = FALSE) %>% as.data.frame() %>% mutate(gene = rownames(.)) %>%
    mutate(p_val_adj = p.adjust(.$p_val, "fdr"))
  comp[["5Pv1P_Test_POSYESmPOSNO"]] <- full_join(comp[["POSNO_5Pv1P"]], comp[["POSYES_5Pv1P"]], by = "gene", suffix = c(".no", "")) %>%
    dplyr::filter((p_val_adj.no > p.th | is.na(p_val_adj.no)==T) & p_val_adj < p.th) %>%
    dplyr::select(!matches("no"))
  comp[["NEGYESvNEGNO_1P"]] <- FindMarkers(integrated, ident.1 = paste0(cluster, "_1P_Test_NEG_YES"), ident.2 = paste0(cluster, "_1P_Test_NEG_NO"), assay = "SCT",  min.pct = minpct, logfc.threshold = lfcT, test.use = test, verbose = FALSE) %>% as.data.frame() %>% mutate(gene = rownames(.)) %>%
    mutate(p_val_adj = p.adjust(.$p_val, "fdr"))
  comp[["NEGYESvNEGNO_5P"]] <- FindMarkers(integrated, ident.1 = paste0(cluster, "_5P_Test_NEG_YES"), ident.2 = paste0(cluster, "_5P_Test_NEG_NO"), assay = "SCT",  min.pct = minpct, logfc.threshold = lfcT, test.use = test, verbose = FALSE) %>% as.data.frame() %>% mutate(gene = rownames(.)) %>%
    mutate(p_val_adj = p.adjust(.$p_val, "fdr"))
  comp[["POSYESvNEGYES_1P"]] <- FindMarkers(integrated, ident.1 = paste0(cluster, "_1P_Test_POS_YES"), ident.2 = paste0(cluster, "_1P_Test_NEG_YES"), assay = "SCT",  min.pct = minpct, logfc.threshold = lfcT, test.use = test, verbose = FALSE) %>% as.data.frame() %>% mutate(gene = rownames(.)) %>%
    mutate(p_val_adj = p.adjust(.$p_val, "fdr"))
  comp[["POSYESvNEGYES_5P"]] <- FindMarkers(integrated, ident.1 = paste0(cluster, "_5P_Test_POS_YES"), ident.2 = paste0(cluster, "_5P_Test_NEG_YES"), assay = "SCT",  min.pct = minpct, logfc.threshold = lfcT, test.use = test, verbose = FALSE) %>% as.data.frame() %>% mutate(gene = rownames(.)) %>%
    mutate(p_val_adj = p.adjust(.$p_val, "fdr"))
  
  wb <- createWorkbook()                                                    
  sheetnames <- names(comp)
  sheets <- lapply(sheetnames, createSheet, wb = wb)
  void <- Map(addDataFrame, comp, sheets)
  saveWorkbook(wb, file = paste0(outwd, cluster, "_POSYESvPOSNO_5Pv1P_Test_lfcT", lfcT, "_minpct", minpct, "_", test, ".xlsx"))
  
  d1 <- full_join(comp[["POSYESvPOSNO_1P"]], comp[["POSYESvPOSNO_5P"]], by = "gene", suffix = c(".posvneg.1P", ".posvneg.5P")) %>%
    full_join(comp[["5Pv1P_Test_POSYESmPOSNO"]], by = "gene") %>%
    #filter(p_val_adj < p.th | p_val_adj.posvneg.1P < p.th | p_val_adj.posvneg.5P < p.th ) %>%
    select(gene, avg_log2FC.posvneg.1P, p_val_adj.posvneg.1P, avg_log2FC.posvneg.5P, p_val_adj.posvneg.5P, avg_log2FC, p_val_adj) %>%
    setNames(c("gene", "lfc_posvneg.1P", "p_posvneg.1P", "lfc_posvneg.5P", "p_posvneg.5P", "lfc_pos.1Pv5P", "p_pos.1Pv5P")) %>%
    mutate(cluster = cluster)
  
  DEGcount.posvneg.1P <- mutate(filter(comp[["POSYESvPOSNO_1P"]], p_val_adj < p.th), updown = as.character(ifelse(avg_log2FC > 0, "up", "down"))) %>% group_by(updown) %>% summarise(count = n()) %>% mutate(cluster = cluster, pairwise = "POSYESvPOSNO_1P")
  DEGcount.posvneg.5P <- mutate(filter(comp[["POSYESvPOSNO_5P"]], p_val_adj < p.th), updown = as.character(ifelse(avg_log2FC > 0, "up", "down"))) %>% group_by(updown) %>% summarise(count = n()) %>% mutate(cluster = cluster, pairwise = "POSYESvPOSNO_5P")
  DEGcount.pos.1Pv5P <- mutate(filter(comp[["POSYES_5Pv1P"]], p_val_adj < p.th), updown = as.character(ifelse(avg_log2FC > 0, "up", "down"))) %>% group_by(updown) %>% summarise(count = n()) %>% mutate(cluster = cluster, pairwise = "POSYES_5Pv1P")
  DEGcounts <- full_join(DEGcount.posvneg.1P, DEGcount.posvneg.5P) %>% full_join(DEGcount.pos.1Pv5P)
  
  output <- list(exported = comp, DEGs = d1, counts = DEGcounts)
  return(output)
}

# For all clusters
comps.TestPYvPN.LR <- list()

for (i in c("D1-MSNs", "D2-MSNs")) {
  comps.TestPYvPN.LR[[i]] <- allDEGs_POSYESvPOSNO_5Pv1P_Test(integrated, cluster = i, lfcT = log2(1.15), minpct = 0.3, test = "LR", p.th = .1)
  #comps.TestPYvPN.LR.nolfcT[[i]] <- allDEGs_POSvNEG_5Pv1P_HC(integrated, cluster = i, lfcT = -Inf, minpct = 0.01, test = "LR", p.th = .05)
}

comps.TestPYvPN.LR.nolfcT <- list()
allDEGs_POSYESvPOSNO_5Pv1P_Test_h1 <- function(integrated, cluster = "D1-MSNs", lfcT = .25, minpct = .25, test = "LR", p.th = .05, outwd = "/Users/arthur/Desktop/") {
  print(paste0("Working on ", cluster))
  Idents(integrated) <- "cluster_pairings_cpp_sorted_Arc"
  comp <- list()
  comp[["POSYESvPOSNO_1P"]] <- FindMarkers(integrated, ident.1 = paste0(cluster, "_1P_Test_POS_YES"), ident.2 = paste0(cluster, "_1P_Test_POS_NO"), assay = "SCT",  min.pct = minpct, logfc.threshold = lfcT, test.use = test, verbose = FALSE) %>% as.data.frame() %>% mutate(gene = rownames(.)) %>%
    mutate(p_val_adj = p.adjust(.$p_val, "fdr"))
  comp[["POSYESvPOSNO_5P"]] <- FindMarkers(integrated, ident.1 = paste0(cluster, "_5P_Test_POS_YES"), ident.2 = paste0(cluster, "_5P_Test_POS_NO"), assay = "SCT",  min.pct = minpct, logfc.threshold = lfcT, test.use = test, verbose = FALSE) %>% as.data.frame() %>% mutate(gene = rownames(.)) %>%
    mutate(p_val_adj = p.adjust(.$p_val, "fdr"))
  comp[["POSNO_5Pv1P"]] <- FindMarkers(integrated, ident.1 = paste0(cluster, "_5P_Test_POS_NO"), ident.2 = paste0(cluster, "_1P_Test_POS_NO"), assay = "SCT",  min.pct = minpct, logfc.threshold = lfcT, test.use = test, verbose = FALSE) %>% as.data.frame() %>% mutate(gene = rownames(.)) %>%
    mutate(p_val_adj = p.adjust(.$p_val, "fdr"))
  comp[["POSYES_5Pv1P"]] <- FindMarkers(integrated, ident.1 = paste0(cluster, "_5P_Test_POS_YES"), ident.2 = paste0(cluster, "_1P_Test_POS_YES"), assay = "SCT",  min.pct = minpct, logfc.threshold = lfcT, test.use = test, verbose = FALSE) %>% as.data.frame() %>% mutate(gene = rownames(.)) %>%
    mutate(p_val_adj = p.adjust(.$p_val, "fdr"))
  comp[["5Pv1P_Test_POSYESmPOSNO"]] <- full_join(comp[["POSNO_5Pv1P"]], comp[["POSYES_5Pv1P"]], by = "gene", suffix = c(".no", "")) %>%
    dplyr::filter((p_val_adj.no > p.th | is.na(p_val_adj.no)==T) & p_val_adj < p.th) %>%
    dplyr::select(!matches("no"))
  
  
  wb <- createWorkbook()                                                    
  sheetnames <- names(comp)
  sheets <- lapply(sheetnames, createSheet, wb = wb)
  void <- Map(addDataFrame, comp, sheets)
  saveWorkbook(wb, file = paste0(outwd, cluster, "_POSYESvPOSNO_5Pv1P_HConly_lfcT", lfcT, "_minpct", minpct, "_", test, ".xlsx"))
  
  d1 <- full_join(comp[["POSYESvPOSNO_1P"]], comp[["POSYESvPOSNO_5P"]], by = "gene", suffix = c(".posvneg.1P", ".posvneg.5P")) %>%
    full_join(comp[["5Pv1P_Test_POSYESmPOSNO"]], by = "gene") %>%
    #filter(p_val_adj < p.th | p_val_adj.posvneg.1P < p.th | p_val_adj.posvneg.5P < p.th ) %>%
    select(gene, avg_log2FC.posvneg.1P, p_val_adj.posvneg.1P, avg_log2FC.posvneg.5P, p_val_adj.posvneg.5P, avg_log2FC, p_val_adj) %>%
    setNames(c("gene", "lfc_posvneg.1P", "p_posvneg.1P", "lfc_posvneg.5P", "p_posvneg.5P", "lfc_pos.1Pv5P", "p_pos.1Pv5P")) %>%
    mutate(cluster = cluster)
  
  DEGcount.posvneg.1P <- mutate(filter(comp[["POSYESvPOSNO_1P"]], p_val_adj < p.th), updown = as.character(ifelse(avg_log2FC > 0, "up", "down"))) %>% group_by(updown) %>% summarise(count = n()) %>% mutate(cluster = cluster, pairwise = "POSYESvPOSNO_1P")
  DEGcount.posvneg.5P <- mutate(filter(comp[["POSYESvPOSNO_5P"]], p_val_adj < p.th), updown = as.character(ifelse(avg_log2FC > 0, "up", "down"))) %>% group_by(updown) %>% summarise(count = n()) %>% mutate(cluster = cluster, pairwise = "POSYESvPOSNO_5P")
  DEGcount.pos.1Pv5P <- mutate(filter(comp[["POSYES_5Pv1P"]], p_val_adj < p.th), updown = as.character(ifelse(avg_log2FC > 0, "up", "down"))) %>% group_by(updown) %>% summarise(count = n()) %>% mutate(cluster = cluster, pairwise = "POSYES_5Pv1P")
  DEGcounts <- full_join(DEGcount.posvneg.1P, DEGcount.posvneg.5P) %>% full_join(DEGcount.pos.1Pv5P)
  
  output <- list(exported = comp, DEGs = d1, counts = DEGcounts)
  return(output)
}
for (i in c("D1-MSNs", "D2-MSNs")) {
  #comps.TestPYvPN.LR[[i]] <- allDEGs_POSYESvPOSNO_5Pv1P_Test(integrated, cluster = i, lfcT = .25, minpct = 0.25, test = "LR", p.th = .05)
  comps.TestPYvPN.LR.nolfcT[[i]] <- allDEGs_POSYESvPOSNO_5Pv1P_Test_h1(integrated, cluster = i, lfcT = -Inf, minpct = 0.01, test = "LR", p.th = .1)
}

## to use, test only comparisons with LR testing

final <- comps.TestPYvPN.LR
counts <- lapply(final, pluck, "counts") %>% data.table::rbindlist(idcol = "Cluster") %>%
  mutate(updown = factor(updown, levels = c("up", "down")), 
         #pairwise = factor(pairwise, levels = c("POSYESvPOSNO_1P", "POSYESvPOSNO_5P", "POSYES_5Pv1P")),
         Cluster = factor(gsub("Grm8-","Grm8-\n", Cluster), levels = c("D1-MSNs", "D2-MSNs")))

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

plot.heatmap <- function(df, ttl) {
  df %>%
    dplyr::select(gene, lfc_posvneg.1P, lfc_posvneg.5P, lfc_pos.1Pv5P) %>%
    pivot_longer(cols = !gene, names_to = "comp", values_to = "log2FC") %>%
    mutate(comp = factor(case_when(comp == "lfc_posvneg.1P" ~ "1PTest - POSYESvPOSNO", comp == "lfc_posvneg.5P" ~ "5PTest - POSYESvPOSNO", comp == "lfc_pos.1Pv5P" ~ "POSYESTest - 5Pv1P"))) %>%
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

hm.d1.1P <- data.frame(gene = dplyr::filter(comps.TestPYvPN.LR$`D1-MSNs`$DEGs, p_posvneg.1P < .1)$gene) %>%
  left_join(comps.TestPYvPN.LR.nolfcT$`D1-MSNs`$DEGs) %>%
  mutate(gene = factor(gene, levels = arrange(., lfc_posvneg.1P)$gene)) %>%
  plot.heatmap(ttl = "D1 - POSYESvPOSNO_1PTest") # + theme(axis.text.x = element_text(size = 5, colour = "black", angle = 90, hjust=0.95,vjust=0.2))
hm.d1.5P <- data.frame(gene = dplyr::filter(comps.TestPYvPN.LR$`D1-MSNs`$DEGs, p_posvneg.5P < .1)$gene) %>%
  left_join(comps.TestPYvPN.LR.nolfcT$`D1-MSNs`$DEGs) %>%
  mutate(gene = factor(gene, levels = arrange(., lfc_posvneg.5P)$gene)) %>%
  plot.heatmap(ttl = "D1 - POSYESvPOSNO_5PTest") # + theme(axis.text.x = element_text(size = 5, colour = "black", angle = 90, hjust=0.95,vjust=0.2))
hm.d1.5Pv1P <- data.frame(gene = dplyr::filter(comps.TestPYvPN.LR$`D1-MSNs`$DEGs, p_pos.1Pv5P < .1)$gene) %>%
  left_join(comps.TestPYvPN.LR.nolfcT$`D1-MSNs`$DEGs) %>%
  mutate(gene = factor(gene, levels = arrange(., lfc_pos.1Pv5P)$gene)) %>%
  plot.heatmap(ttl = "D1 - POSYES_5Pv1PTest")
hm.D2.1P <- data.frame(gene = dplyr::filter(comps.TestPYvPN.LR$`D2-MSNs`$DEGs, p_posvneg.1P < .1)$gene) %>%
  left_join(comps.TestPYvPN.LR.nolfcT$`D2-MSNs`$DEGs) %>%
  mutate(gene = factor(gene, levels = arrange(., lfc_posvneg.1P)$gene)) %>%
  plot.heatmap(ttl = "D2 - POSYESvPOSNO_1PTest") # + theme(axis.text.x = element_text(size = 5, colour = "black", angle = 90, hjust=0.95,vjust=0.2))
hm.D2.5P <- data.frame(gene = dplyr::filter(comps.TestPYvPN.LR$`D2-MSNs`$DEGs, p_posvneg.5P < .1)$gene) %>%
  left_join(comps.TestPYvPN.LR.nolfcT$`D2-MSNs`$DEGs) %>%
  mutate(gene = factor(gene, levels = arrange(., lfc_posvneg.5P)$gene)) %>%
  plot.heatmap(ttl = "D2 - POSYESvPOSNO_5PTest") # + theme(axis.text.x = element_text(size = 5, colour = "black", angle = 90, hjust=0.95,vjust=0.2))
hm.D2.5Pv1P <- data.frame(gene = dplyr::filter(comps.TestPYvPN.LR$`D2-MSNs`$DEGs, p_pos.1Pv5P < .1)$gene) %>%
  left_join(comps.TestPYvPN.LR.nolfcT$`D2-MSNs`$DEGs) %>%
  mutate(gene = factor(gene, levels = arrange(., lfc_pos.1Pv5P)$gene)) %>%
  plot.heatmap(ttl = "D2 - POSYES_5Pv1PTest")
all.heatmaps <- list(hm.d1.1P, hm.d1.5P, hm.d1.5Pv1P, hm.D2.1P, hm.D2.5P, hm.D2.5Pv1P) ; rm(hm.d1.1P, hm.d1.5P, hm.d1.5Pv1P, hm.D2.1P, hm.D2.5P, hm.D2.5Pv1P)
maxgenes <- lapply(all.heatmaps, function(x) {nrow(x$data)/3}) %>% unlist() %>% max()

addextra <- function(t, maxadd) {
  toadd <- maxadd - (nrow(t$data)/3)
  if (toadd > 0) {
    tobind <- data.frame(gene = sort(rep(paste0("extra", seq(1,toadd)), 3)), comp = rep(c("1PTest - POSYESvPOSNO", "5PTest - POSYESvPOSNO", "POSYESTest - 5Pv1P"), toadd), log2FC = rep(NA, toadd*3))
    t$data <- t$data %>% rbind(tobind)
  }
  return(t)
}

all.heatmaps <- lapply(all.heatmaps, addextra, maxadd = 383)

plot_grid(plotlist = all.heatmaps, ncol = 2) %>%
  ggsave(filename = "/Users/arthur/Desktop/yellow_heatmaps.pdf", device = cairo_pdf, width = 21.6, height = 27.9, units = "cm") 


# Subtract activity dependent genes (sig YESvNO in comps.Test.YvN.LR) from DEGs

# DEGs using all reds cells, correcting for pairings & cpp
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

# For all clusters
compsYvN.Test.LR <- list()
for (i in c("D1-MSNs", "D2-MSNs")) {
  compsYvN.Test.LR[[i]] <- allDEGs_YESvNO_split_Test(integrated, cluster = i, lfcT = log2(1.15), minpct = .3, test = "LR", p.th = .1)
}

compsYvN.Test.LR.nolfcT <- list()
compsYvN.Test.LR.nolfcT[["D1-MSNs"]] <- allDEGs_YESvNO_split_Test(integrated, cluster = "D1-MSNs", lfcT = -Inf, minpct = 0.01, test = "LR", p.th = .1)
compsYvN.Test.LR.nolfcT[["D2-MSNs"]] <- allDEGs_YESvNO_split_Test(integrated, cluster = "D2-MSNs", lfcT = -Inf, minpct = 0.01, test = "LR", p.th = .1)











yvn.sig <- compsYvN.Test.LR.nolfcT %>%
  lapply(pluck, "exported") %>%
  lapply(lapply, filter, abs(avg_log2FC) > log2(1.15)) 

pyvpn.sig <- comps.TestPYvPN.LR %>% 
  lapply(pluck, "exported") %>%
  lapply(lapply, filter, p_val_adj < .1)

all.vennys <- list()
for (c in c("D1-MSNs", "D2-MSNs")) {
  for (p in c("1P", "5P")) {
    YG = pyvpn.sig[[c]][[paste0("POSYESvPOSNO_", p)]] %>% setNames(paste0(names(.), ".YG")) %>% setNames(gsub("gene.YG", "gene", names(.)))
    AD = yvn.sig[[c]][[paste0(p, "_Test_YESvNO")]] %>% setNames(paste0(names(.), ".AD")) %>% setNames(gsub("gene.AD", "gene", names(.)))
    YR = pyvpn.sig[[c]][[paste0("POSYESvNEGYES_", p)]] %>% setNames(paste0(names(.), ".YR")) %>% setNames(gsub("gene.YR", "gene", names(.)))
    all.vennys[[paste(c, p, sep = "_")]] <- list(YG = YG$gene, 
                                                 AD = AD$gene, 
                                                 YR = YR$gene) %>%
      venn(show.plot = F) %>%
      attr("intersections") %>%
      setNames(gsub(":","-", names(.))) %>%
      lapply(function(x) { data.frame(gene = x) %>% left_join(YG, by = "gene") %>% left_join(AD, by = "gene") %>% left_join(YR, by = "gene") })
    wb <- createWorkbook()                                                    
    sheetnames <- names(all.vennys[[paste(c, p, sep = "_")]])
    sheets <- lapply(sheetnames, createSheet, wb = wb)
    void <- Map(addDataFrame, all.vennys[[paste(c, p, sep = "_")]], sheets)
    saveWorkbook(wb, file = paste0(out.wd, "Venny_", c, "_",p, "_predict_react_Test", ".xlsx"))
    
  }
}

AD.all <- compsYvN.Test.LR.nolfcT %>% lapply(pluck, "exported") %>%
  lapply(rbindlist, idcol = "pairwise", use.names=T) %>%
  rbindlist(idcol = "cluster") %>%
  separate(pairwise, into = c("pairings", "cpp", "Arc"), sep = "_") %>%
  mutate(cluster_pairings = paste(cluster, pairings, sep = "_")) %>%
  dplyr::select(cluster_pairings, gene, p_val, avg_log2FC, pct.1, pct.2, p_val_adj) %>%
  setNames(c("cluster_pairings", "gene", "p_val.AD", "avg_log2FC.AD", "pct.1.AD", "pct.2.AD", "p_val_adj.AD")) #%>%
#mutate(p_val_adj.AD = NA)

# ALOS REMOVE ANY FC IN AD that is >15% even if non-sig, 
predict.genes <- all.vennys %>% lapply(pluck, "YG-YR") %>%
  rbindlist(idcol = "cluster_pairings") %>%
  select(-matches(".AD$")) %>%
  left_join(AD.all) %>%
  select(cluster_pairings, gene, matches("p_val_adj"), matches("avg_log2FC")) %>%
  
  filter(sign(avg_log2FC.YR)==sign(avg_log2FC.YG)) %>%
  pivot_longer(cols = matches("\\."), names_to = c("stat", "pairwise"), names_sep = "\\.", values_to = "value") %>%
  pivot_wider(names_from = "stat", values_from = "value") %>%
  mutate(label = case_when(p_val_adj < .1 & p_val_adj > .05 ~ "~",
                           p_val_adj < .05 & p_val_adj > .01 ~ "*",
                           p_val_adj < .01 & p_val_adj > .001 ~ "**",
                           p_val_adj < .001 & p_val_adj > .0001 ~ "***",
                           p_val_adj < .0001 ~ "****",
                           T ~ "")) %>%
  arrange(avg_log2FC) %>%
  mutate(gene = factor(gene, levels = unique(.$gene)))

sc = c(-1.52, -.8, 0, .8, 1.52); sccols = c("blue","blue","black","yellow", "yellow")

p.predict <- predict.genes %>%
  ggplot(aes(y = gene, fill = avg_log2FC, x= pairwise) ) +
  facet_wrap(~cluster_pairings, scales = "free", ncol = 1) +
  geom_tile() +
  coord_flip() +
  geom_text(aes(label = label), colour = "white", size = 5/2.8) +
  scale_x_discrete(drop = T) +
  scale_y_discrete(limits=rev, drop = F) +
  scale_fill_gradientn(colours = sccols, na.value = "grey80",
                       values = scales::rescale(sc, from = c(min(sc), max(sc))), 
                       limits = c(min(sc), max(sc)),
                       breaks = c(sc[2], 0, sc[4]),
                       guide = guide_colorbar(title  = "value", title.vjust = 0.9, frame.colour = "black", frame.linewidth = .2*2.8, ticks.linewidth = .2*2.8, ticks.colour = "black", draw.ulim = F, draw.llim = F, title.position = "top", direction = "vertical", barheight = unit(2, "char"), barwidth = unit(.5, "char"))) +
  theme_classic(base_size = 5, base_family = "Arial") +
  theme(plot.background = element_rect(fill="transparent", colour=NA),
        panel.background = element_rect(fill="transparent", colour=NA),
        #aspect.ratio = 1/6,
        plot.margin = margin(.5,.5,.5,.5, "line"),
        plot.title = element_text(size = 6, hjust = .5),
        axis.line = element_blank(), 
        axis.title.y = element_blank(),  axis.title.x = element_blank(), 
        axis.ticks = element_blank(),
        axis.text.x = element_text(size = 5, angle = 45, vjust = 0.8, hjust = 1), axis.text.y = element_text(size = 5, face = "italic"),
        strip.text = element_text(size = 5, color = "black"), strip.background = element_blank(),
        legend.position = "right", legend.justification = c("center", "center"), legend.background = element_rect(fill="transparent"), legend.box.background = element_rect(fill="transparent", color="transparent"),
        legend.text = element_text(size = 5), legend.margin = margin(0,0,0,0), legend.box.margin = margin(0,0,0,0), legend.spacing = unit(-.1, "char"), legend.title = element_text(size=5)) 
ggsave(p.predict, filename = "/Users/arthur/Desktop/react_predict_hm.pdf", device = cairo_pdf, width = 21.6, height = 27.9, units = "cm") 


# Dot plots example genes
Idents(integrated.Test) <- "comb.clusters"
tempD1 <- subset(integrated.Test, idents = "D1-MSNs")
ex.D1 <- DotPlot(tempD1, features = rev(c("Gpm6a", "Adcy9", "Map4k3", "Brd8", "Rab2a", "Mir6236" )), 
                 #cluster.idents = F,
                 assay = "SCT", dot.scale = 6, group.by = "cluster_pairings_cpp_sorted_Arc",
                 cols = c("lightgrey", "#F564E3"),
                 dot.min = .01, col.min = 0, col.max = 1, scale.min = .1, scale.max = 100) + 
  coord_flip() +
  scale_y_discrete(limits=c("D1-MSNs_1P_Test_NEG_NO", "D1-MSNs_1P_Test_NEG_YES", "D1-MSNs_1P_Test_POS_NO", "D1-MSNs_1P_Test_POS_YES",
                            "D1-MSNs_5P_Test_NEG_NO", "D1-MSNs_5P_Test_NEG_YES", "D1-MSNs_5P_Test_POS_NO", "D1-MSNs_5P_Test_POS_YES")) +
  theme_classic(base_family = "Arial", base_size = 5) +
  theme(plot.background = element_rect(fill="transparent", colour=NA),
        panel.background = element_rect(fill="transparent", colour=NA),
        #aspect.ratio = 6/4, 
        axis.title = element_blank(), axis.text = element_text(color = "black"), axis.ticks = element_line(color="black"),
        axis.text.x = element_text(face = "italic", angle = 45, vjust = 0.8, hjust = 1),
        legend.key.size = unit(.5, "line"), legend.position = "bottom", legend.box = "vertical")

tempD2 <- subset(integrated.Test, idents = "D2-MSNs")
ex.D2 <- DotPlot(tempD2, features = rev(c("Miip", "Gpm6a", "Zc3h12b","Ippk",
                                          "Sptbn4",
                                          "Slc35f3")), 
                 #cluster.idents = F,
                 assay = "SCT", dot.scale = 6, group.by = "cluster_pairings_cpp_sorted_Arc",
                 cols = c("lightgrey", "#619CFF"),
                 dot.min = .01, col.min = 0, col.max = 1) + 
  coord_flip() +
  scale_y_discrete(limits=c("D2-MSNs_1P_Test_NEG_NO", "D2-MSNs_1P_Test_NEG_YES", "D2-MSNs_1P_Test_POS_NO", "D2-MSNs_1P_Test_POS_YES",
                            "D2-MSNs_5P_Test_NEG_NO", "D2-MSNs_5P_Test_NEG_YES", "D2-MSNs_5P_Test_POS_NO", "D2-MSNs_5P_Test_POS_YES")) +
  theme_classic(base_family = "Arial", base_size = 5) +
  theme(plot.background = element_rect(fill="transparent", colour=NA),
        panel.background = element_rect(fill="transparent", colour=NA),
        #aspect.ratio = 13/4, 
        axis.title = element_blank(), axis.text = element_text(color = "black"), axis.ticks = element_line(color="black"),
        axis.text.x = element_text(face = "italic", angle = 45, vjust = 0.8, hjust = 1),
        legend.key.size = unit(.5, "line"), legend.position = "bottom", legend.box = "vertical")

ex<-plot_grid(ex.D1, ex.D2, align = "hv", axis = "tblr", ncol=2)
ggsave(ex, filename = "/Users/arthur/Desktop/ex_yellow_DotPlot.pdf", device = cairo_pdf, width = 9, height = 9, units = "cm") 






plot_gene_Test <- function(integrated, gene = "Kcnd2", cluster = "D1-MSNs", y.lim = c(0,2.5)) {
  Idents(integrated) <- "cluster_pairings_cpp"
  P1cells <- subset(integrated, idents = paste0(cluster, "_1P_Test"))
  p1 <- VlnPlot(P1cells, features = c(gene), split.by = "cluster_pairings_cpp_sorted_Arc", assay = "SCT", slot = "data", pt.size = 0) +
    scale_fill_manual(values = c("grey90","red","green2","yellow")) +
    scale_y_continuous(limits = y.lim, expand = c(0,0)) +
    xlab(cluster) +
    ylab("Normalized gene expression") +
    theme_classic(base_size = 5, base_line_size = .2, base_family = "Arial") +
    theme(legend.position = "bottom", legend.key.size = unit(.1, "line"), legend.direction = "vertical",
          aspect.ratio = 1, 
          axis.text = element_text(color = "black"), axis.ticks = element_line(color="black"),
          axis.text.x = element_blank(), axis.ticks.x = element_blank())
  P2cells <- subset(integrated, idents = paste0(cluster, "_5P_Test"))
  p2 <- VlnPlot(P2cells, features = c(gene), split.by = "cluster_pairings_cpp_sorted_Arc", assay = "SCT", slot = "data", pt.size = 0) +
    scale_fill_manual(values = c("grey90","red","green2","yellow")) +
    scale_y_continuous(limits = y.lim, expand = c(0,0)) +
    xlab(cluster) +
    ylab("Normalized gene expression") +
    theme_classic(base_size = 5, base_line_size = .2, base_family = "Arial") +
    theme(legend.position = "bottom", legend.key.size = unit(.1, "line"), legend.direction = "vertical",
          aspect.ratio = 1, 
          axis.text = element_text(color = "black"), axis.ticks = element_line(color="black"),
          axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p <- plot_grid(p1,p2, labels = c("1P", "5P"), ncol = 2)
  return(p)
}

cacna1e <- plot_gene_Test(integrated, "Cacna1e","D1-MSNs")

# final figure
plot_grid(umap.Test.yellow, plot_grid(p.act.Test, p.act.HC, nrow=1, align="hv", axis="tblr"), p.predict, plot_grid(cacna1e, p.counts, nrow=2, align = "hv", axis = "tblr"), nrow = 2, align = "hv", axis = "trbl") %>%
  ggsave(filename = "/Users/arthur/Desktop/yellow.pdf", device = cairo_pdf, width = 21.6, height = 27.9, units = "cm") 


