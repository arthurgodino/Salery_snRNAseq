## Reviewers analysis

# Comparison with Philipps
integrated <- readRDS(file = "/Users/agodino/Desktop/Marine_integrated.rds")

integrated_2 <- integrated
DefaultAssay(integrated_2) <- "integrated"
integrated_2 <- FindNeighbors(integrated_2, dims = 1:16, verbose = TRUE)
integrated_2 <- FindClusters(integrated_2, resolution = 0.2, graph.name = 'integrated_snn', verbose = TRUE)
integrated_2 <- RunUMAP(integrated_2, dims = 1:16, verbose = TRUE)
integrated_2$byPhillips <- case_when(integrated_2$seurat_clusters == 1 ~ 'Drd1-MSN-1',
                                     integrated_2$seurat_clusters == 2 ~ 'Drd1-MSN-1',
                                     integrated_2$seurat_clusters == 5 ~ 'Drd1-MSN-1',
                                     integrated_2$seurat_clusters == 4 ~ 'Drd1-MSN-2',
                                     integrated_2$seurat_clusters == 7 ~ 'Grm8-MSN',
                                     T ~ NA_character_)
Idents(integrated_2) <- "byPhillips"
UMAPPlot(integrated_2, label = T) + NoLegend() & theme(aspect.ratio = 1)

integrated$byPhillips <- integrated_2$byPhillips
Idents(integrated) <- "byPhillips"
UMAPPlot(integrated, label = T) + NoLegend() & theme(aspect.ratio = 1)

D1sub <- subset(integrated, idents = c("Drd1-MSN-1", "Drd1-MSN-2", "Grm8-MSN")) %>%
  DotPlot(features = rev(c("Grm8", "Slc24a2", "Cacng3", "Trpc3", "Itpr2", "Cnr1", "Syt6", "Syt10",  "Cpne8", "Syt4", "Prkca", "Prkcg", "Prkch", "Kcnc1", "Kcng3", "Kcnc2", "Kcnk2", "Kcnd3", "Kcnt2", "Kcnh5", "Kcnn2", "Kcnk13", "Kcnj3", "Kcnip1")), 
                      cluster.idents = F,
                      group.by = "byPhillips", assay = "SCT", dot.scale = 6,
                      cols = c("lightgrey","red"),
                      dot.min = .01, col.min = 0, col.max = 5, scale=F) + 
  scale_y_discrete(limits=rev) +
  coord_flip() +
  theme_classic(base_family = "Arial", base_size = 5) +
  theme(plot.background = element_rect(fill="transparent", colour=NA),
        panel.background = element_rect(fill="transparent", colour=NA),
        aspect.ratio = 20/7, axis.title = element_blank(), axis.text = element_text(color = "black"), axis.ticks = element_line(color="black"),
        axis.text.x = element_text(face = "italic", angle = 45, vjust = 0.8, hjust = 1),
        axis.text.y = element_text(face = "italic"),
        legend.key.size = unit(.5, "line"), legend.position = "bottom", legend.box = "vertical")
ggsave(D1sub, filename = "/Users/agodino/Desktop/philipsgenes.pdf", device = cairo_pdf, width = 6, height = 12, units = "cm") 

D1ieg <- subset(integrated, idents = c("Drd1-MSN-1", "Drd1-MSN-2", "Grm8-MSN")) %>%
  DotPlot(features = rev(c("Arc", "Fos", "Fosb", "Nr4a3", "Fosl2", "Tiparp")), 
          cluster.idents = F,
          split.by = "cpp",
          group.by = "byPhillips", assay = "SCT", dot.scale = 6,
          cols = c("lightgrey","red"),
          dot.min = .01, col.min = 0, col.max = 5, scale=F) + 
  scale_y_discrete(limits=rev) +
  coord_flip() +
  theme_classic(base_family = "Arial", base_size = 5) +
  theme(plot.background = element_rect(fill="transparent", colour=NA),
        panel.background = element_rect(fill="transparent", colour=NA),
        aspect.ratio = 20/7, axis.title = element_blank(), axis.text = element_text(color = "black"), axis.ticks = element_line(color="black"),
        axis.text.x = element_text(face = "italic", angle = 45, vjust = 0.8, hjust = 1),
        axis.text.y = element_text(face = "italic"),
        legend.key.size = unit(.5, "line"), legend.position = "bottom", legend.box = "vertical")
ggsave(D1sub, filename = "/Users/agodino/Desktop/philipsgenes.pdf", device = cairo_pdf, width = 6, height = 12, units = "cm") 

integrated$phil_cpp <- paste(integrated$byPhillips, integrated$cpp, sep="_")
Idents(integrated) <- "phil_cpp"
testvsHC <- list()
for (i in c("Drd1-MSN-1", "Drd1-MSN-2", "Grm8-MSN")) {
  print(i)
testvsHC[[i]] <- FindMarkers(integrated, ident.1 = paste0(i, "_Test"), ident.2 = paste0(i, "_HC"), assay = "SCT",  min.pct = .01, logfc.threshold = -Inf, latent.vars = c("sorted", "pairings"), test.use = "LR", verbose = FALSE) %>% as.data.frame() %>% mutate(gene = rownames(.)) %>%
  mutate(p_val_adj = p.adjust(.$p_val, "fdr"))
}

betweenphilclusts <- list()
  betweenphilclusts[["D1-1vD1-2"]] <- FindMarkers(integrated, ident.1 = "Drd1-MSN-1_Test", ident.2 = "Drd1-MSN-2_Test", assay = "SCT",  min.pct = .01, logfc.threshold = -Inf, latent.vars = c("sorted", "pairings"), test.use = "LR", verbose = FALSE) %>% as.data.frame() %>% mutate(gene = rownames(.)) %>%
  mutate(p_val_adj = p.adjust(.$p_val, "fdr"))
  betweenphilclusts[["D1-1vDGrm"]] <- FindMarkers(integrated, ident.1 = "Drd1-MSN-1_Test", ident.2 = "Grm8-MSN_Test", assay = "SCT",  min.pct = .01, logfc.threshold = -Inf, latent.vars = c("sorted", "pairings"), test.use = "LR", verbose = FALSE) %>% as.data.frame() %>% mutate(gene = rownames(.)) %>%
    mutate(p_val_adj = p.adjust(.$p_val, "fdr"))
  betweenphilclusts[["D1-2vGrm"]] <- FindMarkers(integrated, ident.1 = "Drd1-MSN-2_Test", ident.2 = "Grm8-MSN_Test", assay = "SCT",  min.pct = .01, logfc.threshold = -Inf, latent.vars = c("sorted", "pairings"), test.use = "LR", verbose = FALSE) %>% as.data.frame() %>% mutate(gene = rownames(.)) %>%
    mutate(p_val_adj = p.adjust(.$p_val, "fdr"))

write.csv(testvsHC$`Drd1-MSN-1`, file = "/Users/agodino/Desktop/testvHC_D1-1.csv")
write.csv(testvsHC$`Drd1-MSN-2`, file = "/Users/agodino/Desktop/testvHC_D1-2.csv")
write.csv(testvsHC$`Grm8-MSN`, file = "/Users/agodino/Desktop/testvHC_Grm8.csv")

write.csv(betweenphilclusts$`D1-1vD1-2`, file = "/Users/agodino/Desktop/1v2.csv")
write.csv(betweenphilclusts$`D1-1vDGrm`, file = "/Users/agodino/Desktop/1vG.csv")
write.csv(betweenphilclusts$`D1-2vGrm`, file = "/Users/agodino/Desktop/2vG.csv")

  

Idents(integrated) <- "cpp"
test <- subset(integrated, idents = "Test")
hc <- subset(integrated, idents = "HC")
Idents(test) <- "byPhillips"
Idents(hc) <- "byPhillips"

matrix.Test.1 <- as.matrix(subset(test, idents = "Drd1-MSN-1")@assays$SCT@data) %>% t() %>% as.data.frame() %>% select(Arc, Fos, Fosb, Fosl2, Nr4a3, Tiparp) %>% mutate(cluster = "Drd1-MSN-1", cpp = "Test")
matrix.Test.2 <- as.matrix(subset(test, idents = "Drd1-MSN-2")@assays$SCT@data) %>% t() %>% as.data.frame() %>% select(Arc, Fos, Fosb,Fosl2, Nr4a3, Tiparp) %>% mutate(cluster = "Drd1-MSN-2", cpp = "Test")
matrix.Test.8 <- as.matrix(subset(test, idents = "Grm8-MSN")@assays$SCT@data) %>% t() %>% as.data.frame() %>% select(Arc, Fos, Fosb,Fosl2, Nr4a3, Tiparp) %>% mutate(cluster = "Grm8-MSN", cpp = "Test")
matrix.HC.1 <- as.matrix(subset(hc, idents = "Drd1-MSN-1")@assays$SCT@data) %>% t() %>% as.data.frame() %>% select(Arc, Fos, Fosb,Fosl2, Nr4a3, Tiparp) %>% mutate(cluster = "Drd1-MSN-1", cpp = "HC")
matrix.HC.2 <- as.matrix(subset(hc, idents = "Drd1-MSN-2")@assays$SCT@data) %>% t() %>% as.data.frame() %>% select(Arc, Fos, Fosb,Fosl2, Nr4a3, Tiparp) %>% mutate(cluster = "Drd1-MSN-2", cpp = "HC")
matrix.HC.8 <- as.matrix(subset(hc, idents = "Grm8-MSN")@assays$SCT@data) %>% t() %>% as.data.frame() %>% select(Arc, Fos, Fosb,Fosl2, Nr4a3, Tiparp)%>% mutate(cluster = "Grm8-MSN", cpp = "HC")

comb <- rbindlist(list(matrix.Test.1, matrix.Test.2, matrix.Test.8, matrix.HC.1, matrix.HC.2, matrix.HC.8))

msd <- comb %>% 
       pivot_longer(cols = c("Arc", "Fos", "Fosb", "Fosl2", "Nr4a3", "Tiparp"), names_to = "gene", values_to = "expr") %>%
        group_by(gene, cluster, cpp) %>%
                summarize(sem = sd(expr)/sqrt(n()), expr = mean(expr))

plotIEG <- msd %>%
           ggplot(aes(x=cluster, y=expr, fill=cpp)) +
           facet_wrap(~gene, nrow =2, scales = "free") +
           geom_col(position = position_dodge()) +
           geom_errorbar(aes(ymin=expr-sem, ymax=expr+sem), position = position_dodge() ) +
           theme_classic() +
           theme(aspect.ratio = 1)
ggsave(plotIEG, filename = "/Users/agodino/Desktop/CPP_IEGs.pdf", device = cairo_pdf, width = 12, height = 12, units = "cm") 

for (g in c("Arc", "Fos", "Fosb", "Fosl2", "Nr4a3", "Tiparp")) {
  tmp <- comb %>% select(!!sym(g), cluster, cpp) %>% mutate(expr=!!sym(g)) 
  model <- lm(expr ~ cluster*cpp, data=tmp)
  anova <- car::Anova(model, type="III")
  capture.output(anova, file=paste0("/Users/agodino/Desktop/", g, "_anova.txt"))
  emm <- emmeans::emmeans(model, c("cluster", "cpp")) %>% pairs(simple="each")
  ph <- rbind(emm[[1]], emm[[2]], adjust="fdr")
  capture.output(ph, file=paste0("/Users/agodino/Desktop/", g, "_fdr.txt"))
}






## Tagged cells

GFP_Phil <- table(integrated$byPhillips, integrated$sorted_pairings) %>%
    t() %>% 
    as.data.frame() %>%
    group_by(Var1) %>%
    mutate(Total = sum(Freq)) %>%
    ungroup() %>%
    mutate(perc = 100*Freq/Total)

p.AP <- GFP_Phil %>%
                ggplot(aes(x=Var2, y=perc, color=Var1)) +
                geom_line(aes(group=Var1)) +
                geom_point() +
                scale_y_continuous(limits=c(0,100), expand=c(0,0)) +
                scale_color_manual(values=c("grey85","grey75","green2","green4")) +
                theme_classic() +
                theme(aspect.ratio=1)
ggsave(p.AP, filename = "/Users/agodino/Desktop/ArcGFP_byPhillips.pdf", device = cairo_pdf, width = 12, height = 12, units = "cm") 


cont1P <- GFP_Phil %>%
          filter(grepl("1P", Var1)==T) %>%
          mutate(sorted = gsub("_1P", "", Var1)) %>%
          select(sorted, Var2, Freq) %>%
          pivot_wider(names_from = Var2, values_from = Freq) %>%
          column_to_rownames("sorted") %>% as.matrix() %>% as.table()
X2.c1P <- chisq.test(cont1P)
capture.output(X2.c1P, file = "/Users/agodino/Desktop/X2.c1P.txt")

cont5P <- GFP_Phil %>%
  filter(grepl("5P", Var1)==T) %>%
  mutate(sorted = gsub("_5P", "", Var1)) %>%
  select(sorted, Var2, Freq) %>%
  pivot_wider(names_from = Var2, values_from = Freq) %>%
  column_to_rownames("sorted") %>% as.matrix() %>% as.table()
X2.c5P <- chisq.test(cont5P)
capture.output(X2.c5P, file = "/Users/agodino/Desktop/X2.c5P.txt")


pair.p.G <- data.frame(p.1P = 2*pnorm(colMeans(abs(X2.c1P$residuals)), mean=0, sd=1, lower.tail=F),
                       p.5P = 2*pnorm(colMeans(abs(X2.c5P$residuals)), mean=0, sd=1, lower.tail=F)) %>%
  tibble::rownames_to_column(var="cluster") %>%
  pivot_longer(cols = starts_with("p"), names_to = "pairing", values_to = "p") %>%
  mutate(p.adj = p.adjust(p, method = "fdr")) %>%
  mutate(pairing = gsub("p.","",pairing)) %>%
  mutate(sig = case_when(p.adj < 0.001 ~ "***", 
                         p.adj < 0.01 & p.adj > 0.001 ~ "**",
                         p.adj < 0.05 & p.adj > 0.01 ~ "*",
                         T ~ ""))
write.csv(pair.p.G, "/Users/agodino/Desktop/pairwise.Green.csv")

#####
#####

p1 <- FeaturePlot(integrated, reduction = "umap", features = c("nCount_RNA")) & theme(aspect.ratio = 1)
p2 <- FeaturePlot(integrated, reduction = "umap", features = c("nCount_SCT")) & theme(aspect.ratio = 1)
p3 <- FeaturePlot(integrated, reduction = "umap", features = c("nFeature_RNA")) & theme(aspect.ratio = 1)
p4 <- FeaturePlot(integrated, reduction = "umap", features = c("nFeature_SCT")) & theme(aspect.ratio = 1)
p5 <- FeaturePlot(integrated_2, reduction = "umap", features = c("Grm8"), min.cutoff = 0, max.cutoff = "q90") & theme(aspect.ratio = 1)
FeaturePlot(integrated, reduction = "umap", features = c("Slc24a2"), min.cutoff = 0, max.cutoff = "q90") & theme(aspect.ratio = 1)
FeaturePlot(integrated, reduction = "umap", features = c("Cacng3"), min.cutoff = 0, max.cutoff = "q90") & theme(aspect.ratio = 1)
FeaturePlot(integrated, reduction = "umap", features = c("Trpc3"), min.cutoff = 0, max.cutoff = "q90") & theme(aspect.ratio = 1)
FeaturePlot(integrated, reduction = "umap", features = c("Itpr2"), min.cutoff = 0, max.cutoff = "q90") & theme(aspect.ratio = 1)
FeaturePlot(integrated, reduction = "umap", features = c("Cnr1"), min.cutoff = 0, max.cutoff = "q90") & theme(aspect.ratio = 1)
FeaturePlot(integrated, reduction = "umap", features = c("Syt6"), min.cutoff = 0, max.cutoff = "q90") & theme(aspect.ratio = 1)
FeaturePlot(integrated, reduction = "umap", features = c("Syt10"), min.cutoff = 0, max.cutoff = "q90") & theme(aspect.ratio = 1)
FeaturePlot(integrated, reduction = "umap", features = c("Cpne8"), min.cutoff = 0, max.cutoff = "q90") & theme(aspect.ratio = 1)
FeaturePlot(integrated, reduction = "umap", features = c("Syt4"), min.cutoff = 0, max.cutoff = "q90") & theme(aspect.ratio = 1)
FeaturePlot(integrated, reduction = "umap", features = c("Prkca"), min.cutoff = 0, max.cutoff = "q90") & theme(aspect.ratio = 1)
FeaturePlot(integrated, reduction = "umap", features = c("Prkcg"), min.cutoff = 0, max.cutoff = "q90") & theme(aspect.ratio = 1)
FeaturePlot(integrated, reduction = "umap", features = c("Prkch"), min.cutoff = 0, max.cutoff = "q90") & theme(aspect.ratio = 1)
FeaturePlot(integrated, reduction = "umap", features = c("Kcnc1"), min.cutoff = 0, max.cutoff = "q90") & theme(aspect.ratio = 1)
FeaturePlot(integrated, reduction = "umap", features = c("Kcng3"), min.cutoff = 0, max.cutoff = "q90") & theme(aspect.ratio = 1)
FeaturePlot(integrated, reduction = "umap", features = c("Kcnk2"), min.cutoff = 0, max.cutoff = "q90") & theme(aspect.ratio = 1)
FeaturePlot(integrated, reduction = "umap", features = c("Kcnd3"), min.cutoff = 0, max.cutoff = "q90") & theme(aspect.ratio = 1)
FeaturePlot(integrated, reduction = "umap", features = c("Kcnt2"), min.cutoff = 0, max.cutoff = "q90") & theme(aspect.ratio = 1)
p6 <- FeaturePlot(integrated, reduction = "umap", features = c("Kcnj3"), min.cutoff = 0, max.cutoff = "q90") & theme(aspect.ratio = 1)
ggsave(cowplot::plot_grid(p1,p2,p3,p4,p5,p6), filename = "/Users/agodino/Desktop/extraumaps.pdf", device = cairo_pdf, width = 50, height = 50, units = "cm") 


############ GREEN #######################
## Background lists GREEN
integrated <- readRDS(file = "/Users/agodino/Desktop/Marine_integrated.rds")

dir.create("/Users/agodino/Desktop/forIPA_wBG")
setwd("/Users/agodino/Desktop/forIPA_wBG")
library(readxl)
library(tidyverse)

D1 <- subset(integrated, idents = "D1-MSNs") %>%
      GetAssayData(assay="SCT", layer="counts") %>%
      as.data.frame() %>%
      mutate_all(~ ifelse(. != 0, 1, 0))
ex.per <- 100*Matrix::rowSums(D1)/ncol(D1)
D1.BG <- data.frame(gene = names(ex.per), exp.percent.D1 = ex.per)

D2 <- subset(integrated, idents = "D2-MSNs") %>%
  GetAssayData(assay="SCT", layer="counts") %>%
  as.data.frame() %>%
  mutate_all(~ ifelse(. != 0, 1, 0))
ex.per <- 100*Matrix::rowSums(D2)/ncol(D2)
D2.BG <- data.frame(gene = names(ex.per), exp.percent.D2 = ex.per)

BG <- left_join(D1.BG, D2.BG)

D1.1P.PvN <- read_excel("/Users/agodino/Desktop/inputs/D1_1P_POSvNEG.xlsx") %>%
              select(gene, pct.1, pct.2, avg_log2FC, p_val_adj) %>%
              setNames(c("gene", "pct.1.D1.1P.PvN", "pct.2.D1.1P.PvN", "log2FC.D1.1P.PvN", "padj.D1.1P.PvN"))
D1.5P.PvN <- read_excel("/Users/agodino/Desktop/inputs/D1_5P_POSvNEG.xlsx") %>%
  select(gene, pct.1, pct.2, avg_log2FC, p_val_adj) %>%
  setNames(c("gene", "pct.1.D1.5P.PvN", "pct.2.D1.5P.PvN", "log2FC.D1.5P.PvN", "padj.D1.5P.PvN"))
D1.5v1.P <- read_excel("/Users/agodino/Desktop/inputs/D1_POS_5Pv1P.xlsx") %>%
  select(gene, pct.1, pct.2, avg_log2FC, p_val_adj) %>%
  setNames(c("gene", "pct.1.D1.5v1.P", "pct.2.D1.5v1.P", "log2FC.D1.5v1.P", "padj.D1.5v1.P"))

D2.1P.PvN <- read_excel("/Users/agodino/Desktop/inputs/D2_1P_POSvNEG.xlsx") %>%
  select(gene, pct.1, pct.2, avg_log2FC, p_val_adj) %>%
  setNames(c("gene", "pct.1.D2.1P.PvN", "pct.2.D2.1P.PvN", "log2FC.D2.1P.PvN", "padj.D2.1P.PvN"))
D2.5P.PvN <- read_excel("/Users/agodino/Desktop/inputs/D2_5P_POSvNEG.xlsx") %>%
  select(gene, pct.1, pct.2, avg_log2FC, p_val_adj) %>%
  setNames(c("gene", "pct.1.D2.5P.PvN", "pct.2.D2.5P.PvN", "log2FC.D2.5P.PvN", "padj.D2.5P.PvN"))
D2.5v1.P <- read_excel("/Users/agodino/Desktop/inputs/D2_POS_5Pv1P.xlsx") %>%
  select(gene, pct.1, pct.2, avg_log2FC, p_val_adj) %>%
  setNames(c("gene", "pct.1.D2.5v1.P", "pct.2.D2.5v1.P", "log2FC.D2.5v1.P", "padj.D2.5v1.P"))

all <- BG %>%
       left_join(D1.1P.PvN) %>%
  left_join(D1.5P.PvN) %>%
  left_join(D1.5v1.P) %>%
  left_join(D2.1P.PvN) %>%
  left_join(D2.5P.PvN) %>%
  left_join(D2.5v1.P) %>%
  mutate(across(contains("log2FC"), ~ ifelse(is.na(.), 0, .))) %>%
  mutate(across(contains("padj"), ~ ifelse(is.na(.), 1, .)))
#write.csv(all, "all_w_BG.csv")  

allD1 <- all %>% filter(exp.percent.D1>0) %>%
                 select(matches("gene|D1")) %>%
  select(!matches("pct"))
allD2 <- all %>% filter(exp.percent.D2>0) %>%
  select(matches("gene|D2")) %>%
  select(!matches("pct"))
write.csv(allD1, file="allD1.csv", row.names=F)
write.csv(allD2, file="allD2.csv", row.names=F)


for (g in c("1P.PvN", "5P.PvN", "5v1.P")) {
allD1 %>%
          select(matches(paste0("gene|",g))) %>%
          write.csv(file = paste0("D1_", g, ".csv"), row.names=F)
allD2 %>%
    select(matches(paste0("gene|",g))) %>%
    write.csv(file = paste0("D2_", g, ".csv"), row.names=F)
}


## IPA upstream GREEN
setwd("/Users/agodino/Desktop/IPA_results/")
allUpsReg <- lapply(list.files(pattern="UpRegs"), read_excel, skip=1)
names(allUpsReg) <- list.files(pattern="UpRegs")
df.UpsReg <- data.table::rbindlist(allUpsReg, idcol = "File", fill=T) %>%
  mutate(comp = gsub(".xls", "", File), File=NULL) %>%
  mutate(comp = gsub("UpRegs_green_", "", comp)) %>%
  dplyr::select(comp, 1, 3,5, 6, 8,9, 10) %>%
  setNames(c("comp", "UpsReg", "type", "biascorr.zscore", "zscore", "pval", "pval_adj", "targets")) %>%
  mutate(ngenes = 1+str_count(targets, ",")) %>%
  pivot_wider(id_cols = c("UpsReg", "type"), names_sep = "__", names_from = "comp", values_from = c("biascorr.zscore", "zscore", "pval", "pval_adj", "targets", "ngenes")) %>%
  mutate(n.NA = rowSums(is.na(.[,15:20]))) %>%
  mutate(ngenes_D1_1P.PvN = ifelse(is.na(ngenes__D1_1P.PvN)==T, 0, ngenes__D1_1P.PvN), 
         ngenes_D1_5P.PvN = ifelse(is.na(ngenes__D1_5P.PvN)==T, 0, ngenes__D1_5P.PvN),
         ngenes_D1_5v1.P = ifelse(is.na(ngenes__D1_5v1.P)==T, 0, ngenes__D1_5v1.P),
         ngenes_D2_1P.PvN = ifelse(is.na(ngenes__D2_1P.PvN)==T, 0, ngenes__D2_1P.PvN), 
         ngenes_D2_5P.PvN = ifelse(is.na(ngenes__D2_5P.PvN)==T, 0, ngenes__D2_5P.PvN),
         ngenes_D2_5v1.P = ifelse(is.na(ngenes__D2_5v1.P)==T, 0, ngenes__D2_5v1.P)) %>%
  rowwise() %>%
  mutate(ordersig = sum(pval_adj__D1_1P.PvN, pval_adj__D1_5P.PvN, pval_adj__D1_5v1.P, pval_adj__D2_1P.PvN, pval_adj__D2_5P.PvN, pval_adj__D2_5v1.P, na.rm =T),
         avg.ngenes = sum(ngenes__D1_1P.PvN, ngenes__D1_5P.PvN, ngenes__D1_5v1.P, ngenes__D2_1P.PvN, ngenes__D2_5P.PvN, ngenes__D2_5v1.P, na.rm =T)/6) %>%
  arrange(n.NA, ordersig) %>%
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

cowplot::plot_grid(p.chems, p.genes, align="hv", axis="tblr", nrow=2) %>%
  ggsave(filename = "green_UpsReg.pdf", device = cairo_pdf, width = 15, height = 15, units = "cm")







############ RED#######################
## Background lists RED
integrated <- readRDS(file = "/Users/agodino/Desktop/inputs/Marine_integrated.rds")

dir.create("/Users/agodino/Desktop/forIPA_wBG")
setwd("/Users/agodino/Desktop/forIPA_wBG")
library(readxl)
library(tidyverse)

D1 <- subset(integrated, idents = "D1-MSNs") %>%
  GetAssayData(assay="SCT", layer="counts") %>%
  as.data.frame() %>%
  mutate_all(~ ifelse(. != 0, 1, 0))
ex.per <- 100*Matrix::rowSums(D1)/ncol(D1)
D1.BG <- data.frame(gene = names(ex.per), exp.percent.D1 = ex.per)

D2 <- subset(integrated, idents = "D2-MSNs") %>%
  GetAssayData(assay="SCT", layer="counts") %>%
  as.data.frame() %>%
  mutate_all(~ ifelse(. != 0, 1, 0))
ex.per <- 100*Matrix::rowSums(D2)/ncol(D2)
D2.BG <- data.frame(gene = names(ex.per), exp.percent.D2 = ex.per)

BG <- left_join(D1.BG, D2.BG)

D1.1P.YvN <- read.csv("/Users/agodino/Desktop/inputs/D1_1P_Test_YESvNO.csv") %>%
  select(gene, pct.1, pct.2, avg_log2FC, p_val_adj) %>%
  setNames(c("gene", "pct.1.D1.1P.YvN", "pct.2.D1.1P.YvN", "log2FC.D1.1P.YvN", "padj.D1.1P.YvN"))
D1.5P.YvN <- read.csv("/Users/agodino/Desktop/inputs/D1_5P_Test_YESvNO.csv") %>%
  select(gene, pct.1, pct.2, avg_log2FC, p_val_adj) %>%
  setNames(c("gene", "pct.1.D1.5P.YvN", "pct.2.D1.5P.YvN", "log2FC.D1.5P.YvN", "padj.D1.5P.YvN"))
D1.5v1.Y <- read.csv("/Users/agodino/Desktop/inputs/D1_5Pv1P_Test_YES.csv") %>%
  select(gene, pct.1, pct.2, avg_log2FC, p_val_adj) %>%
  setNames(c("gene", "pct.1.D1.5v1.Y", "pct.2.D1.5v1.Y", "log2FC.D1.5v1.Y", "padj.D1.5v1.Y"))

D2.1P.YvN <- read.csv("/Users/agodino/Desktop/inputs/D2_1P_Test_YESvNO.csv") %>%
  select(gene, pct.1, pct.2, avg_log2FC, p_val_adj) %>%
  setNames(c("gene", "pct.1.D2.1P.YvN", "pct.2.D2.1P.YvN", "log2FC.D2.1P.YvN", "padj.D2.1P.YvN"))
D2.5P.YvN <- read.csv("/Users/agodino/Desktop/inputs/D2_5P_Test_YESvNO.csv") %>%
  select(gene, pct.1, pct.2, avg_log2FC, p_val_adj) %>%
  setNames(c("gene", "pct.1.D2.5P.YvN", "pct.2.D2.5P.YvN", "log2FC.D2.5P.YvN", "padj.D2.5P.YvN"))
D2.5v1.Y <- read.csv("/Users/agodino/Desktop/inputs/D2_5Pv1P_Test_YES.csv") %>%
  select(gene, pct.1, pct.2, avg_log2FC, p_val_adj) %>%
  setNames(c("gene", "pct.1.D2.5v1.Y", "pct.2.D2.5v1.Y", "log2FC.D2.5v1.Y", "padj.D2.5v1.Y"))

all <- BG %>%
  left_join(D1.1P.YvN) %>%
  left_join(D1.5P.YvN) %>%
  left_join(D1.5v1.Y) %>%
  left_join(D2.1P.YvN) %>%
  left_join(D2.5P.YvN) %>%
  left_join(D2.5v1.Y) %>%
  mutate(across(contains("log2FC"), ~ ifelse(is.na(.), 0, .))) %>%
  mutate(across(contains("padj"), ~ ifelse(is.na(.), 1, .)))
#write.csv(all, "all_w_BG.csv")  

allD1 <- all %>% filter(exp.percent.D1>0) %>%
  select(matches("gene|D1")) %>%
  select(!matches("pct"))
allD2 <- all %>% filter(exp.percent.D2>0) %>%
  select(matches("gene|D2")) %>%
  select(!matches("pct"))
write.csv(allD1, file="red_D1.csv", row.names=F)
write.csv(allD2, file="red_D2.csv", row.names=F)


for (g in c("1P.YvN", "5P.YvN", "5v1.Y")) {
  allD1 %>%
    select(matches(paste0("gene|",g))) %>%
    write.csv(file = paste0("D1_", g, ".csv"), row.names=F)
  allD2 %>%
    select(matches(paste0("gene|",g))) %>%
    write.csv(file = paste0("D2_", g, ".csv"), row.names=F)
}


## IPA upstream RED
setwd("/Users/agodino/Desktop/IPA_results/")
allUpsReg <- lapply(list.files(pattern="UpRegs_red"), read_excel, skip=1)
names(allUpsReg) <- list.files(pattern="UpRegs_red")
df.UpsReg <- data.table::rbindlist(allUpsReg, idcol = "File", fill=T) %>%
  mutate(comp = gsub(".xls", "", File), File=NULL) %>%
  mutate(comp = gsub("UpRegs_red_", "", comp)) %>%
  dplyr::select(comp, 1, 3,5, 6, 8,9, 10) %>%
  setNames(c("comp", "UpsReg", "type", "biascorr.zscore", "zscore", "pval", "pval_adj", "targets")) %>%
  mutate(ngenes = 1+str_count(targets, ",")) %>%
  pivot_wider(id_cols = c("UpsReg", "type"), names_sep = "__", names_from = "comp", values_from = c("biascorr.zscore", "zscore", "pval", "pval_adj", "targets", "ngenes")) %>%
  mutate(n.NA = rowSums(is.na(.[,15:20]))) %>%
  mutate(ngenes_D1_1P.YvN = ifelse(is.na(ngenes__D1_1P.YvN)==T, 0, ngenes__D1_1P.YvN), 
         ngenes_D1_5P.YvN = ifelse(is.na(ngenes__D1_5P.YvN)==T, 0, ngenes__D1_5P.YvN),
         ngenes_D1_5v1.Y = ifelse(is.na(ngenes__D1_5v1.Y)==T, 0, ngenes__D1_5v1.Y),
         ngenes_D2_1P.YvN = ifelse(is.na(ngenes__D2_1P.YvN)==T, 0, ngenes__D2_1P.YvN), 
         ngenes_D2_5P.YvN = ifelse(is.na(ngenes__D2_5P.YvN)==T, 0, ngenes__D2_5P.YvN),
         ngenes_D2_5v1.Y = ifelse(is.na(ngenes__D2_5v1.Y)==T, 0, ngenes__D2_5v1.Y)) %>%
  rowwise() %>%
  mutate(ordersig = sum(pval_adj__D1_1P.YvN, pval_adj__D1_5P.YvN, pval_adj__D1_5v1.Y, pval_adj__D2_1P.YvN, pval_adj__D2_5P.YvN, pval_adj__D2_5v1.Y, na.rm =T),
         avg.ngenes = sum(ngenes__D1_1P.YvN, ngenes__D1_5P.YvN, ngenes__D1_5v1.Y, ngenes__D2_1P.YvN, ngenes__D2_5P.YvN, ngenes__D2_5v1.Y, na.rm =T)/6) %>%
  arrange(n.NA, ordersig) %>%
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
                        guide = guide_colorbar(title  = "Activation z-score", title.vjust = 0.9, frame.colour = "black", frame.linewidth = .2*2.8, ticks.linewidth = .2*2.8, ticks.colour = "black", draw.ulim = F, draw.llim = F, title.Yosition = "top", direction = "vertical", barheight = unit(2, "char"), barwidth = unit(.5, "char"))) +
  theme_classic(base_size = 6, base_family = "Arial") +
  theme(plot.background = element_rect(fill="transparent", colour=NA),
        panel.background = element_rect(fill="transparent", colour=NA),
        #aspect.ratio = 1/1,
        plot.margin = margin(.5,.5,.5,.5, "line"),
        axis.title.y = element_blank(),  axis.title.x = element_blank(), 
        axis.text = element_text(size = 5, colour = "black"),
        strip.text = element_text(size = 5, color = "black"), strip.background = element_blank(),
        legend.Yosition = "right", legend.justification = c("center", "center"), legend.background = element_rect(fill="transparent"), legend.box.background = element_rect(fill="transparent", color="transparent"),
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
                        guide = guide_colorbar(title  = "Activation z-score", title.vjust = 0.9, frame.colour = "black", frame.linewidth = .2*2.8, ticks.linewidth = .2*2.8, ticks.colour = "black", draw.ulim = F, draw.llim = F, title.Yosition = "top", direction = "vertical", barheight = unit(2, "char"), barwidth = unit(.5, "char"))) +
  theme_classic(base_size = 6, base_family = "Arial") +
  theme(plot.background = element_rect(fill="transparent", colour=NA),
        panel.background = element_rect(fill="transparent", colour=NA),
        #aspect.ratio = 1/1,
        plot.margin = margin(.5,.5,.5,.5, "line"),
        axis.title.y = element_blank(),  axis.title.x = element_blank(), 
        axis.text = element_text(size = 5, colour = "black"),
        strip.text = element_text(size = 5, color = "black"), strip.background = element_blank(),
        legend.Yosition = "right", legend.justification = c("center", "center"), legend.background = element_rect(fill="transparent"), legend.box.background = element_rect(fill="transparent", color="transparent"),
        legend.text = element_text(size = 5), legend.margin = margin(0,0,0,0), legend.box.margin = margin(0,0,0,0), legend.spacing = unit(-.1, "char"), legend.title = element_text(size=5))

cowplot::plot_grid(p.chems, p.genes, align="hv", axis="tblr", nrow=2) %>%
  ggsave(filename = "red_UpsReg.pdf", device = cairo_pdf, width = 15, height = 15, units = "cm")

## YELLOW 5v1
integrated <- readRDS(file = "/Users/agodino/Desktop/inputs/Marine_integrated.rds")

dir.create("/Users/agodino/Desktop/forIPA_wBG")
setwd("/Users/agodino/Desktop/forIPA_wBG")
library(readxl)
library(tidyverse)

D1 <- subset(integrated, idents = "D1-MSNs") %>%
  GetAssayData(assay="SCT", layer="counts") %>%
  as.data.frame() %>%
  mutate_all(~ ifelse(. != 0, 1, 0))
ex.per <- 100*Matrix::rowSums(D1)/ncol(D1)
D1.BG <- data.frame(gene = names(ex.per), exp.percent.D1 = ex.per)

D2 <- subset(integrated, idents = "D2-MSNs") %>%
  GetAssayData(assay="SCT", layer="counts") %>%
  as.data.frame() %>%
  mutate_all(~ ifelse(. != 0, 1, 0))
ex.per <- 100*Matrix::rowSums(D2)/ncol(D2)
D2.BG <- data.frame(gene = names(ex.per), exp.percent.D2 = ex.per)

BG <- left_join(D1.BG, D2.BG)

D1.5v1.PY <- read_excel("/Users/agodino/Desktop/inputs/D1_POSYES_5Pv1P.xlsx") %>%
  select(gene, pct.1, pct.2, avg_log2FC, p_val_adj) %>%
  setNames(c("gene", "pct.1.D1.5v1.PY", "pct.2.D1.5v1.PY", "log2FC.D1.5v1.PY", "padj.D1.5v1.PY"))
D2.5v1.PY <- read_excel("/Users/agodino/Desktop/inputs/D2_POSYES_5Pv1P.xlsx") %>%
  select(gene, pct.1, pct.2, avg_log2FC, p_val_adj) %>%
  setNames(c("gene", "pct.1.D2.5v1.PY", "pct.2.D2.5v1.PY", "log2FC.D2.5v1.PY", "padj.D2.5v1.PY"))

all <- BG %>%
  left_join(D1.5v1.PY) %>%
  left_join(D2.5v1.PY) %>%
  mutate(across(contains("log2FC"), ~ ifelse(is.na(.), 0, .))) %>%
  mutate(across(contains("padj"), ~ ifelse(is.na(.), 1, .)))
#write.csv(all, "all_w_BG.csv")  

allD1 <- all %>% filter(exp.percent.D1>0) %>%
  select(matches("gene|D1")) %>%
  select(!matches("pct"))
allD2 <- all %>% filter(exp.percent.D2>0) %>%
  select(matches("gene|D2")) %>%
  select(!matches("pct"))
write.csv(allD1, file="yellow_D1.csv", row.names=F)
write.csv(allD2, file="yellow_D2.csv", row.names=F)

## CONTROL cluster marker genes

D1.markers <- read_excel("/Users/agodino/Desktop/inputs/D1_markers.xlsx") %>%
  select(gene, pct.1, pct.2, avg_log2FC, p_val_adj) %>%
  setNames(c("gene", "pct.1.D1", "pct.2.D1", "log2FC.D1", "padj.D1"))
D2.markers <- read_excel("/Users/agodino/Desktop/inputs/D2_markers.xlsx") %>%
  select(gene, pct.1, pct.2, avg_log2FC, p_val_adj) %>%
  setNames(c("gene", "pct.1.D2", "pct.2.D2", "log2FC.D2", "padj.D2"))

all <- BG %>%
  left_join(D1.markers) %>%
  left_join(D2.markers) %>%
  mutate(across(contains("log2FC"), ~ ifelse(is.na(.), 0, .))) %>%
  mutate(across(contains("padj"), ~ ifelse(is.na(.), 1, .)))
#write.csv(all, "all_w_BG.csv")  

allD1 <- all %>% filter(exp.percent.D1>0) %>%
  select(matches("gene|D1")) %>%
  select(!matches("pct"))
allD2 <- all %>% filter(exp.percent.D2>0) %>%
  select(matches("gene|D2")) %>%
  select(!matches("pct"))
write.csv(allD1, file="markers_D1.csv", row.names=F)
write.csv(allD2, file="markers_D2.csv", row.names=F)


## Comaprison with other IEGs R2
integrated <- readRDS(file = "/Users/agodino/Desktop/Marine_integrated.rds")

integrated$Egr1 <- ifelse(GetAssayData(object = integrated, assay = "SCT", slot = "data")["Egr1",] > 0, "YES", "NO")
integrated$Egr3 <- ifelse(GetAssayData(object = integrated, assay = "SCT", slot = "data")["Egr3",] > 0, "YES", "NO")
integrated$Fos <- ifelse(GetAssayData(object = integrated, assay = "SCT", slot = "data")["Fos",] > 0, "YES", "NO")
integrated$Fosb <- ifelse(GetAssayData(object = integrated, assay = "SCT", slot = "data")["Fosb",] > 0, "YES", "NO")
integrated$Fosl2 <- ifelse(GetAssayData(object = integrated, assay = "SCT", slot = "data")["Fosl2",] > 0, "YES", "NO")
integrated$Homer1 <- ifelse(GetAssayData(object = integrated, assay = "SCT", slot = "data")["Homer1",] > 0, "YES", "NO")
integrated$Junb <- ifelse(GetAssayData(object = integrated, assay = "SCT", slot = "data")["Junb",] > 0, "YES", "NO")
integrated$Npas4 <- ifelse(GetAssayData(object = integrated, assay = "SCT", slot = "data")["Npas4",] > 0, "YES", "NO")
integrated$Nr4a1 <- ifelse(GetAssayData(object = integrated, assay = "SCT", slot = "data")["Nr4a1",] > 0, "YES", "NO")
integrated$Nr4a3 <- ifelse(GetAssayData(object = integrated, assay = "SCT", slot = "data")["Nr4a3",] > 0, "YES", "NO")

all.IEGs <- data.frame(celltype=integrated$comb.clusters,
                       Arc=integrated$Arc.expr,
                       Egr1=integrated$Egr1,
                       Egr3=integrated$Egr3,
                       Fos=integrated$Fos,
                       Fosb=integrated$Fosb,
                       Fosl2=integrated$Fosl2,
                       Homer1=integrated$Homer1,
                       Junb=integrated$Junb,
                       Npas4=integrated$Npas4,
                       Nr4a1=integrated$Nr4a1,
                       Nr4a3=integrated$Nr4a3                       )
all.IEGs$yescount <- rowSums(all.IEGs=="YES")
arcplusany <- all.IEGs %>% filter(Arc=="YES" & yescount >= 2)
                    
sub <- all.IEGs %>% 
       filter(Arc=="YES"|Egr1=="YES"|Fos=="YES"|Npas4=="YES") %>%
       mutate(forVenny = paste0("Arc", Arc, "_Fos", Fos, "_Npas4", Npas4)) #"_Egr1", Egr1, 
table(sub$forVenny)



# Classifier
integrated <- readRDS(file = "/Users/agodino/Desktop/inputs/Marine_integrated.rds") #integrated@reductions$pca@cell.embeddings
df <-  integrated@assays$SCT$counts %>% t() %>% as.data.frame() %>% rownames_to_column("cellID") 
meta <- integrated$cluster_sorted_pairings_cpp %>% as.data.frame() %>% setNames("x") %>% rownames_to_column("cellID") %>%
        separate(x, into=c("cluster", "sorted", "pairings", "cpp"), sep = "_", remove=T) %>% as.data.frame()
df <- left_join(df, meta)

library(readxl)
D1.1P.PvN <- read_excel("/Users/agodino/Desktop/inputs/D1_1P_POSvNEG.xlsx") %>%
  select(gene, pct.1, pct.2, avg_log2FC, p_val_adj) %>%
  setNames(c("gene", "pct.1.D1.1P.PvN", "pct.2.D1.1P.PvN", "log2FC.D1.1P.PvN", "padj.D1.1P.PvN")) %>%
  filter(padj.D1.1P.PvN < 0)
D1.5P.PvN <- read_excel("/Users/agodino/Desktop/inputs/D1_5P_POSvNEG.xlsx") %>%
  select(gene, pct.1, pct.2, avg_log2FC, p_val_adj) %>%
  setNames(c("gene", "pct.1.D1.5P.PvN", "pct.2.D1.5P.PvN", "log2FC.D1.5P.PvN", "padj.D1.5P.PvN")) %>%
  filter(padj.D1.5P.PvN < 0)
D1.5v1.P <- read_excel("/Users/agodino/Desktop/inputs/D1_POS_5Pv1P.xlsx") %>%
  select(gene, pct.1, pct.2, avg_log2FC, p_val_adj) %>%
  setNames(c("gene", "pct.1.D1.5v1.P", "pct.2.D1.5v1.P", "log2FC.D1.5v1.P", "padj.D1.5v1.P"))%>%
  filter(padj.D1.5v1.P < .1)
D2.1P.PvN <- read_excel("/Users/agodino/Desktop/inputs/D2_1P_POSvNEG.xlsx") %>%
  select(gene, pct.1, pct.2, avg_log2FC, p_val_adj) %>%
  setNames(c("gene", "pct.1.D2.1P.PvN", "pct.2.D2.1P.PvN", "log2FC.D2.1P.PvN", "padj.D2.1P.PvN"))%>%
  filter(padj.D2.1P.PvN < 0)
D2.5P.PvN <- read_excel("/Users/agodino/Desktop/inputs/D2_5P_POSvNEG.xlsx") %>%
  select(gene, pct.1, pct.2, avg_log2FC, p_val_adj) %>%
  setNames(c("gene", "pct.1.D2.5P.PvN", "pct.2.D2.5P.PvN", "log2FC.D2.5P.PvN", "padj.D2.5P.PvN"))%>%
  filter(padj.D2.5P.PvN < 0)
D2.5v1.P <- read_excel("/Users/agodino/Desktop/inputs/D2_POS_5Pv1P.xlsx") %>%
  select(gene, pct.1, pct.2, avg_log2FC, p_val_adj) %>%
  setNames(c("gene", "pct.1.D2.5v1.P", "pct.2.D2.5v1.P", "log2FC.D2.5v1.P", "padj.D2.5v1.P"))%>%
  filter(padj.D2.5v1.P < .1)

D1VF <- full_join(D1.1P.PvN, D1.5P.PvN) %>%
        full_join(D1.5v1.P)
D1gf <- c(D1VF$gene, "pairings")

D2VF <- full_join(D2.1P.PvN, D2.5P.PvN) %>%
  full_join(D2.5v1.P)
D2gf <- c(D2VF$gene, "pairings")


Idents(integrated) <- "cluster_sorted_pairings_cpp"
D1.1P <- subset(integrated, idents = "D1-MSNs_POS_1P_HC") %>%
  GetAssayData(assay="SCT", layer="counts") %>%
  as.data.frame() %>%
  mutate_all(~ ifelse(. != 0, 1, 0))
ex.per <- 100*Matrix::rowSums(D1.1P)/ncol(D1.1P)
D1.1P.BG <- data.frame(gene = names(ex.per), exp.percent.D1.1P = ex.per)
D1.5P <- subset(integrated, idents = "D1-MSNs_POS_5P_HC") %>%
  GetAssayData(assay="SCT", layer="counts") %>%
  as.data.frame() %>%
  mutate_all(~ ifelse(. != 0, 1, 0))
ex.per <- 100*Matrix::rowSums(D1.5P)/ncol(D1.5P)
D1.5P.BG <- data.frame(gene = names(ex.per), exp.percent.D1.5P = ex.per)

D2.1P <- subset(integrated, idents = "D2-MSNs_POS_1P_HC") %>%
  GetAssayData(assay="SCT", layer="counts") %>%
  as.data.frame() %>%
  mutate_all(~ ifelse(. != 0, 1, 0))
ex.per <- 100*Matrix::rowSums(D2.1P)/ncol(D2.1P)
D2.1P.BG <- data.frame(gene = names(ex.per), exp.percent.D2.1P = ex.per)
D2.5P <- subset(integrated, idents = "D2-MSNs_POS_5P_HC") %>%
  GetAssayData(assay="SCT", layer="counts") %>%
  as.data.frame() %>%
  mutate_all(~ ifelse(. != 0, 1, 0))
ex.per <- 100*Matrix::rowSums(D2.5P)/ncol(D2.5P)
D2.5P.BG <- data.frame(gene = names(ex.per), exp.percent.D2.5P = ex.per)

BG <- full_join(D1.1P.BG, D1.5P.BG) %>% full_join(D2.1P.BG) %>% full_join(D2.5P.BG)

## Machine learning using SVM to predict outcome based on behavior
#Using all genes expressed in at least 30% of cells in one group
datasetD1 <- df %>% 
  filter(cluster=="D1-MSNs" & cpp=="HC" & sorted == "POS") %>%
  select(-cellID, -cluster, -cpp, -sorted) %>%
  mutate(pairings = factor(pairings)) %>%
  select(all_of(c(filter(BG, exp.percent.D1.1P > 30 | exp.percent.D1.5P > 30)$gene, "pairings")))
datasetD2 <- df %>% 
  filter(cluster=="D2-MSNs" & cpp=="HC" & sorted == "POS") %>%
  select(-cellID, -cluster, -cpp, -sorted) %>%
  mutate(pairings = factor(pairings)) %>%
  select(all_of(c(filter(BG, exp.percent.D2.1P > 30 | exp.percent.D2.5P > 30)$gene, "pairings")))

# Splitting the dataset into the Training set and Test set
options(expressions = 5e5)
customSVM <- function(ds, sr = .75, k = 5) {
  library(caTools)
  split = sample.split(ds$pairings, SplitRatio = sr)
  training_set = subset(ds, split == TRUE)
  test_set = subset(ds, split == FALSE)
  
  #training_set[-ncol(ds)] = scale(training_set[-ncol(ds)])
  #test_set[-ncol(ds)] = scale(test_set[-ncol(ds)])
  
  # Fitting SVM to the Training set
  library(e1071)
  tuned = tune.svm(pairings ~ ., cost = 10^(-3:1), type = 'C-classification', kernel = 'linear', data = training_set, tunecontrol= tune.control(cross=k))
  #classifier = svm(formula = AvorExp ~ ., data = training_set, type = 'C-classification', kernel = 'linear')
  svmfit = tuned$best.model
  #table(training_set[,c("AvorExp")], predict(svmfit))
  
  # Predicting the Test set results
  y_pred = predict(svmfit, newdata = select(test_set, -pairings))
  # Making the Confusion Matrix
  cm = table(test_set[,c("pairings")], y_pred) 
  return(as.data.frame(cm))
}

set.seed(124)
D1.cm.list <- list()
D2.cm.list <- list()
for (i in 1:10) {
  print(i)
  D1.cm.list[[i]] <- customSVM(datasetD1, sr = .75, k = 5) %>% mutate(iteration = i, CellType = "D1") %>% mutate(prop = 100*Freq/sum(.$Freq))
  print(i)
  D2.cm.list[[i]] <- customSVM(datasetD2, sr = .75, k = 5) %>% mutate(iteration = i, CellType = "D2") %>% mutate(prop = 100*Freq/sum(.$Freq))
}

accuracy <- rbind(rbindlist(D1.cm.list), rbindlist(D2.cm.list)) %>%
  mutate(pred.type = ifelse(Var1==y_pred, "concordant", "discordant")) %>%
  group_by(CellType, pred.type, iteration) %>% summarise(accuracy = sum(prop)) %>%
  filter(pred.type == "concordant")
msd <- accuracy %>% group_by(CellType) %>% summarize(sd = sd(accuracy, na.rm =T), sem = sd(accuracy, na.rm =T)/sqrt(length(accuracy)), accuracy=mean(accuracy, na.rm =T))

library(ggbeeswarm)
p.acu <- accuracy %>% ggplot(aes(x=CellType, y=accuracy)) + 
  geom_hline(yintercept = 50, linetype = 3) +
  geom_crossbar(data=msd, aes(ymin=accuracy, ymax=accuracy), color = "black", width = .4, fatten = 1) +
  geom_errorbar(data=msd, aes(ymin=accuracy-sem, ymax=accuracy+sem), color = "black", width = 0) + 
  geom_beeswarm(aes(color = CellType), size = .2) +
  scale_y_continuous(limits = c(0, 100), expand = c(0,0)) +
  scale_color_manual(values=c("#d55e00","#009e73")) +
  theme_classic(base_family = "Helvetica", base_size = 6) + 
  theme(plot.background = element_rect(fill="transparent", colour=NA),
        panel.background = element_rect(fill="transparent", colour=NA),
        plot.title = element_text(size = 6, hjust = .5),
        #aspect.ratio = 2/1,
        plot.margin = margin(.5,.5,.5,.5, "line"),
        axis.line = element_line(colour = "black", lineend = "square"),
        axis.title.y = element_text(margin = margin(0,.5,0,0, "line")),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(size = 5, colour = "black"),
        strip.background = element_rect(color = "black", fill = "transparent", size = 1),
        strip.text = element_text(size=6, colour = "black"),
        legend.position = c(.02, .98), legend.justification = c("left","top"), legend.key.size = unit(.5, "line"), 
        legend.text = element_text(size = 5), legend.box = "vertical", legend.title = element_blank())

#Using only 5P v 1P DEGs
datasetD1 <- df %>% 
  filter(cluster=="D1-MSNs" & cpp=="HC" & sorted == "POS") %>%
  select(-cellID, -cluster, -cpp, -sorted) %>%
  mutate(pairings = factor(pairings)) %>%
  select(all_of(D1gf))
table(datasetD1$pairings)
datasetD2 <- df %>% 
  filter(cluster=="D2-MSNs" & cpp=="HC" & sorted == "POS") %>%
  select(-cellID, -cluster, -cpp, -sorted) %>%
  mutate(pairings = factor(pairings)) %>%
  select(all_of(D2gf))
table(datasetD2$pairings)

# Splitting the dataset into the Training set and Test set
options(expressions = 5e5)
customSVM <- function(ds, sr = .75, k = 5) {
  library(caTools)
  split = sample.split(ds$pairings, SplitRatio = sr)
  training_set = subset(ds, split == TRUE)
  test_set = subset(ds, split == FALSE)
  
  #training_set[-ncol(training_set)] = scale(training_set[-ncol(training_set)])
  #test_set[-ncol(training_set)] = scale(test_set[-ncol(training_set)])
  
  # Fitting SVM to the Training set
  library(e1071)
  tuned = tune.svm(pairings ~ ., cost = 10^(-3:1), type = 'C-classification', kernel = 'linear', data = training_set, tunecontrol= tune.control(cross=k))
  #classifier = svm(formula = AvorExp ~ ., data = training_set, type = 'C-classification', kernel = 'linear')
  svmfit = tuned$best.model
  #table(training_set[,c("AvorExp")], predict(svmfit))
  
  # Predicting the Test set results
  y_pred = predict(svmfit, newdata = select(test_set, -pairings))
  # Making the Confusion Matrix
  cm = table(test_set[,c("pairings")], y_pred) 
  return(as.data.frame(cm))
}

set.seed(124)
D1.cm.list <- list()
D2.cm.list <- list()
for (i in 1:10) {
  print(i)
  D1.cm.list[[i]] <- customSVM(datasetD1, sr = .75, k = 5) %>% mutate(iteration = i, CellType = "D1") %>% mutate(prop = 100*Freq/sum(.$Freq))
  print(i)
  D2.cm.list[[i]] <- customSVM(datasetD2, sr = .75, k = 5) %>% mutate(iteration = i, CellType = "D2") %>% mutate(prop = 100*Freq/sum(.$Freq))
}

accuracy <- rbind(rbindlist(D1.cm.list), rbindlist(D2.cm.list)) %>%
  mutate(pred.type = ifelse(Var1==y_pred, "concordant", "discordant")) %>%
  group_by(CellType, pred.type, iteration) %>% summarise(accuracy = sum(prop)) %>%
  filter(pred.type == "concordant")
write.csv(accuracy, "/Users/agodino/Desktop/accuracy.5Pv1PDEGs.csv")

Fscore <- rbind(rbindlist(D1.cm.list), rbindlist(D2.cm.list)) %>%
  mutate(term = case_when(Var1=="1P" & y_pred=="1P" ~ "TN",
                          Var1=="5P" & y_pred=="1P" ~ "FN",
                          Var1=="1P" & y_pred=="5P" ~ "FP",
                          Var1=="5P" & y_pred=="5P" ~ "TP")) %>%
  pivot_wider(id_cols = c(CellType, iteration), names_from = term, values_from = Freq) %>%
  group_by(CellType, iteration) %>%
  summarize(precision = TP/(TP+FP), 
            recall = TP/(TP+FN)) %>%
  mutate(F1 = 2*(precision*recall)/(precision+recall)) %>%
  mutate(F1 = ifelse(is.na(F1), 0, F1))
write.csv(Fscore, "/Users/agodino/Desktop/Fscore.5Pv1PDEGs.csv")

msd <- Fscore %>% group_by(CellType) %>% summarize(sd = sd(F1, na.rm =T), sem = sd(F1, na.rm =T)/sqrt(length(F1)), F1=mean(F1, na.rm =T))

library(ggbeeswarm)
p.F1 <- Fscore %>% ggplot(aes(x=CellType, y=F1)) + 
  geom_hline(yintercept = .5, linetype = 3) +
  geom_crossbar(data=msd, aes(ymin=F1, ymax=F1), color = "black", width = .4, fatten = 1) +
  geom_errorbar(data=msd, aes(ymin=F1-sem, ymax=F1+sem), color = "black", width = 0) + 
  geom_beeswarm(aes(color = CellType), size = .2) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=c("#d55e00","#009e73")) +
  theme_classic(base_family = "Helvetica", base_size = 6) + 
  theme(plot.background = element_rect(fill="transparent", colour=NA),
        panel.background = element_rect(fill="transparent", colour=NA),
        plot.title = element_text(size = 6, hjust = .5),
        #aspect.ratio = 2/1,
        plot.margin = margin(.5,.5,.5,.5, "line"),
        axis.line = element_line(colour = "black", lineend = "square"),
        axis.title.y = element_text(margin = margin(0,.5,0,0, "line")),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(size = 5, colour = "black"),
        strip.background = element_rect(color = "black", fill = "transparent", size = 1),
        strip.text = element_text(size=6, colour = "black"),
        legend.position = c(.02, .98), legend.justification = c("left","top"), legend.key.size = unit(.5, "line"), 
        legend.text = element_text(size = 5), legend.box = "vertical", legend.title = element_blank())
p.F1

msd <- accuracy %>% group_by(CellType) %>% summarize(sd = sd(accuracy, na.rm =T), sem = sd(accuracy, na.rm =T)/sqrt(length(accuracy)), accuracy=mean(accuracy, na.rm =T))

library(ggbeeswarm)
p.acu <- accuracy %>% ggplot(aes(x=CellType, y=accuracy)) + 
  geom_hline(yintercept = 50, linetype = 3) +
  geom_crossbar(data=msd, aes(ymin=accuracy, ymax=accuracy), color = "black", width = .4, fatten = 1) +
  geom_errorbar(data=msd, aes(ymin=accuracy-sem, ymax=accuracy+sem), color = "black", width = 0) + 
  geom_beeswarm(aes(color = CellType), size = .2) +
  scale_y_continuous(limits = c(0, 100), expand = c(0,0)) +
  scale_color_manual(values=c("#d55e00","#009e73")) +
  theme_classic(base_family = "Helvetica", base_size = 6) + 
  theme(plot.background = element_rect(fill="transparent", colour=NA),
        panel.background = element_rect(fill="transparent", colour=NA),
        plot.title = element_text(size = 6, hjust = .5),
        #aspect.ratio = 2/1,
        plot.margin = margin(.5,.5,.5,.5, "line"),
        axis.line = element_line(colour = "black", lineend = "square"),
        axis.title.y = element_text(margin = margin(0,.5,0,0, "line")),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(size = 5, colour = "black"),
        strip.background = element_rect(color = "black", fill = "transparent", size = 1),
        strip.text = element_text(size=6, colour = "black"),
        legend.position = c(.02, .98), legend.justification = c("left","top"), legend.key.size = unit(.5, "line"), 
        legend.text = element_text(size = 5), legend.box = "vertical", legend.title = element_blank())
p.acu

# As controls, random selection of the same number of genes that are not DEGs
# D1 n DEGs = 87, D2 = 19

D1gf
D2gf
D1randuniverse <- BG %>% filter(!gene %in% D1gf)
D2randuniverse <- BG %>% filter(!gene %in% D2gf)
set.seed(126)
D1shuffle <- sample(D1randuniverse$gene, length(D1gf))
D2shuffle <- sample(D2randuniverse$gene, length(D2gf))

datasetD1 <- df %>% 
  filter(cluster=="D1-MSNs" & cpp=="HC" & sorted == "POS") %>%
  select(-cellID, -cluster, -cpp, -sorted) %>%
  mutate(pairings = factor(pairings)) %>%
  select(all_of(c(D1shuffle, "pairings")))
datasetD2 <- df %>% 
  filter(cluster=="D2-MSNs" & cpp=="HC" & sorted == "POS") %>%
  select(-cellID, -cluster, -cpp, -sorted) %>%
  mutate(pairings = factor(pairings)) %>%
  select(all_of(c(D2shuffle, "pairings")))

customSVM <- function(ds, sr = .75, k = 5) {
  library(caTools)
  split = sample.split(ds$pairings, SplitRatio = sr)
  training_set = subset(ds, split == TRUE)
  test_set = subset(ds, split == FALSE)
  
  #training_set[-ncol(ds)] = scale(training_set[-ncol(ds)])
  #test_set[-ncol(ds)] = scale(test_set[-ncol(ds)])
  
  # Fitting SVM to the Training set
  library(e1071)
  tuned = tune.svm(pairings ~ ., cost = 10^(-3:1), type = 'C-classification', kernel = 'linear', data = training_set, tunecontrol= tune.control(cross=k))
  #classifier = svm(formula = AvorExp ~ ., data = training_set, type = 'C-classification', kernel = 'linear')
  svmfit = tuned$best.model
  #table(training_set[,c("AvorExp")], predict(svmfit))
  
  # Predicting the Test set results
  y_pred = predict(svmfit, newdata = select(test_set, -pairings))
  # Making the Confusion Matrix
  cm = table(test_set[,c("pairings")], y_pred) 
  return(as.data.frame(cm))
}

set.seed(124)
D1.cm.list <- list()
D2.cm.list <- list()
for (i in 1:10) {
  print(i)
  D1.cm.list[[i]] <- customSVM(datasetD1, sr = .75, k = 5) %>% mutate(iteration = i, CellType = "D1") %>% mutate(prop = 100*Freq/sum(.$Freq))
  print(i)
  D2.cm.list[[i]] <- customSVM(datasetD2, sr = .75, k = 5) %>% mutate(iteration = i, CellType = "D2") %>% mutate(prop = 100*Freq/sum(.$Freq))
}

accuracy <- rbind(rbindlist(D1.cm.list), rbindlist(D2.cm.list)) %>%
  mutate(pred.type = ifelse(Var1==y_pred, "concordant", "discordant")) %>%
  group_by(CellType, pred.type, iteration) %>% summarise(accuracy = sum(prop)) %>%
  filter(pred.type == "concordant")
write.csv(accuracy, "/Users/agodino/Desktop/accuracy.5Pv1PDEGs_shuffle.csv")

Fscore <- rbind(rbindlist(D1.cm.list), rbindlist(D2.cm.list)) %>%
  mutate(term = case_when(Var1=="1P" & y_pred=="1P" ~ "TN",
                          Var1=="5P" & y_pred=="1P" ~ "FN",
                          Var1=="1P" & y_pred=="5P" ~ "FP",
                          Var1=="5P" & y_pred=="5P" ~ "TP")) %>%
  pivot_wider(id_cols = c(CellType, iteration), names_from = term, values_from = Freq) %>%
  group_by(CellType, iteration) %>%
  summarize(precision = TP/(TP+FP), 
            recall = TP/(TP+FN)) %>%
  mutate(F1 = 2*(precision*recall)/(precision+recall)) %>%
  mutate(F1 = ifelse(is.na(F1), 0, F1))
write.csv(Fscore, "/Users/agodino/Desktop/Fscore.5Pv1P_shuffle.csv")

msd <- Fscore %>% group_by(CellType) %>% summarize(sd = sd(F1, na.rm =T), sem = sd(F1, na.rm =T)/sqrt(length(F1)), F1=mean(F1, na.rm =T))

library(ggbeeswarm)
p.F1 <- Fscore %>% ggplot(aes(x=CellType, y=F1)) + 
  geom_hline(yintercept = .5, linetype = 3) +
  geom_crossbar(data=msd, aes(ymin=F1, ymax=F1), color = "black", width = .4, fatten = 1) +
  geom_errorbar(data=msd, aes(ymin=F1-sem, ymax=F1+sem), color = "black", width = 0) + 
  geom_beeswarm(aes(color = CellType), size = .2) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=c("#d55e00","#009e73")) +
  theme_classic(base_family = "Helvetica", base_size = 6) + 
  theme(plot.background = element_rect(fill="transparent", colour=NA),
        panel.background = element_rect(fill="transparent", colour=NA),
        plot.title = element_text(size = 6, hjust = .5),
        #aspect.ratio = 2/1,
        plot.margin = margin(.5,.5,.5,.5, "line"),
        axis.line = element_line(colour = "black", lineend = "square"),
        axis.title.y = element_text(margin = margin(0,.5,0,0, "line")),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(size = 5, colour = "black"),
        strip.background = element_rect(color = "black", fill = "transparent", size = 1),
        strip.text = element_text(size=6, colour = "black"),
        legend.position = c(.02, .98), legend.justification = c("left","top"), legend.key.size = unit(.5, "line"), 
        legend.text = element_text(size = 5), legend.box = "vertical", legend.title = element_blank())
p.F1

msd <- accuracy %>% group_by(CellType) %>% summarize(sd = sd(accuracy, na.rm =T), sem = sd(accuracy, na.rm =T)/sqrt(length(accuracy)), accuracy=mean(accuracy, na.rm =T))

library(ggbeeswarm)
p.acu <- accuracy %>% ggplot(aes(x=CellType, y=accuracy)) + 
  geom_hline(yintercept = 50, linetype = 3) +
  geom_crossbar(data=msd, aes(ymin=accuracy, ymax=accuracy), color = "black", width = .4, fatten = 1) +
  geom_errorbar(data=msd, aes(ymin=accuracy-sem, ymax=accuracy+sem), color = "black", width = 0) + 
  geom_beeswarm(aes(color = CellType), size = .2) +
  scale_y_continuous(limits = c(0, 100), expand = c(0,0)) +
  scale_color_manual(values=c("#d55e00","#009e73")) +
  theme_classic(base_family = "Helvetica", base_size = 6) + 
  theme(plot.background = element_rect(fill="transparent", colour=NA),
        panel.background = element_rect(fill="transparent", colour=NA),
        plot.title = element_text(size = 6, hjust = .5),
        #aspect.ratio = 2/1,
        plot.margin = margin(.5,.5,.5,.5, "line"),
        axis.line = element_line(colour = "black", lineend = "square"),
        axis.title.y = element_text(margin = margin(0,.5,0,0, "line")),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(size = 5, colour = "black"),
        strip.background = element_rect(color = "black", fill = "transparent", size = 1),
        strip.text = element_text(size=6, colour = "black"),
        legend.position = c(.02, .98), legend.justification = c("left","top"), legend.key.size = unit(.5, "line"), 
        legend.text = element_text(size = 5), legend.box = "vertical", legend.title = element_blank())
p.acu

## YELLOW, predict reactivation ie YELLOW over all GREEN

integrated <- readRDS(file = "/Users/agodino/Desktop/inputs/Marine_integrated.rds")
df <-  integrated@assays$SCT$counts %>% t() %>% as.data.frame() %>% rownames_to_column("cellID") 
meta <- integrated$cluster_pairings_cpp_sorted_Arc %>% as.data.frame() %>% setNames("x") %>% rownames_to_column("cellID") %>%
  separate(x, into=c("cluster", "pairings", "cpp", "sorted", "ArcYN"), sep = "_", remove=T) %>% as.data.frame()
df <- left_join(df, meta)

D1.1P.react <- c("Mdh1", "Adcy9","Myo5a","Gm28791","Skil","Gpm6a","Cxadr","Nol4","Cacng3","Cul1","Prkca","Sp3","Tnpo1","Nin","Smarca4","Klf12","Kmt2c","Ptk2","Cdh18","Zfp609","Lingo2","Ptchd4")
D1.5P.react <- c("Map4k3","Nvl","Tbc1d5","Josd1","Cnnm2","Map3k3","Brd8","Lars2","Thsd7a","5330438D12Rik","Ogdhl","Dnm3","Apbb2","Rlf","Ctnna2","Rab2a","C78859","Zfp638","Clpx","Mtmr2","Klhl3","Vat1l","Cnot6l","Ptchd1","Mir6236","Hace1","Rasa1")
D2.1P.react <- c("Mcu","Ehbp1","Dennd1a","Lzts1","Med13l","Cenpv","Zc3h12b","Sptbn4","Slc35f3","Zcchc7","Lpp","5330438D12Rik","Gpm6a")
D2.5P.react <- c("Miip")

datasetD1.1P <- df %>% 
  filter(cluster=="D1-MSNs" & cpp=="Test" & pairings == "1P" & sorted == "POS") %>%
  select(-cellID, -cluster, -cpp, -pairings, -sorted) %>%
  mutate(ArcYN = factor(ArcYN)) %>%
  select(all_of(c(D1.1P.react, "ArcYN")))
table(datasetD1.1P$ArcYN)
datasetD1.5P <- df %>% 
  filter(cluster=="D1-MSNs" & cpp=="Test" & pairings == "5P" & sorted == "POS") %>%
  select(-cellID, -cluster, -cpp, -pairings, -sorted) %>%
  mutate(ArcYN = factor(ArcYN)) %>%
  select(all_of(c(D1.5P.react, "ArcYN")))
table(datasetD1.5P$ArcYN)
datasetD2.1P <- df %>% 
  filter(cluster=="D2-MSNs" & cpp=="Test" & pairings == "1P" & sorted == "POS") %>%
  select(-cellID, -cluster, -cpp, -pairings, -sorted) %>%
  mutate(ArcYN = factor(ArcYN)) %>%
  select(all_of(c(D2.1P.react, "ArcYN")))
table(datasetD2.1P$ArcYN)
datasetD2.5P <- df %>% 
  filter(cluster=="D2-MSNs" & cpp=="Test" & pairings == "5P" & sorted == "POS") %>%
  select(-cellID, -cluster, -cpp, -pairings, -sorted) %>%
  mutate(ArcYN = factor(ArcYN)) %>%
  select(all_of(c(D2.5P.react, "ArcYN")))
table(datasetD2.5P$ArcYN)

# Splitting the dataset into the Training set and Test set
options(expressions = 5e5)
customSVM <- function(ds, sr = .75, k = 5) {
  library(caTools)
  split = sample.split(ds$ArcYN, SplitRatio = sr)
  training_set = subset(ds, split == TRUE)
  test_set = subset(ds, split == FALSE)
  
  #training_set[-ncol(training_set)] = scale(training_set[-ncol(training_set)])
  #test_set[-ncol(training_set)] = scale(test_set[-ncol(training_set)])
  
  # Fitting SVM to the Training set
  library(e1071)
  tuned = tune.svm(ArcYN ~ ., cost = 10^(-3:1), type = 'C-classification', kernel = 'linear', data = training_set, tunecontrol= tune.control(cross=k))
  #classifier = svm(formula = AvorExp ~ ., data = training_set, type = 'C-classification', kernel = 'linear')
  svmfit = tuned$best.model
  #table(training_set[,c("AvorExp")], predict(svmfit))
  
  # Predicting the Test set results
  y_pred = predict(svmfit, newdata = select(test_set, -ArcYN))
  # Making the Confusion Matrix
  cm = table(test_set[,c("ArcYN")], y_pred) 
  return(as.data.frame(cm))
}

set.seed(124)
D1.1P.cm.list <- list()
D2.1P.cm.list <- list()
D1.5P.cm.list <- list()
D2.5P.cm.list <- list()
for (i in 1:10) {
  print(i)
  D1.1P.cm.list[[i]] <- customSVM(datasetD1.1P, sr = .75, k = 5) %>% mutate(iteration = i, CellType = "D1", pairings = "1P") %>% mutate(prop = 100*Freq/sum(.$Freq))
  print(i)
  D2.1P.cm.list[[i]] <- customSVM(datasetD2.1P, sr = .75, k = 5) %>% mutate(iteration = i, CellType = "D2", pairings = "1P") %>% mutate(prop = 100*Freq/sum(.$Freq))
  print(i)
  D1.5P.cm.list[[i]] <- customSVM(datasetD1.5P, sr = .75, k = 5) %>% mutate(iteration = i, CellType = "D1", pairings = "5P") %>% mutate(prop = 100*Freq/sum(.$Freq))
  print(i)
  D2.5P.cm.list[[i]] <- customSVM(datasetD2.5P, sr = .75, k = 5) %>% mutate(iteration = i, CellType = "D2", pairings = "5P") %>% mutate(prop = 100*Freq/sum(.$Freq))
}

accuracy <- rbind(rbindlist(D1.1P.cm.list), rbindlist(D2.1P.cm.list), rbindlist(D1.5P.cm.list), rbindlist(D2.5P.cm.list)) %>%
  mutate(pred.type = ifelse(Var1==y_pred, "concordant", "discordant")) %>%
  group_by(CellType, pairings, pred.type, iteration) %>% summarise(accuracy = sum(prop)) %>%
  filter(pred.type == "concordant") %>%
  mutate(celltype_pairing = paste(CellType, pairings, sep="_"))
write.csv(accuracy, "/Users/agodino/Desktop/accuracy.react.csv")


Fscore <- rbind(rbindlist(D1.1P.cm.list), rbindlist(D2.1P.cm.list), rbindlist(D1.5P.cm.list), rbindlist(D2.5P.cm.list)) %>%
  mutate(term = case_when(Var1=="NO" & y_pred=="NO" ~ "TN",
                          Var1=="YES" & y_pred=="NO" ~ "FN",
                          Var1=="NO" & y_pred=="YES" ~ "FP",
                          Var1=="YES" & y_pred=="YES" ~ "TP")) %>%
  pivot_wider(id_cols = c(CellType, pairings, iteration), names_from = term, values_from = Freq) %>%
  group_by(CellType, pairings, iteration) %>%
  summarize(precision = TP/(TP+FP), 
            recall = TP/(TP+FN)) %>%
  mutate(F1 = 2*(precision*recall)/(precision+recall)) %>%
  mutate(F1 = ifelse(is.na(F1), 0, F1)) %>%
  mutate(celltype_pairing = paste(CellType, pairings, sep="_"))
write.csv(accuracy, "/Users/agodino/Desktop/Fscore.react.csv")

msd <- Fscore %>% group_by(CellType, celltype_pairing) %>% summarize(sd = sd(F1, na.rm =T), sem = sd(F1, na.rm =T)/sqrt(length(F1)), F1=mean(F1, na.rm =T))

library(ggbeeswarm)
p.F1 <- Fscore %>% ggplot(aes(x=celltype_pairing, y=F1)) + 
  geom_hline(yintercept = .5, linetype = 3) +
  geom_crossbar(data=msd, aes(ymin=F1, ymax=F1), color = "black", width = .4, fatten = 1) +
  geom_errorbar(data=msd, aes(ymin=F1-sem, ymax=F1+sem), color = "black", width = 0) + 
  geom_beeswarm(aes(color = CellType), size = .2) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=c("#d55e00","#009e73")) +
  theme_classic(base_family = "Helvetica", base_size = 6) + 
  theme(plot.background = element_rect(fill="transparent", colour=NA),
        panel.background = element_rect(fill="transparent", colour=NA),
        plot.title = element_text(size = 6, hjust = .5),
        #aspect.ratio = 2/1,
        plot.margin = margin(.5,.5,.5,.5, "line"),
        axis.line = element_line(colour = "black", lineend = "square"),
        axis.title.y = element_text(margin = margin(0,.5,0,0, "line")),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(size = 5, colour = "black"),
        strip.background = element_rect(color = "black", fill = "transparent", size = 1),
        strip.text = element_text(size=6, colour = "black"),
        legend.position = c(.02, .98), legend.justification = c("left","top"), legend.key.size = unit(.5, "line"), 
        legend.text = element_text(size = 5), legend.box = "vertical", legend.title = element_blank())
p.F1

msd <- accuracy %>% group_by(CellType, celltype_pairing) %>% summarize(sd = sd(accuracy, na.rm =T), sem = sd(accuracy, na.rm =T)/sqrt(length(accuracy)), accuracy=mean(accuracy, na.rm =T))
library(ggbeeswarm)
p.acu <- accuracy %>% ggplot(aes(x=celltype_pairing, y=accuracy)) + 
  geom_hline(yintercept = 50, linetype = 3) +
  geom_crossbar(data=msd, aes(ymin=accuracy, ymax=accuracy), color = "black", width = .4, fatten = 1) +
  geom_errorbar(data=msd, aes(ymin=accuracy-sem, ymax=accuracy+sem), color = "black", width = 0) + 
  geom_beeswarm(aes(color = CellType), size = .2) +
  scale_y_continuous(limits = c(0, 100), expand = c(0,0)) +
  scale_color_manual(values=c("#d55e00","#009e73")) +
  theme_classic(base_family = "Helvetica", base_size = 6) + 
  theme(plot.background = element_rect(fill="transparent", colour=NA),
        panel.background = element_rect(fill="transparent", colour=NA),
        plot.title = element_text(size = 6, hjust = .5),
        #aspect.ratio = 2/1,
        plot.margin = margin(.5,.5,.5,.5, "line"),
        axis.line = element_line(colour = "black", lineend = "square"),
        axis.title.y = element_text(margin = margin(0,.5,0,0, "line")),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(size = 5, colour = "black"),
        strip.background = element_rect(color = "black", fill = "transparent", size = 1),
        strip.text = element_text(size=6, colour = "black"),
        legend.position = c(.02, .98), legend.justification = c("left","top"), legend.key.size = unit(.5, "line"), 
        legend.text = element_text(size = 5), legend.box = "vertical", legend.title = element_blank())
p.acu



# random shuffle

Idents(integrated) <- "cluster_sorted_pairings_cpp"
D1.1P <- subset(integrated, idents = "D1-MSNs_POS_1P_Test") %>%
  GetAssayData(assay="SCT", layer="counts") %>%
  as.data.frame() %>%
  mutate_all(~ ifelse(. != 0, 1, 0))
ex.per <- 100*Matrix::rowSums(D1.1P)/ncol(D1.1P)
D1.1P.BG <- data.frame(gene = names(ex.per), exp.percent.D1.1P = ex.per)
D1.5P <- subset(integrated, idents = "D1-MSNs_POS_5P_Test") %>%
  GetAssayData(assay="SCT", layer="counts") %>%
  as.data.frame() %>%
  mutate_all(~ ifelse(. != 0, 1, 0))
ex.per <- 100*Matrix::rowSums(D1.5P)/ncol(D1.5P)
D1.5P.BG <- data.frame(gene = names(ex.per), exp.percent.D1.5P = ex.per)

D2.1P <- subset(integrated, idents = "D2-MSNs_POS_1P_Test") %>%
  GetAssayData(assay="SCT", layer="counts") %>%
  as.data.frame() %>%
  mutate_all(~ ifelse(. != 0, 1, 0))
ex.per <- 100*Matrix::rowSums(D2.1P)/ncol(D2.1P)
D2.1P.BG <- data.frame(gene = names(ex.per), exp.percent.D2.1P = ex.per)
D2.5P <- subset(integrated, idents = "D2-MSNs_POS_5P_Test") %>%
  GetAssayData(assay="SCT", layer="counts") %>%
  as.data.frame() %>%
  mutate_all(~ ifelse(. != 0, 1, 0))
ex.per <- 100*Matrix::rowSums(D2.5P)/ncol(D2.5P)
D2.5P.BG <- data.frame(gene = names(ex.per), exp.percent.D2.5P = ex.per)

BG <- full_join(D1.1P.BG, D1.5P.BG) %>% full_join(D2.1P.BG) %>% full_join(D2.5P.BG)

D1.1Pranduniverse <- BG %>% filter(!gene %in% D1.1P.react)
D2.1Pranduniverse <- BG %>% filter(!gene %in% D2.1P.react)
D1.5Pranduniverse <- BG %>% filter(!gene %in% D1.5P.react)
D2.5Pranduniverse <- BG %>% filter(!gene %in% D2.5P.react)
set.seed(126)
D1.1Pshuffle <- sample(D1.1Pranduniverse$gene, length(D1.1P.react))
D2.1Pshuffle <- sample(D2.1Pranduniverse$gene, length(D2.1P.react))
D1.5Pshuffle <- sample(D1.5Pranduniverse$gene, length(D1.5P.react))
D2.5Pshuffle <- sample(D2.5Pranduniverse$gene, length(D2.5P.react))


datasetD1.1P <- df %>% 
  filter(cluster=="D1-MSNs" & cpp=="Test" & pairings == "1P" & sorted == "POS") %>%
  select(-cellID, -cluster, -cpp, -pairings, -sorted) %>%
  mutate(ArcYN = factor(ArcYN)) %>%
  select(all_of(c(D1.1Pshuffle, "ArcYN")))
datasetD1.5P <- df %>% 
  filter(cluster=="D1-MSNs" & cpp=="Test" & pairings == "5P" & sorted == "POS") %>%
  select(-cellID, -cluster, -cpp, -pairings, -sorted) %>%
  mutate(ArcYN = factor(ArcYN)) %>%
  select(all_of(c(D1.5Pshuffle, "ArcYN")))
datasetD2.1P <- df %>% 
  filter(cluster=="D2-MSNs" & cpp=="Test" & pairings == "1P" & sorted == "POS") %>%
  select(-cellID, -cluster, -cpp, -pairings, -sorted) %>%
  mutate(ArcYN = factor(ArcYN)) %>%
  select(all_of(c(D2.1Pshuffle, "ArcYN")))
datasetD2.5P <- df %>% 
  filter(cluster=="D2-MSNs" & cpp=="Test" & pairings == "5P" & sorted == "POS") %>%
  select(-cellID, -cluster, -cpp, -pairings, -sorted) %>%
  mutate(ArcYN = factor(ArcYN)) %>%
  select(all_of(c(D2.5Pshuffle, "ArcYN")))


# Splitting the dataset into the Training set and Test set
options(expressions = 5e5)
customSVM <- function(ds, sr = .75, k = 5) {
  library(caTools)
  split = sample.split(ds$ArcYN, SplitRatio = sr)
  training_set = subset(ds, split == TRUE)
  test_set = subset(ds, split == FALSE)
  
  #training_set[-ncol(training_set)] = scale(training_set[-ncol(training_set)])
  #test_set[-ncol(training_set)] = scale(test_set[-ncol(training_set)])
  
  # Fitting SVM to the Training set
  library(e1071)
  tuned = tune.svm(ArcYN ~ ., cost = 10^(-3:1), type = 'C-classification', kernel = 'linear', data = training_set, tunecontrol= tune.control(cross=k))
  #classifier = svm(formula = AvorExp ~ ., data = training_set, type = 'C-classification', kernel = 'linear')
  svmfit = tuned$best.model
  #table(training_set[,c("AvorExp")], predict(svmfit))
  
  # Predicting the Test set results
  y_pred = predict(svmfit, newdata = select(test_set, -ArcYN))
  # Making the Confusion Matrix
  cm = table(test_set[,c("ArcYN")], y_pred) 
  return(as.data.frame(cm))
}

set.seed(124)
D1.1P.cm.list <- list()
D2.1P.cm.list <- list()
D1.5P.cm.list <- list()
D2.5P.cm.list <- list()
for (i in 1:10) {
  print(i)
  D1.1P.cm.list[[i]] <- customSVM(datasetD1.1P, sr = .75, k = 5) %>% mutate(iteration = i, CellType = "D1", pairings = "1P") %>% mutate(prop = 100*Freq/sum(.$Freq))
  print(i)
  D2.1P.cm.list[[i]] <- customSVM(datasetD2.1P, sr = .75, k = 5) %>% mutate(iteration = i, CellType = "D2", pairings = "1P") %>% mutate(prop = 100*Freq/sum(.$Freq))
  print(i)
  D1.5P.cm.list[[i]] <- customSVM(datasetD1.5P, sr = .75, k = 5) %>% mutate(iteration = i, CellType = "D1", pairings = "5P") %>% mutate(prop = 100*Freq/sum(.$Freq))
  print(i)
  D2.5P.cm.list[[i]] <- customSVM(datasetD2.5P, sr = .75, k = 5) %>% mutate(iteration = i, CellType = "D2", pairings = "5P") %>% mutate(prop = 100*Freq/sum(.$Freq))
}

accuracy <- rbind(rbindlist(D1.1P.cm.list), rbindlist(D2.1P.cm.list), rbindlist(D1.5P.cm.list), rbindlist(D2.5P.cm.list)) %>%
  mutate(pred.type = ifelse(Var1==y_pred, "concordant", "discordant")) %>%
  group_by(CellType, pairings, pred.type, iteration) %>% summarise(accuracy = sum(prop)) %>%
  filter(pred.type == "concordant") %>%
  mutate(celltype_pairing = paste(CellType, pairings, sep="_"))
write.csv(accuracy, "/Users/agodino/Desktop/accuracy.react.shuffle.csv")

Fscore <- rbind(rbindlist(D1.1P.cm.list), rbindlist(D2.1P.cm.list), rbindlist(D1.5P.cm.list), rbindlist(D2.5P.cm.list)) %>%
  mutate(term = case_when(Var1=="NO" & y_pred=="NO" ~ "TN",
                 Var1=="YES" & y_pred=="NO" ~ "FN",
                 Var1=="NO" & y_pred=="YES" ~ "FP",
                 Var1=="YES" & y_pred=="YES" ~ "TP")) %>%
  pivot_wider(id_cols = c(CellType, pairings, iteration), names_from = term, values_from = Freq) %>%
  group_by(CellType, pairings, iteration) %>%
  summarize(precision = TP/(TP+FP), 
            recall = TP/(TP+FN)) %>%
  mutate(F1 = 2*(precision*recall)/(precision+recall)) %>%
  mutate(F1 = ifelse(is.na(F1), 0, F1)) %>%
  mutate(celltype_pairing = paste(CellType, pairings, sep="_"))
write.csv(accuracy, "/Users/agodino/Desktop/Fscore.react.shuffle.csv")


Fscore <- rbind(rbindlist(D1.1P.cm.list), rbindlist(D2.1P.cm.list), rbindlist(D1.5P.cm.list), rbindlist(D2.5P.cm.list)) %>%
  mutate(term = case_when(Var1=="NO" & y_pred=="NO" ~ "TN",
                          Var1=="YES" & y_pred=="NO" ~ "FN",
                          Var1=="NO" & y_pred=="YES" ~ "FP",
                          Var1=="YES" & y_pred=="YES" ~ "TP")) %>%
  pivot_wider(id_cols = c(CellType, pairings, iteration), names_from = term, values_from = Freq) %>%
  group_by(CellType, pairings, iteration) %>%
  summarize(precision = TP/(TP+FP), 
            recall = TP/(TP+FN)) %>%
  mutate(F1 = 2*(precision*recall)/(precision+recall)) %>%
  mutate(F1 = ifelse(is.na(F1), 0, F1)) %>%
  mutate(celltype_pairing = paste(CellType, pairings, sep="_"))
write.csv(accuracy, "/Users/agodino/Desktop/Fscore.react.shuffle.csv")

msd <- Fscore %>% group_by(CellType, celltype_pairing) %>% summarize(sd = sd(F1, na.rm =T), sem = sd(F1, na.rm =T)/sqrt(length(F1)), F1=mean(F1, na.rm =T))

library(ggbeeswarm)
p.F1 <- Fscore %>% ggplot(aes(x=celltype_pairing, y=F1)) + 
  geom_hline(yintercept = .5, linetype = 3) +
  geom_crossbar(data=msd, aes(ymin=F1, ymax=F1), color = "black", width = .4, fatten = 1) +
  geom_errorbar(data=msd, aes(ymin=F1-sem, ymax=F1+sem), color = "black", width = 0) + 
  geom_beeswarm(aes(color = CellType), size = .2) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(values=c("#d55e00","#009e73")) +
  theme_classic(base_family = "Helvetica", base_size = 6) + 
  theme(plot.background = element_rect(fill="transparent", colour=NA),
        panel.background = element_rect(fill="transparent", colour=NA),
        plot.title = element_text(size = 6, hjust = .5),
        #aspect.ratio = 2/1,
        plot.margin = margin(.5,.5,.5,.5, "line"),
        axis.line = element_line(colour = "black", lineend = "square"),
        axis.title.y = element_text(margin = margin(0,.5,0,0, "line")),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(size = 5, colour = "black"),
        strip.background = element_rect(color = "black", fill = "transparent", size = 1),
        strip.text = element_text(size=6, colour = "black"),
        legend.position = c(.02, .98), legend.justification = c("left","top"), legend.key.size = unit(.5, "line"), 
        legend.text = element_text(size = 5), legend.box = "vertical", legend.title = element_blank())
p.F1


# Dopamine receptors
integrated <- readRDS(file = "/Users/agodino/Desktop/inputs/Marine_integrated.rds")
Idents(integrated) <- "cpp"
integrated.HC <- subset(integrated, idents = "HC")
Idents(integrated.HC) <- "comb.clusters"

tempD1 <- subset(integrated.HC, idents = "D1-MSNs")
ex <- DotPlot(tempD1, features = rev(c("Drd1", "Drd2", "Drd3", "Drd4", "Drd5")), 
              #cluster.idents = F,
              assay = "SCT", dot.scale = 6, group.by = "cluster_sorted_pairings_cpp",
              cols = c("lightgrey", "#F564E3"),
              dot.min = .01, col.min = 0, col.max = 1, #scale.min = .1
              ) + 
  coord_flip() +
  #scale_y_discrete(labels = levels(integrated$comb.clusters)) +
  theme_classic(base_family = "Arial", base_size = 5) +
  theme(plot.background = element_rect(fill="transparent", colour=NA),
        panel.background = element_rect(fill="transparent", colour=NA),
        aspect.ratio = 8/4, axis.title = element_blank(), axis.text = element_text(color = "black"), axis.ticks = element_line(color="black"),
        axis.text.x = element_text(face = "italic", angle = 45, vjust = 0.8, hjust = 1),
        legend.key.size = unit(.5, "line"), legend.position = "bottom", legend.box = "vertical")
ggsave(ex, filename = "/Users/agodino/Desktop/Drd_green_D1_DotPlot.pdf", device = cairo_pdf, width = 9, height = 9, units = "cm") 

tempD2 <- subset(integrated.HC, idents = "D2-MSNs")
ex <- DotPlot(tempD2, features = rev(c("Drd1", "Drd2", "Drd3", "Drd4", "Drd5")), 
              #cluster.idents = F,
              assay = "SCT", dot.scale = 6, group.by = "cluster_sorted_pairings_cpp",
              cols = c("lightgrey", "#619CFF"),
              dot.min = .01, col.min = 0, col.max = 1, #scale.min = .1
) + 
  coord_flip() +
  #scale_y_discrete(labels = levels(integrated$comb.clusters)) +
  theme_classic(base_family = "Arial", base_size = 5) +
  theme(plot.background = element_rect(fill="transparent", colour=NA),
        panel.background = element_rect(fill="transparent", colour=NA),
        aspect.ratio = 8/4, axis.title = element_blank(), axis.text = element_text(color = "black"), axis.ticks = element_line(color="black"),
        axis.text.x = element_text(face = "italic", angle = 45, vjust = 0.8, hjust = 1),
        legend.key.size = unit(.5, "line"), legend.position = "bottom", legend.box = "vertical")
ggsave(ex, filename = "/Users/agodino/Desktop/Drd_green_D2_DotPlot.pdf", device = cairo_pdf, width = 9, height = 9, units = "cm") 




# Dot plots example genes
integrated.Test <- subset(integrated, idents = "Test")
Idents(integrated.Test) <- "comb.clusters"
tempD1 <- subset(integrated.Test, idents = "D1-MSNs")
ex.D1 <- DotPlot(tempD1, features = rev(c("Drd1", "Drd2", "Drd3", "Drd4", "Drd5" )), 
                 #cluster.idents = F,
                 assay = "SCT", dot.scale = 6, group.by = "cluster_pairings_cpp_Arc",
                 cols = c("lightgrey", "#F564E3"),
                 dot.min = .01, col.min = 0, col.max = 1, #scale.min = .1, scale.max = 100
                 ) + 
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
ex.D2 <- DotPlot(tempD2, features = rev(c("Drd1", "Drd2", "Drd3", "Drd4", "Drd5" )), 
                 #cluster.idents = F,
                 assay = "SCT", dot.scale = 6, group.by = "cluster_pairings_cpp_Arc",
                 cols = c("lightgrey", "#619CFF"),
                 dot.min = .01, col.min = 0, col.max = 1#, scale.min = .1, scale.max = 100
                 ) + 
  coord_flip() +
  scale_y_discrete(limits=c("D2-MSNs_1P_Test_NO", "D2-MSNs_5P_Test_NO", "D2-MSNs_1P_Test_YES", "D2-MSNs_5P_Test_YES")) +
  theme_classic(base_family = "Arial", base_size = 5) +
  theme(plot.background = element_rect(fill="transparent", colour=NA),
        panel.background = element_rect(fill="transparent", colour=NA),
        #aspect.ratio = 13/4, 
        axis.title = element_blank(), axis.text = element_text(color = "black"), axis.ticks = element_line(color="black"),
        axis.text.x = element_text(face = "italic", angle = 45, vjust = 0.8, hjust = 1),
        legend.key.size = unit(.5, "line"), legend.position = "bottom", legend.box = "vertical")

ex<-cowplot::plot_grid(ex.D1, ex.D2, align = "hv", axis = "tblr", ncol=2)
ggsave(ex, filename = "/Users/agodino/Desktop/ex_red_DotPlot.pdf", device = cairo_pdf, width = 9, height = 9, units = "cm") 


# Dot plots example genes
Idents(integrated.Test) <- "comb.clusters"
tempD1 <- subset(integrated.Test, idents = "D1-MSNs")
ex.D1 <- DotPlot(tempD1, features = rev(c("Drd1", "Drd2", "Drd3", "Drd4", "Drd5")), 
                 #cluster.idents = F,
                 assay = "SCT", dot.scale = 6, group.by = "cluster_pairings_cpp_sorted_Arc",
                 cols = c("lightgrey", "#F564E3"),
                 dot.min = .01, col.min = 0, col.max = 1#, scale.min = .1, scale.max = 100
                 ) + 
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
ex.D2 <- DotPlot(tempD2, features = rev(c("Drd1", "Drd2", "Drd3", "Drd4", "Drd5")), 
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

ex<-cowplot::plot_grid(ex.D1, ex.D2, align = "hv", axis = "tblr", ncol=2)
ggsave(ex, filename = "/Users/agodino/Desktop/ex_yellow_DotPlot.pdf", device = cairo_pdf, width = 9, height = 9, units = "cm") 
