## Seurat new snRNA-Seq Marine

## Integration all in one

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

wd <- "/Users/arthur/Library/CloudStorage/GoogleDrive-arthur.godino@icahn.mssm.edu/My Drive/Collaborations/Marine/Marine_snRNAseq_2/Cell Ranger Output/"
cellranger_import <- function(wd, pattern = "[[:digit:]].tar.gz") {
  samplenames <- list.files(wd, pattern = pattern) %>% gsub(".tar.gz", "", .)
  l <- list()
  for (i in samplenames) {
    print(paste0("Importing ", i))
    l[[i]] <- Read10X(data.dir = paste0(wd, i, "/filtered_feature_bc_matrix/")) %>%
      CreateSeuratObject(counts = ., project = i, min.cells = 0, min.features = 0)
  }
  return(l)
}
raw.list <- cellranger_import(wd)

# Merging
raw <- merge(raw.list[[1]], y = c(raw.list[[2]],raw.list[[3]],raw.list[[4]],raw.list[[5]],raw.list[[6]],raw.list[[7]],raw.list[[8]]), add.cell.ids = names(raw.list), project = "all")
rm(raw.list)
raw[["percent.mt"]] <- PercentageFeatureSet(raw, pattern = "^mt-")
raw <- subset(raw, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 1 & nCount_RNA > 900)

# Integration
raw.list.filt <- SplitObject(raw, split.by = "ident")
for (i in 1:length(raw.list.filt)) { raw.list.filt[[i]] <- SCTransform(raw.list.filt[[i]], vars.to.regress = c("percent.mt"), verbose = TRUE) }
features <- SelectIntegrationFeatures(object.list = raw.list.filt, nfeatures = 3000)
raw.list.filt <- PrepSCTIntegration(object.list = raw.list.filt, anchor.features = features, verbose =  TRUE)
anchors <- FindIntegrationAnchors(object.list = raw.list.filt, normalization.method = "SCT", anchor.features = features, verbose = TRUE)
integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = TRUE)  #integration of only 3000 anchors

## PCA
DefaultAssay(integrated) <- "integrated"
integrated <- RunPCA(integrated, verbose = TRUE)
ElbowPlot(integrated, ndims = 50)

# Series clustering, ie UMAP run on cluster graph
DefaultAssay(integrated) <- "integrated"
integrated <- FindNeighbors(integrated, dims = 1:16, verbose = TRUE)
integrated <- FindClusters(integrated, resolution = 0.1, graph.name = 'integrated_snn', verbose = TRUE) #, graph.name = 'integrated_snn'
integrated <- RunUMAP(integrated, graph = 'integrated_snn', verbose = TRUE) #, graph = 'integrated_snn'

#integrated <- readRDS("/Users/arthur/Library/CloudStorage/GoogleDrive-arthur.godino@icahn.mssm.edu/My Drive/Marine_snRNAseq_2/Results/integrated_17dims.rds")
# UMAP by group
integrated$Group <- factor(case_when(integrated$orig.ident == "1plus-scRNA-0622" ~ "POS_1P_Test",
                              integrated$orig.ident == "1minus-scRNA-0622" ~ "NEG_1P_Test",
                              integrated$orig.ident == "1plus-scRNA-0627" ~ "POS_1P_HC",
                              integrated$orig.ident == "1minus-scRNA-0627" ~ "NEG_1P_HC",
                              integrated$orig.ident == "2plus-scRNA-0622" ~ "POS_5P_Test",
                              integrated$orig.ident == "2minus-scRNA-0622" ~ "NEG_5P_Test",
                              integrated$orig.ident == "2plus-scRNA-0627" ~ "POS_5P_HC",
                              integrated$orig.ident == "2minus-scRNA-0627" ~ "NEG_5P_HC"))
integrated$sorted <- factor(str_extract(integrated$Group, "^[[:alpha:]]+"))
integrated$pairings <- factor(str_extract(integrated$Group, "(?<=_)[[:alnum:]]+(?=_)"))
integrated$cpp <- factor(str_extract(integrated$Group, "[[:alnum:]]+$"))

UMAPPlot(integrated, label = T) + NoLegend() & theme(aspect.ratio = 1) 

# Combine 2 D2 and 3 D1 (including Grm8) clusters for clarity, as well as poly and oligodendrocytes together
integrated$seurat_clusters <- case_when(integrated$seurat_clusters == 0 ~ 0,
                                        integrated$seurat_clusters == 1 ~ 1,
                                        integrated$seurat_clusters == 2 ~ 0,
                                        integrated$seurat_clusters == 3 ~ 3,
                                        integrated$seurat_clusters == 4 ~ 4,
                                        integrated$seurat_clusters == 5 ~ 5,
                                        integrated$seurat_clusters == 6 ~ 6,
                                        integrated$seurat_clusters == 7 ~ 1,
                                        integrated$seurat_clusters == 8 ~ 0,
                                        integrated$seurat_clusters == 9 ~ 5)
UMAPPlot(integrated, group.by = "seurat_clusters", label = T) + NoLegend() & theme(aspect.ratio = 1) 

VlnPlot(integrated, features = c("nCount_RNA", "nFeature_RNA"),group.by = "seurat_clusters", pt.size = .1)
FeaturePlot(integrated, reduction = "umap", features = c("nCount_RNA", "nFeature_RNA")) & theme(aspect.ratio = 1) 

DefaultAssay(integrated) <- "SCT"
#FeaturePlot(integrated, reduction = "umap", features = c("Drd1"), min.cutoff = 0, max.cutoff = "q90") & theme(aspect.ratio = 1) 
#FeaturePlot(integrated, reduction = "umap", features = c("Drd3"), min.cutoff = 0, max.cutoff = "q90") & theme(aspect.ratio = 1) 
#FeaturePlot(integrated, reduction = "umap", features = c("Drd2"), min.cutoff = 0, max.cutoff = "q90") & theme(aspect.ratio = 1) 
#FeaturePlot(integrated, reduction = "umap", features = c("Grm8"), min.cutoff = 0, max.cutoff = "q90") & theme(aspect.ratio = 1) 
#FeaturePlot(integrated, reduction = "umap", features = c("Elavl2"), min.cutoff = 0, max.cutoff = "q90") & theme(aspect.ratio = 1) 
#FeaturePlot(integrated, reduction = "umap", features = c("Sst"), min.cutoff = 0, max.cutoff = "q90") & theme(aspect.ratio = 1) 
#FeaturePlot(integrated, reduction = "umap", features = c("Kit"), min.cutoff = 0, max.cutoff = "q90") & theme(aspect.ratio = 1) 
#FeaturePlot(integrated, reduction = "umap", features = c("Ppp1r1b"), min.cutoff = 0, max.cutoff = "q90") & theme(aspect.ratio = 1) 
#FeaturePlot(integrated, reduction = "umap", features = c("Foxp2"), min.cutoff = 0, max.cutoff = "q90") & theme(aspect.ratio = 1) 
#FeaturePlot(integrated, reduction = "umap", features = c("Bcl11b"), min.cutoff = 0, max.cutoff = "q90") & theme(aspect.ratio = 1) 
#FeaturePlot(integrated, reduction = "umap", features = c("Mbp"), min.cutoff = 0, max.cutoff = "q90") & theme(aspect.ratio = 1)
#FeaturePlot(integrated, reduction = "umap", features = c("Gja1"), min.cutoff = 0, max.cutoff = "q90") & theme(aspect.ratio = 1)
#FeaturePlot(integrated, reduction = "umap", features = c("Pdgfra"), min.cutoff = 0, max.cutoff = "q90") & theme(aspect.ratio = 1)
#FeaturePlot(integrated, reduction = "umap", features = c("Arhgap15"), min.cutoff = 0, max.cutoff = "q90") & theme(aspect.ratio = 1)
#FeaturePlot(integrated, reduction = "umap", features = c("Ebf1"), min.cutoff = 0, max.cutoff = "q90") & theme(aspect.ratio = 1)
#FeaturePlot(integrated, reduction = "umap", features = c("Slc17a7"), min.cutoff = 0, max.cutoff = "q90") & theme(aspect.ratio = 1)
#FeaturePlot(integrated, reduction = "umap", features = c("Tac2"), min.cutoff = 0, max.cutoff = "q90") & theme(aspect.ratio = 1)
#FeaturePlot(integrated, reduction = "umap", features = c("Penk"), min.cutoff = 0, max.cutoff = "q90") & theme(aspect.ratio = 1)
#FeaturePlot(integrated, reduction = "umap", features = c("Pdyn"), min.cutoff = 0, max.cutoff = "q90") & theme(aspect.ratio = 1)
#FeaturePlot(integrated, reduction = "umap", features = c("Pdyn"), min.cutoff = 0, max.cutoff = "q90") & theme(aspect.ratio = 1)

integrated$comb.clusters <- case_when(integrated$seurat_clusters == 0 ~ "D1-MSNs",
                                      integrated$seurat_clusters == 1 ~ "D2-MSNs",
                                      integrated$seurat_clusters == 3 ~ "Astrocytes",
                                      integrated$seurat_clusters == 4 ~ "Interneurons",
                                      integrated$seurat_clusters == 5 ~ "Oligodendrocytes",
                                      integrated$seurat_clusters == 6 ~ "Glutamatergic")
integrated$comb.clusters <- factor(integrated$comb.clusters, levels = rev(c("D1-MSNs", "D2-MSNs", "Interneurons", "Astrocytes", "Oligodendrocytes", "Glutamatergic")))

integrated <- SetIdent(integrated, value='comb.clusters')
umap <- UMAPPlot(integrated, label = T, pt.size = .1, label.size = 5/2.8, cols = setNames(scales::hue_pal()(6), levels(integrated$comb.clusters))) + 
  theme_classic(base_size = 5, base_family = "Arial") +
  scale_x_continuous(limits = c(-14,21), expand = c(0,0)) +
  scale_y_continuous(limits = c(-21,14), expand = c(0,0)) +
  theme(plot.background = element_rect(fill="transparent", colour=NA),
        panel.background = element_rect(fill="transparent", colour=NA),
        plot.margin = margin(.5,.5,.5,.5, "char"),
        aspect.ratio = 1, axis.text = element_text(color = "black"),
        legend.position = "none")
ggsave(umap, filename = "/Users/arthur/Desktop/clustering_UMAP.pdf", device = cairo_pdf, width = 15, height = 15, units = "cm") 

# Raw marker genes
# Normalize for sequencing depth in SCT array
integrated <- PrepSCTFindMarkers(integrated)

Idents(integrated) <- "comb.clusters"
raw.markers <- FindAllMarkers(integrated, only.pos = TRUE, min.pct = .3, logfc.threshold = log2(1.15), assay = "SCT", test.use = "LR", latent.vars = "orig.ident")

export_genelist <- function(rm, f = "/Users/arthur/Desktop/raw.markers_integrated_dim16_res0.1_lr.xlsx") {
t <- rm %>% split(., .$cluster)
wb <- createWorkbook()                                                    
sheetnames <- names(t)
sheets <- lapply(sheetnames, createSheet, wb = wb)
void <- Map(addDataFrame, t, sheets)
saveWorkbook(wb, file = f) }
export_genelist(raw.markers, "/Users/arthur/Desktop/raw.markers_integrated_lr_SCT.xlsx")

# Dot plot of marker genes
dp.markers <- DotPlot(integrated, features = c("Rbfox3", "Gad1", "Ppp1r1b", "Foxp2", "Bcl11b", "Drd1", "Tac1", "Pdyn",  "Drd2", "Adora2a", "Penk", "Grm8", "Elavl2", "Kit", "Sst", "Gja1", "Gli3", "Mbp", "Mog", "Slc17a7"), 
                      cluster.idents = F,
                      group.by = "comb.clusters", assay = "SCT", dot.scale = 6, split.by = "comb.clusters",
                      cols = setNames(scales::hue_pal()(6), levels(integrated$comb.clusters)),
                      dot.min = .01, col.min = 0, col.max = 1) + 
  scale_y_discrete(labels = levels(integrated$comb.clusters)) +
  theme_classic(base_family = "Arial", base_size = 5) +
  theme(plot.background = element_rect(fill="transparent", colour=NA),
        panel.background = element_rect(fill="transparent", colour=NA),
        aspect.ratio = 7/20, axis.title = element_blank(), axis.text = element_text(color = "black"), axis.ticks = element_line(color="black"),
        axis.text.x = element_text(face = "italic", angle = 45, vjust = 0.8, hjust = 1),
        legend.key.size = unit(.5, "line"), legend.position = "bottom", legend.box = "vertical")
ggsave(dp.markers, filename = "/Users/arthur/Desktop/clustering_DotPlot.pdf", device = cairo_pdf, width = 9, height = 6, units = "cm") 

# Annotate metadata
integrated$cluster_sorted <- paste(integrated$comb.clusters, integrated$sorted, sep = "_")
integrated$cluster_pairings <- paste(integrated$comb.clusters, integrated$pairings, sep = "_")
integrated$cluster_cpp <- paste(integrated$comb.clusters, integrated$cpp, sep = "_")
integrated$cluster_pairings_cpp <- paste(integrated$comb.clusters, integrated$pairings, integrated$cpp, sep = "_")
integrated$cluster_sorted_pairings <- paste(integrated$comb.clusters, integrated$sorted, integrated$pairings, sep = "_")
integrated$cluster_sorted_pairings_cpp <- paste(integrated$comb.clusters, integrated$sorted, integrated$pairings, integrated$cpp, sep = "_")
integrated$sorted_pairings <- paste(integrated$sorted, integrated$pairings, sep = "_")
integrated$sorted_pairings_cpp <- paste(integrated$sorted, integrated$pairings, integrated$cpp, sep = "_")
integrated$Arc.expr <- ifelse(GetAssayData(object = integrated, assay = "SCT", slot = "data")["Arc",] > 0, "YES", "NO")
integrated$pairings_cpp_Arc <- paste(integrated$pairings, integrated$cpp, integrated$Arc.expr, sep = "_")
integrated$react <- factor(paste(integrated$sorted, integrated$Arc.expr, sep = "_"), levels = c("NEG_NO", "POS_NO", "NEG_YES", "POS_YES"))
integrated$cluster_pairings_cpp_sorted_Arc <- paste(integrated$comb.clusters, integrated$pairings, integrated$cpp, integrated$sorted, integrated$Arc.expr, sep = "_")
integrated$cluster_pairings_cpp_Arc <- paste(integrated$comb.clusters, integrated$pairings, integrated$cpp, integrated$Arc.expr, sep = "_")
integrated$cluster_cpp_react <- paste(integrated$comb.clusters, integrated$cpp, integrated$react, sep = "_")
integrated$cluster_Arc <- paste(integrated$comb.clusters, integrated$Arc.expr, sep = "_")
integrated$pairings_cpp_sorted_Arc <- factor(paste(integrated$pairings, integrated$cpp, integrated$sorted, integrated$Arc.expr, sep = "_"), levels = c("5P_Test_POS_YES", "1P_Test_POS_YES", "5P_Test_NEG_YES", "1P_Test_NEG_YES", "5P_Test_POS_NO", "1P_Test_POS_NO", "5P_Test_NEG_NO", "1P_Test_NEG_NO",
                                                                                                                                                       "5P_HC_POS_YES", "1P_HC_POS_YES", "5P_HC_NEG_YES", "1P_HC_NEG_YES", "5P_HC_POS_NO", "1P_HC_POS_NO", "5P_HC_NEG_NO", "1P_HC_NEG_NO"))
integrated$pairings_sorted_Arc <- paste(integrated$pairings, integrated$sorted, integrated$Arc.expr, sep = "_")

# Remove all glutamatergic cells
Idents(integrated) <- "comb.clusters"
integrated <- subset(integrated, idents = "Glutamatergic", invert = T)
integrated <- PrepSCTFindMarkers(integrated)   # renormalize for library size now that glutamatergic are removed

saveRDS(integrated, "/Users/arthur/Desktop/Marine_integrated.rds")
