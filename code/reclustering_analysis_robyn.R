#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 08/21/2024
# Description: Comparative analysis of immune and tumor compartments between
# TGFb high/low
#==============================================================================#

#==============================================================================#
# Libraries, helper functions, and environment variables
#==============================================================================#

library(Seurat)
library(dplyr)
library(tidyverse)
library(ggrepel)
library(patchwork)
library(scProportionTest)
library(escape)

source("/home/hnatri/CART_TGFb/code/CART_plot_functions.R")
source("/home/hnatri/CART_TGFb/code/13384_tumor_ms_themes.R")
source("/home/hnatri/CART_TGFb/code/utilities.R")

setwd("/home/hnatri/CART_TGFb/")

set.seed(1234)
options(future.globals.maxSize = 30000 * 1024^2)
reduction <- "integrated_sct_umap"

# A function for reclustering Seurat data
recluster_seurat <- function(seurat_object){
  DefaultAssay(seurat_object) <- "RNA"
  seurat_object <- SCTransform(seurat_object,
                               vars.to.regress = c("percent.mt"))
  seurat_object <- RunPCA(seurat_object,
                          reduction.name = "pca",
                          verbose = F)
  pcs <- get_pcs(seurat_object, reduction_name = "pca")
  message(pcs)
  seurat_object <- RunUMAP(seurat_object,
                           reduction = "pca",
                           reduction.name = "umap",
                           dims = 1:min(pcs),
                           return.model = TRUE)
  seurat_object <- FindNeighbors(seurat_object,
                                 reduction = "pca",
                                 dims = 1:min(pcs),
                                 graph.name = c("snn",
                                                "snn"))
  seurat_object <- FindClusters(seurat_object,
                                resolution = c(0.1,0.2,0.3,0.5,0.8,1),
                                graph.name = "snn")
  
  return(seurat_object)
}


#==============================================================================#
# Import data
#==============================================================================#

tumors <- readRDS("/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/tumor_integrated_UPN109pre_noUPN208_soupX_snn_metadata_no4_7_DoubletFinder.rds")
immune_fibro <- readRDS("/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/13384_tumors_immune_fibro_GSEA.rds")

# Inspecting clustering and cell type annotations
DimPlot(tumors,
        group.by = "celltype",
        cols = tumor_celltype_col,
        reduction = reduction,
        label = T,
        label.box = T,
        label.size = 3,
        repel = T,
        raster = T,
        raster.dpi = c(1024, 1024),
        pt.size = 3) +
  ggtitle("") +
  theme_classic() +
  NoLegend() +
  NoAxes() +
  coord_fixed(1)

DimPlot(immune_fibro,
        group.by = "celltype",
        cols = immune_fibro_celltype_col,
        reduction = "umap",
        label = T,
        label.box = T,
        label.size = 3,
        repel = T,
        raster = T,
        raster.dpi = c(1024, 1024),
        pt.size = 3) +
  ggtitle("") +
  theme_classic() +
  NoLegend() +
  NoAxes() +
  coord_fixed(1)

# CD45 expression
FeaturePlot(tumors,
            features = c("PTPRC"),
            cols = c("gray89", "tomato3"),
            ncol = 1,
            raster = T,
            raster.dpi = c(2048, 1024),
            pt.size = 4) &
  coord_fixed(ratio = 1) &
  theme_classic() &
  NoAxes()

immune_robyn <- readRDS("/scratch/hnatri/CART/fromRobyn/IMMUNEminus618_TGFbELISA+TGFbELISA2.rds")
#tumor_robyn <- readRDS("/scratch/hnatri/CART/fromRobyn/TumorPLUS_29+0.5.rds")

# Fetching TGFb status for each sample from Robyn's object
tgfb_status <- immune_robyn@meta.data %>%
  as.data.frame() %>%
  dplyr::select("UPN", "TGFbELISA", "TGFbELISA2") %>%
  distinct()

# Numbers of samples in each group
table(tgfb_status$TGFbELISA)
table(tgfb_status$TGFbELISA2)

# Adding TGFb status to the annotated data
immune_fibro$TGFbELISA <- mapvalues(x = immune_fibro$UPN,
                                    from = tgfb_status$UPN,
                                    to = tgfb_status$TGFbELISA)
tumors$TGFbELISA <- mapvalues(x = tumors$UPN,
                              from = tgfb_status$UPN,
                              to = tgfb_status$TGFbELISA)

immune_fibro$TGFbELISA2 <- mapvalues(x = immune_fibro$UPN,
                                     from = tgfb_status$UPN,
                                     to = tgfb_status$TGFbELISA2)
tumors$TGFbELISA2 <- mapvalues(x = tumors$UPN,
                               from = tgfb_status$UPN,
                               to = tgfb_status$TGFbELISA2)


#==============================================================================#
# Splitting and reclustering
#==============================================================================#

# Immune clusters
immune_seurat <- subset(immune_fibro, subset = celltype %in% c(paste0("M", seq(1, 9)),
                                                               paste0("L", seq(1, 10)),
                                                               "B1", "N1"))
# Tumor clusters
tumor_seurat <- subset(tumors, subset = celltype %in% c("Olig1", "Olig2", "Fibroblast",
                                                        paste0("Tumor", seq(1, 6))))

# Reclustering
immune_seurat <- recluster_seurat(immune_seurat)
tumor_seurat <- recluster_seurat(tumor_seurat)

# Saving
#saveRDS(immune_seurat, "/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/immune_reclustered_TGFb.rds")
#saveRDS(tumor_seurat, "/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/tumor_reclustered_TGFb.rds")

# Inspecting the clustering
#immune_seurat <- readRDS("/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/immune_reclustered_TGFb.rds")
#tumor_seurat <- readRDS("/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/tumor_reclustered_TGFb.rds")

DimPlot(tumor_seurat,
        group.by = "celltype",
        cols = tumor_celltype_col,
        reduction = "umap",
        label = T,
        label.box = T,
        label.size = 3,
        repel = T,
        raster = T,
        raster.dpi = c(1024, 1024),
        pt.size = 3) +
  ggtitle("") +
  theme_classic() +
  NoLegend() +
  NoAxes() +
  coord_fixed(1)

DimPlot(immune_seurat,
        group.by = "celltype",
        cols = immune_fibro_celltype_col,
        reduction = "umap",
        label = T,
        label.box = T,
        label.size = 3,
        repel = T,
        raster = T,
        raster.dpi = c(1024, 1024),
        pt.size = 3) +
  ggtitle("") +
  theme_classic() +
  NoLegend() +
  NoAxes() +
  coord_fixed(1)

FeaturePlot(tumor_seurat,
            features = c("PTPRC"),
            cols = c("gray89", "tomato3"),
            ncol = 1,
            raster = T,
            raster.dpi = c(2048, 1024),
            pt.size = 4) &
  coord_fixed(ratio = 1) &
  theme_classic() &
  NoAxes()

FeaturePlot(immune_seurat,
            features = c("PTPRC"),
            cols = c("gray89", "tomato3"),
            #order = T,
            #keep.scale = "all",
            ncol = 1,
            raster = T,
            raster.dpi = c(2048, 1024),
            pt.size = 4) &
  coord_fixed(ratio = 1) &
  theme_classic() &
  NoAxes()

DotPlot(immune_seurat,
        features = c("PTPRC", "CD4", "CD8A", "ITGAX", "CD33"),
        group.by = "celltype",
        cols = c("azure", "tomato3")) +
  RotatedAxis()

DotPlot(tumor_seurat,
        features = c("PTPRC", "CD4", "CD8A", "ITGAX", "CD33"),
        group.by = "celltype",
        cols = c("azure", "tomato3")) +
  RotatedAxis()

#==============================================================================#
# Reassigning Tumor6
#==============================================================================#

# Subsetting, merging with the immune object, and reclustering
tumor6 <- subset(tumor_seurat, subset = celltype == "Tumor6")
immune_tumor6_seurat <- merge(immune_seurat, tumor6)
immune_tumor6_seurat <- recluster_seurat(immune_tumor6_seurat)

#saveRDS(immune_tumor6_seurat, "/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/immune_tumor6_seurat.rds")

# All clusters except Tumor6 for comparison
tumor_wout6 <- subset(tumor_seurat, subset = celltype == "Tumor6", invert = T)

##saveRDS(tumor_wout6, "/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/tumor_wout6_seurat.rds")
#immune_tumor6_seurat <- readRDS("/tgen_labs/banovich/BCTCSF/Heini/13384_Tumor/immune_tumor6_seurat.rds")

# Plotting immune cell types with Tumor6
DimPlot(immune_tumor6_seurat,
        group.by = "celltype",
        cols = c(tumor_celltype_col, immune_fibro_celltype_col),
        reduction = "umap",
        label = T,
        label.box = T,
        label.size = 3,
        repel = T,
        raster = T,
        raster.dpi = c(1024, 1024),
        pt.size = 3) +
  ggtitle("") +
  theme_classic() +
  NoLegend() +
  NoAxes() +
  coord_fixed(1)

FeaturePlot(immune_tumor6_seurat,
            features = c("PTPRC"),
            cols = c("gray89", "tomato3"),
            #order = T,
            #keep.scale = "all",
            ncol = 1,
            raster = T,
            raster.dpi = c(2048, 1024),
            pt.size = 4) &
  coord_fixed(ratio = 1) &
  theme_classic() &
  NoAxes()

#==============================================================================#
# Cell type proportions between TGFb high and low
#==============================================================================#

# Analysis-ready object without Tumor6
DimPlot(immune_seurat,
        group.by = "celltype",
        cols = immune_fibro_celltype_col,
        reduction = "umap",
        label = T,
        label.box = T,
        label.size = 3,
        repel = T,
        raster = T,
        raster.dpi = c(1024, 1024),
        pt.size = 3) +
  ggtitle("") +
  theme_classic() +
  NoLegend() +
  NoAxes() +
  coord_fixed(1)

# Cell type proportions by TGFbELISA2
table(immune_seurat$TGFbELISA2)
immune_seurat_tgfb <- subset(immune_seurat, subset = TGFbELISA2 %in% c("HighTGFb", "LowTGFb"))

table(immune_seurat_tgfb$TGFbELISA2, immune_seurat_tgfb$celltype)

create_clusterpropplot(seurat_object = immune_seurat_tgfb,
                       group_var = "TGFbELISA2",
                       group2 = "HighTGFb",
                       group1 = "LowTGFb",
                       plot_var = "celltype",
                       plot_colors = immune_fibro_celltype_col,
                       var_names = c("TGFb low", "TGFb high"),
                       legend_title = "")

# Using scProportionTest
prop_test <- sc_utils(immune_seurat_tgfb)

# Permutation testing and bootstrapping
prop_test <- permutation_test(
  prop_test, cluster_identity = "celltype",
  sample_1 = "HighTGFb", sample_2 = "LowTGFb",
  sample_identity = "TGFbELISA2")

# HighTGFb gets negative values
perm_plot <- permutation_plot(prop_test)

perm_plot + scale_colour_manual(values = c("tomato", "azure2")) + NoLegend()
