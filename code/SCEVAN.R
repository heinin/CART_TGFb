

#==============================================================================
# Libraries and environment variables
#==============================================================================

#library(devtools)
#install_github("miccec/yaGST")
#install_github("AntonioDeFalco/SCEVAN")
#
#packageurl <- "http://cran.r-project.org/src/contrib/tidytree_0.4.6.tar.gz"
#install.packages(packageurl, repos=NULL, type="source")

library(Seurat)
library(SCEVAN)
library(tidyverse)
library(RColorBrewer)

set.seed(1234)
options(future.globals.maxSize = 300000 * 1024^2)

setwd("/scratch/hnatri/SCEVAN")

#==============================================================================
# Load data
#==============================================================================

seurat_data <- readRDS("/scratch/hnatri/CART/TGFb_scRNAseq_integrated_seurat.rds")
DefaultAssay(seurat_data) <- "RNA"

unique(seurat_data$Sample)
#seurat_data <- subset(seurat_data, subset = Sample == "UCPIN_KO", invert = T)

unique(seurat_data$KO_WT)
#seurat_data <- subset(seurat_data, subset = KO_WT == "KO")

####count_mtx <- LayerData(seurat_data, assay = "RNA", layer = "counts.1")
####
#####==============================================================================
##### Run SCEVAN on whole dataset
#####==============================================================================
####
##### Using WT cells as reference
#####norm_cells <- seurat_data@meta.data %>%
#####  filter(KO_WT == "WT") %>% rownames()
####
####results <- pipelineCNA(count_mtx,
####                       #norm_cell = norm_cells,
####                       SUBCLONES = T,
####                       sample = "all_samples_KO")
####
####saveRDS(results, "/scratch/hnatri/CART/SCEVANres.rds")
#####results <- readRDS("/scratch/hnatri/CART/SCEVANres.rds")
####
##### Add SCEVAN info to an existing Seurat object
####seurat_data <- Seurat::AddMetaData(seurat_data, metadata = results)
####
####saveRDS(seurat_data, "/scratch/hnatri/CART/TGFb_scRNAseq_integrated_seurat_SCEVAN_KOonly_bclones.rds")
####
####q(save = "no")
####seurat_data <- readRDS("/scratch/hnatri/CART/TGFb_scRNAseq_integrated_seurat_SCEVAN_KOonly_bclones.rds")
####
####head(seurat_data@meta.data)

#==============================================================================
# Per sample analysis
#==============================================================================

HD680 <- subset(seurat_data, subset = CellLine == "HD680")
HD562 <- subset(seurat_data, subset = CellLine == "HD562")
HD624 <- subset(seurat_data, subset = CellLine == "HD624")
UCPIN <- subset(seurat_data, subset = CellLine == "UCPIN")

HD680ko <- subset(HD680, subset = KO_WT == "KO")
HD562ko <- subset(HD562, subset = KO_WT == "KO")
HD624ko <- subset(HD624, subset = KO_WT == "KO")
UCPINko <- subset(HD624, subset = KO_WT == "KO")

HD680wt <- subset(HD680, subset = KO_WT == "WT")
HD562wt <- subset(HD562, subset = KO_WT == "WT")
HD624wt <- subset(HD624, subset = KO_WT == "WT")

HD680_ko <- LayerData(HD680ko, assay = "RNA", layer = "counts.1")
HD562_ko <- LayerData(HD562ko, assay = "RNA", layer = "counts.1")
HD624_ko <- LayerData(HD624ko, assay = "RNA", layer = "counts.1")
UCPIN_ko <- LayerData(HD624ko, assay = "RNA", layer = "counts.1")

HD680_wt <- LayerData(HD680wt, assay = "RNA", layer = "counts.1")
HD562_wt <- LayerData(HD562wt, assay = "RNA", layer = "counts.1")
HD624_wt <- LayerData(HD624wt, assay = "RNA", layer = "counts.1")


results <- pipelineCNA(UCPIN_ko,
                       SUBCLONES = T,
                       sample = "UCPIN_KO")

saveRDS(results, "/scratch/hnatri/CART/TGFbRa2KOWT_multisample_UCPIN_KO_res.rds")

q(save = "no")

results <- pipelineCNA(HD680_wt,
                       SUBCLONES = T,
                       sample = "HD680_WT")

saveRDS(results, "/scratch/hnatri/CART/TGFbRa2KOWT_multisample_HD680_WT_res.rds")

#listCountMtx <- list(HD562_ko = HD562_ko,
#                     HD562_wt = HD562_wt)

#results <- SCEVAN::multiSampleComparisonClonalCN(listCountMtx,
#                                                 analysisName = "TGFbRa2KO_HD562",
#                                                 organism = "human",
#                                                 par_cores = 5)
results <- pipelineCNA(HD562_wt,
                       SUBCLONES = T,
                       sample = "HD562_WT")

saveRDS(results, "/scratch/hnatri/CART/TGFbRa2KOWT_multisample_HD562_WT_res.rds")

#listCountMtx <- list(HD624_ko = HD624_ko,
#                     HD624_wt = HD624_wt)
#
#results <- SCEVAN::multiSampleComparisonClonalCN(listCountMtx,
#                                                 analysisName = "TGFbRa2KO_HD624",
#                                                 organism = "human",
#                                                 par_cores = 5)

results <- pipelineCNA(HD624_wt,
                       SUBCLONES = T,
                       sample = "HD624_WT")

saveRDS(results, "/scratch/hnatri/CART/TGFbRa2KOWT_multisample_HD624_WT_res.rds")

q(save = "no")

#==============================================================================
# Without UCPIN
#==============================================================================

seurat_data <- readRDS("/scratch/hnatri/CART/TGFb_scRNAseq_integrated_seurat.rds")
DefaultAssay(seurat_data) <- "RNA"

unique(seurat_data$Sample)
seurat_data <- subset(seurat_data, subset = Sample == "UCPIN_KO", invert = T)

ko <- subset(seurat_data, subset = KO_WT == "KO")
wt <- subset(seurat_data, subset = KO_WT == "WT")

ko <- LayerData(ko, assay = "RNA", layer = "counts.1")
wt <- LayerData(wt, assay = "RNA", layer = "counts.1")

listCountMtx <- list(ko = ko,
                     wt = wt)

results <- SCEVAN::multiSampleComparisonClonalCN(listCountMtx,
                                                 analysisName = "TGFbRa2KOWT_multisample_woutUCPIN",
                                                 organism = "human",
                                                 par_cores = 5)

saveRDS(results, "/scratch/hnatri/CART/TGFbRa2KOWT_multisample_woutUCPIN_res.rds")

