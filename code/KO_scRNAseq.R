#==============================================================================
# Author(s) : Heini M. Natri, hnatri@tgen.org
# Date: 07/26/2024
# Description: TGFbRa2KO, CAR/TGFb sc-RNAseq data analysis
#==============================================================================

suppressMessages({library(Seurat)
  library(ggplot2)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(googlesheets4)})

#==============================================================================
# Helper functions
#==============================================================================

source("/home/hnatri/13384_CART/CART_plot_functions.R")
source("/home/hnatri/SingleCellBestPractices/scripts/preprocessing_qc_module.R")
source("/home/hnatri/SingleCellBestPractices/scripts/integration_module.R")

#==============================================================================
# Environment variables
#==============================================================================

set.seed(1234)
options(future.globals.maxSize = 30000 * 1024^2)
reduction <- "integratedSCTumap"

#==============================================================================
# Import data
#==============================================================================

# Metadata
gs4_deauth()
TGFbRa2_KO_metadata  <- gs4_get("https://docs.google.com/spreadsheets/d/1T3oB6dS_rbrijBWKAlPUCOWNAFc2Zz1yHDsUJRSdpOg/edit?usp=sharing")
sheet_names(TGFbRa2_KO_metadata)
scRNAseq_samples <- read_sheet(TGFbRa2_KO_metadata, sheet = "scRNAseq")
scRNAseq_samples$Sample_Name_Run <- paste0(scRNAseq_samples$Sample_Name, "_",
                                           scRNAseq_samples$Run)

# For SoupX
path_list <- list("Run1" = "/tgen_labs/banovich/SingleCell/CellRanger/8_0_1/Projects/BCTCSF/CITE/F07624-GEX_F07628-CAR_F07626-FB/",
                  "Run2" = "/tgen_labs/banovich/SingleCell/CellRanger/8_0_1/Projects/BCTCSF/CITE/F07625-GEX_F07629-CAR_F07627-FB/")

# Importing  and processing data
# Two GEM runs, same samples, same loading
# F07624-GEX_F07626-FB_F07628-CAR
# F07625-GEX_F07627-FB_F07629-CAR
# F07630-VDJ
# F07631-VDJ

# Cell hashing antibody names
hash_antibodies <- c("TotalSeqC0257_Hashtag7",
                     "TotalSeqC0258_Hashtag8",
                     "TotalSeqC0259_Hashtag9",
                     "TotalSeqC0260_Hashtag10",
                     "TotalSeqC0262_Hashtag12",
                     "TotalSeqC0263_Hashtag13",
                     "TotalSeqC0264_Hashtag14")

# Demultiplexing
seurat_data_list <- prep_seurat_list_multiplexed(metadata = scRNAseq_samples,
                                                 batch_ID = "Run",
                                                 cellRanger_path = "CellRanger_path",
                                                 cell_ID_prefix = "Sample_Name_Run",
                                                 CellHashing_Ab = "Hash_ID",
                                                 Hash_Abs = hash_antibodies,
                                                 pos_quant = 0.99)

table(seurat_data_list[[1]]$Sample_Name)
table(seurat_data_list[[2]]$Sample_Name)

# Merging for plotting
seurat_data_merged <- merge(seurat_data_list[[1]], seurat_data_list[[2]])
seurat_data_merged <- JoinLayers(seurat_data_merged, assay = "RNA")

table(seurat_data_merged$Sample_Name) %>% as.data.frame() %>%
  ggplot(aes(x = Var1, y = Freq)) +
  geom_bar(position="dodge", stat="identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab("Sample") +
  ylab("# cells")

table(seurat_data_merged$hash.ID) %>% as.data.frame() %>%
  ggplot(aes(x = Var1, y = Freq)) +
  geom_bar(position="dodge", stat="identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab("Sample") +
  ylab("# cells")

DefaultAssay(seurat_data_merged) <- "RNA"

#==============================================================================
# Plotting CAR enrichment
#==============================================================================

# The name of the construct is "CLTX-Construct"
layer_data <- LayerData(seurat_data_merged,
                        assay = "RNA",
                        layer = "counts")
layer_data <- as.data.frame(t(layer_data))
layer_data$Sample_Name <- plyr::mapvalues(x = rownames(layer_data),
                                          from = rownames(seurat_data_merged@meta.data),
                                          to = as.character(seurat_data_merged@meta.data$Sample_Name))
layer_data$hash.ID <- plyr::mapvalues(x = rownames(layer_data),
                                      from = rownames(seurat_data_merged@meta.data),
                                      to = as.character(seurat_data_merged@meta.data$hash.ID))

layer_data %>% ggplot(aes(x = Sample_Name, y = log2(`IL13OP`))) +
  geom_violin() +
  theme_classic()

layer_data %>% ggplot(aes(x = hash.ID, y = log2(`IL13OP`))) +
  geom_violin() +
  theme_classic()

layer_data %>% ggplot(aes(x = log2(`IL13OP`))) +
  geom_histogram() +
  facet_wrap(~Sample_Name, scales = "free") +
  theme_classic()

layer_data %>% ggplot(aes(x = log2(`IL13OP`))) +
  geom_histogram() +
  facet_wrap(~hash.ID, scales = "free") +
  theme_classic()

# What % of CAR library mapped to the construct?
dim(layer_data)
ncol(layer_data)

unique(layer_data$hash.ID)

layer_data %>%
  dplyr::mutate(sum = rowSums(across(1:(ncol(layer_data)-2)))) %>%
  #dplyr::summarise(CAR_counts = sum(`CLTX-Construct`),
  #                 lib_size = sum(1:(ncol(layer_data)-2)),
  #                 pct_CAR = (CAR_counts/lib_size)*100) %>%
  group_by(hash.ID) %>%
  dplyr::summarise(CAR_counts = sum(`IL13OP`),
                   lib_size = sum(sum),
                   pct_CAR = (CAR_counts/lib_size)*100) %>%
  ungroup()

# Inspecting CAR enrichment libraries only

#==============================================================================
# Filtering, SoupX, integration, dimensionality reduction, and clustering
#==============================================================================

# Filtering to retain singlets
seurat_singlets_list <- lapply(seurat_data_list, function(xx){
  xx <- subset(xx, subset = HTO_classification.global == "Singlet")
  
  xx
})
names(seurat_singlets_list) <- names(seurat_data_list)

# Saving demultiplexed data in 10x format
for (i in 1:length(seurat_singlets_list)){
  print(paste0("Converting sample ", names(seurat_singlets_list[i])))
  obj.sub <- seurat_singlets_list[[i]] 
  
  counts <- LayerData(obj.sub, assay = "RNA", layer = "counts")
  DropletUtils::write10xCounts(path = paste0("/scratch/hnatri/CART/SoupX/demultiplexed_", names(seurat_singlets_list[i])), 
                               x = counts, 
                               barcodes = colnames(counts), # cell names
                               gene.id = rownames(counts), # Gene identifiers, one per row of X
                               type = "sparse", version="3", overwrite = T)
}

# Running SoupX
seurat_soupx_list <- sapply(names(seurat_singlets_list), function(xx){
  print(xx)
  
  # Read in count and droplet data
  # Converted demultiplexed counts
  message("Reading converted demultiplexed counts")
  d10x_toc <- Read10X(paste0("/scratch/hnatri/CART/SoupX/demultiplexed_", xx))
  
  #batch_id <- seurat_singlets_list[[xx]]$Run
  
  # Need to read in batch specific empty droplet file
  message("Reading empty droplet data")
  d10x_tod <- Read10X(paste0(path_list[[xx]], "/outs/raw_feature_bc_matrix/"))
  
  if(length(names(d10x_tod))>1){
    d10x_tod_gex <- d10x_tod[[1]]
  } else {
    d10x_tod_gex <- d10x_tod
  }
  
  colnames(d10x_tod_gex) <- paste0(xx, "_", colnames(d10x_tod_gex))
  
  # Some batches only have features that are expressed in at least one cell;
  # need to fix feature order
  #rownames(d10x_toc) <- gsub("_", "-", rownames(d10x_toc))
  #rownames(d10x_tod_gex) <- gsub("_", "-", rownames(d10x_tod_gex))
  d10x_tod_gex <- d10x_tod_gex[rownames(d10x_toc),]
  
  # Run SoupX
  sc <- SoupChannel(d10x_tod_gex, d10x_toc, calcSoupProfile = FALSE) 
  
  sc <- estimateSoup(sc)
  toc_seu <- CreateSeuratObject(d10x_toc)
  toc_seu <- SCTransform(toc_seu, vst.flavor = "v2")
  if(ncol(toc_seu)<50){
    toc_seu <- RunPCA(toc_seu, npcs = ncol(toc_seu)-1)
  } else{
    toc_seu <- RunPCA(toc_seu)
  }
  toc_seu <- RunUMAP(toc_seu, dims = 1:15)
  toc_seu <- FindNeighbors(toc_seu, dims = 1:15)
  toc_seu <- FindClusters(toc_seu, resolution = 1)
  
  # Add meta data to soupX object
  sc <- setClusters(sc, setNames(toc_seu$seurat_clusters, rownames(toc_seu@meta.data)))
  # Estimate contamination (automated method)
  message(paste0("Getting autoEstCont for: ", i))
  sc <- autoEstCont(sc, tfidfMin = 0.6, soupQuantile = 0.7, forceAccept = TRUE, doPlot = F) 
  out <- adjustCounts(sc)
  
  # Create Seurat object using corrected data
  d10x_seu <- CreateSeuratObject(out, assay = "SoupX_RNA")
  d10x_seu[["RNA"]] <- toc_seu@assays[["RNA"]]
  d10x_seu <- PercentageFeatureSet(d10x_seu, pattern = "^MT-", col.name = "percent.mt_RNA", assay = "RNA")
  d10x_seu <- PercentageFeatureSet(d10x_seu, pattern = "^RP[SL]|^MRP[SL]", col.name = "percent.ribo_RNA", assay = "RNA")
  d10x_seu <- PercentageFeatureSet(d10x_seu, pattern = "^MT-", col.name = "percent.mt_SoupX_RNA", assay = "SoupX_RNA")
  d10x_seu <- PercentageFeatureSet(d10x_seu, pattern = "^RP[SL]|^MRP[SL]", col.name = "percent.ribo_SoupX_RNA", assay = "SoupX_RNA")
  
  # Add sample metadata
  d10x_seu$Sample <- paste0(sapply(strsplit(colnames(d10x_seu), "_"), `[`, 1), "_",
                            sapply(strsplit(colnames(d10x_seu), "_"), `[`, 2))
  
  d10x_seu$CellLine <-sapply(strsplit(colnames(d10x_seu), "_"), `[`, 1)
  d10x_seu$KO_WT <-sapply(strsplit(colnames(d10x_seu), "_"), `[`, 2)
  d10x_seu$Batch <- sapply(strsplit(colnames(d10x_seu), "_"), `[`, 3)
  d10x_seu$cellname <- colnames(d10x_seu)
  
  #d10x_seu <- RenameCells(d10x_seu,
  #                        new.names = paste0(xx, "_", colnames(d10x_seu)))
  
  d10x_seu
})

names(seurat_soupx_list) <- names(seurat_singlets_list)

# Merge
seurat_merged <- merge(x = seurat_soupx_list[[1]], y = seurat_soupx_list[[2]])

saveRDS(seurat_soupx_list, "/scratch/hnatri/CART/SoupX/seurat_soupx_list.rds")
saveRDS(seurat_merged, "/scratch/hnatri/CART/SoupX/seurat_merged.rds")
#q(save = "no")
#seurat_merged <- readRDS("/scratch/hnatri/CART/SoupX/seurat_merged.rds")
#seurat_soupx_list <- readRDS("/scratch/hnatri/CART/SoupX/seurat_soupx_list.rds")

# Plotting to determine filtering thresholds
par(mfrow=c(2,3))

# PLOTS 1 & 2: nCount vs. nFeature
smoothScatter(log2(seurat_merged$nCount_SoupX_RNA), log2(seurat_merged$nCount_SoupX_RNA),
              xlab = "log2(nCount_SoupX_RNA)", ylab = "log2(nFeature_SoupX_RNA)")

smoothScatter(seurat_merged$nCount_SoupX_RNA, seurat_merged$nCount_SoupX_RNA,
              xlab = "nCount_SoupX_RNA", ylab = "nFeature_SoupX_RNA")

# PLOTS 3 & 4: nCount vs. percent.mt_RNA
smoothScatter(seurat_merged$percent.mt_SoupX_RNA, log2(seurat_merged$nCount_SoupX_RNA),
              xlab = "% MT", ylab = "log2(nCount_SoupX_RNA)")

smoothScatter(seurat_merged$percent.mt_SoupX_RNA, seurat_merged$nCount_SoupX_RNA,
              xlab = "% MT", ylab = "nCount_SoupX_RNA")
abline(v = 10, h = 1000, 
       lty = "dashed", lwd = 1.25, col = "red")

# PLOTS 5 & 6: nFeature vs. percent.mt_RNA
smoothScatter(seurat_merged$percent.mt_SoupX_RNA, log2(seurat_merged$nFeature_SoupX_RNA),
              xlab = "% MT", ylab = "log2(nFeature_RNA)")

smoothScatter(seurat_merged$percent.mt_SoupX_RNA, seurat_merged$nFeature_SoupX_RNA,
              xlab = "% MT", ylab = "nFeature_RNA")
abline(v = 10, h = 500, 
       lty = "dashed", lwd = 1.25, col = "red")

#seurat_filtered <- subset(seurat_merged, subset = percent.mt_SoupX_RNA < 20 &
#                            nFeature_SoupX_RNA > 500 & nCount_SoupX_RNA >1000)

seurat_soupx_list <- filter_manual_mt_rb_genes(sample_seurat_list = seurat_soupx_list,
                                               pt_mt = 20,
                                               nFeature = 500,
                                               nCount = 1000)

saveRDS(seurat_soupx_list, "/scratch/hnatri/CART/TGFb_scRNAseq_list_filtered.rds")
#seurat_soupx_list <- readRDS("/scratch/hnatri/CART/TGFb_scRNAseq_list_filtered.rds")

# Adding cell cycle scores
seurat_soupx_list <- add_cell_cycle_score(seurat_soupx_list)

head(seurat_soupx_list[[1]]@meta.data)

# Integration
seurat_soupx_list <- run_sctransform(seurat_list = seurat_soupx_list,
                                     n_variable_features = 1000,
                                     vars_to_regress = c("percent.mt_RNA",
                                                         "percent.ribo_RNA",
                                                         "S.Score",
                                                         "G2M.Score"))

saveRDS(seurat_soupx_list, "/scratch/hnatri/CART/TGFb_scRNAseq_list_filtered.rds")
#seurat_soupx_list <- readRDS("/scratch/hnatri/CART/TGFb_scRNAseq_list_filtered.rds")

integrated_seurat <- sct_rpca_integration(seurat_list = seurat_soupx_list)

saveRDS(integrated_seurat, "/scratch/hnatri/CART/TGFb_scRNAseq_integrated_seurat.rds")
#integrated_seurat <- readRDS("/scratch/hnatri/CART/TGFb_scRNAseq_integrated_seurat.rds")

#==============================================================================
# Visualization
#==============================================================================

colnames(integrated_seurat@meta.data)
head(integrated_seurat@meta.data)
table(integrated_seurat$seurat_clusters)

# Cluster colors
clusters <- as.factor(c(0, seq(1:22)))
cluster_col <- colorRampPalette(brewer.pal(11, "Paired"))(length(clusters))
names(cluster_col) <- levels(cluster_col)

DimPlot(integrated_seurat,
        group.by = "seurat_clusters",
        cols = cluster_col,
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

DimPlot(integrated_seurat,
        group.by = "seurat_clusters",
        split.by = "KO_WT",
        #cols = cluster_col,
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

table(integrated_seurat$Sample)

# Reads mapping to the construct
# The name of the construct is "IL13OP"
integrated_seurat <- JoinLayers(integrated_seurat, assay = "RNA")
layer_data_RNA <- LayerData(integrated_seurat,
                            assay = "RNA",
                            layer = "counts")
layer_data_RNA <- as.data.frame(t(layer_data_RNA))

range(layer_data_RNA[,"IL13OP"])
range(layer_data_RNA)

p1 <- layer_data_RNA %>% ggplot(aes(x = IL13OP)) +
  geom_histogram() +
  theme_classic()

p2 <- layer_data_RNA %>% ggplot(aes(x = log2(IL13OP))) +
  geom_histogram() +
  theme_classic()

p1 + p2

# What % of the library mapped to the construct?
layer_data_RNA %>%
  dplyr::mutate(sum = rowSums(across(1:(ncol(layer_data_RNA))))) %>%
  dplyr::summarise(CAR_counts = sum(IL13OP),
                   lib_size = sum(sum),
                   pct_CAR = (CAR_counts/lib_size)*100) %>%
  ungroup()

layer_data_sum <- layer_data_RNA %>%
  dplyr::mutate(sum = rowSums(across(1:(ncol(layer_data_RNA)))))

layer_data_sum[1:10, (ncol(layer_data_sum)-10):ncol(layer_data_sum)]

# What % of cells reach a threshold for CAR reads?
nrow(layer_data_RNA)
layer_data_RNA %>% filter(`IL13OP` >= 1) %>% nrow()
layer_data_RNA %>% filter(`IL13OP` >= 2) %>% nrow()
layer_data_RNA %>% filter(`IL13OP` >= 3) %>% nrow()
