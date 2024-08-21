#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 2023/07/13
# Description: infercnv analysis on 13384 CAR T trial tumor samples
#==============================================================================#

#==============================================================================#
# Loading libraries
#==============================================================================#

#.libPaths("/home/hnatri/R/library")
#.libPaths("/home/hnatri/R/4.1_libs")
# install.packages("/home/hnatri/infercnv/", repos=NULL, type="source")

library(Seurat)
library(infercnv) # inferCNV

#==============================================================================#
# Environment variables
#==============================================================================#

set.seed(1234)
options(future.globals.maxSize = 30000 * 1024^2)

#==============================================================================#
# Prepare files for infercnv
#==============================================================================#

seurat_data <- readRDS("/scratch/hnatri/CART/TGFb_scRNAseq_integrated_seurat.rds")
Idents(seurat_data) <- seurat_data$KO_WT

seurat_data <- RenameCells(object = seurat_data,
                           new.names = gsub("-", "_", rownames(seurat_data@meta.data)))

# Need to rename cells to avoid ones starting with numbers
#seurat_data <- RenameCells(object = seurat_data, add.cell.id = "cell")
#seurat_data <- RenameCells(object = seurat_data, new.names = paste0("tumor_", colnames(seurat_data)))
head(colnames(seurat_data))

seurat_data$Sample_seurat_clusters <- paste0(seurat_data$Sample, "_", seurat_data$seurat_clusters)
#Idents(seurat_data) <- seurat_data$Sample_seurat_clusters

#==============================================================================#
# infercnv analysis
#==============================================================================#

# A function for creating the inputs for infercnv for each cluster
make_data <- function(seurat_object){
  message("Creating inputs")

  #selected <- WhichCells(seurat_data, idents = cluster)
  #seurat_subset <- subset(seurat_data, cells = selected)
  
  # Extract out the count matrix
  singleCell.counts.matrix <- LayerData(seurat_object, layer = "counts.1", assay = "RNA")
  
  # Change hyphen to underscore (infercnv will give erorrs if there's hyphen in cellbarcode)
  colnames(singleCell.counts.matrix) <- gsub("-", "_", colnames(singleCell.counts.matrix))
  
  # Create a cellAnnotations file
  cellAnnotations <- as.data.frame(cbind(rownames(seurat_object@meta.data),
                                         seurat_object@meta.data$KO_WT))
  
  colnames(cellAnnotations) <- NULL
  cellAnnotations[,1] <- gsub("-", "_", cellAnnotations[,1])
  rownames(cellAnnotations) <- cellAnnotations[,1]
  cellAnnotations[,1] <- NULL
  
  #setdiff(rownames(cellAnnotations), colnames(singleCell.counts.matrix))
  #setdiff(colnames(singleCell.counts.matrix), rownames(cellAnnotations))
  
  cellAnnotations <- cellAnnotations[colnames(singleCell.counts.matrix),]
  
  write.table(cellAnnotations, file=paste0("/scratch/hnatri/CART/infercnv/seurat_data_cellAnnotations_KOWT.txt"),
              sep= "\t", quote = F)
  
  write.table(round(singleCell.counts.matrix, digits=3), 
              file=paste0("/scratch/hnatri/CART/infercnv/seurat_data_singleCell.counts.matrix_KOWT"), quote=F, sep="\t")
}

# A function for running infercnv for each cluster
run_infercnv <- function(){
  message("Running infercnv")
  # Create the infercnv object
  # Note: need to save all the related files in the working directory
  infercnv_obj <- CreateInfercnvObject(raw_counts_matrix=paste0("/scratch/hnatri/CART/infercnv/seurat_data_singleCell.counts.matrix_KOWT"),
                                       annotations_file=paste0("/scratch/hnatri/CART/infercnv/seurat_data_cellAnnotations_KOWT.txt"),
                                       delim="\t",
                                       gene_order_file="/home/hnatri/13384_CART/13384_Tumors/gene_ordering_file.txt", # this file has geneID, chr, start, end
                                       ref_group_names="WT",
                                       min_max_counts_per_cell=c(100,Inf))
  
  # Perform infercnv operations to reveal CNV signal
  infercnv_obj <- infercnv::run(infercnv_obj,
                                min_cells_per_gene = 10, #default is 3
                                cutoff=0.1,  # expressed genes, use 1 for smart-seq, 0.1 for 10x-genomics
                                out_dir=paste0("/scratch/hnatri/CART/infercnv/seurat_data_HMMi3_samplelevel_output_refgroup/"),  # dir is auto-created for storing outputs
                                cluster_by_groups=T,   # group cells by type
                                denoise=T,
                                num_threads=8,
                                analysis_mode='subclusters', # c("samples", "subclusters", "cells")
                                no_plot=T,
                                HMM=T,
                                HMM_type="i3", # 3 states (neutral, depletion, duplication)
                                hclust_method='ward.D2',
                                sd_amplifier=3,  # sets midpoint for logistic
                                noise_logistic=T, # turns gradient filtering on
                                #num_ref_groups=2, # no normal cells defined
                                tumor_subcluster_partition_method='random_trees')
  #k_obs_groups = 5) # split observation groups into 5 to find malignant + normal cells
  
  # Saving the infercnv result object
  message("Saving the infercnv result object")
  saveRDS(infercnv_obj, file = paste0("/scratch/hnatri/CART/infercnv/seurat_data_HMMi3_samplelevel_refgroup_subclusters.rds"))
  
  # Plotting
  message("Plotting")
  setwd("/scratch/hnatri/CART/infercnv/")
  plot_cnv(infercnv_obj,
           out_dir=".",
           obs_title="Observations (Cells)",
           ref_title="References (Cells)",
           cluster_by_groups=T,
           x.center=1,
           x.range="auto",
           hclust_method='ward.D2',
           # color_safe_pal=TRUE, #using a color blindness safe palette
           output_filename=paste0("seurat_data_HMMi3_samplelevel_refgroup_subclusters_res"),
           output_format="pdf",
           #k_obs_groups = 5,
           dynamic_resize=0)
  
  message("Adding results to the Seurat object")
  
  # Adjust cell barcode in seurat object to match with inferncv (inferncv uses "_", not "-")
  annot <- read.table(paste0("/scratch/hnatri/CART/infercnv/seurat_data_cellAnnotations_KOWT.txt"))
  #selected <- WhichCells(seurat_data, idents = cluster)
  
  #rownames(seurat_data@meta.data) <- gsub("-", "_", rownames(seurat_data@meta.data))
  seurat_subset <- subset(seurat_data, cells = annot$V1)
  colnames(seurat_subset@assays$RNA$counts.1) <- rownames(seurat_subset@meta.data)
  #colnames(seurat_subset@assays$RNA$data) <- rownames(seurat_subset@meta.data)
  
  # Add HMM info into Seurat object 
  seurat_subset <- infercnv::add_to_seurat(infercnv_output_path=paste0("/scratch/hnatri/CART/infercnv/seurat_data_HMMi3_subclusters_refgroup_output/"),
                                           seurat_obj=seurat_subset, 
                                           top_n=10)
  
  # Adjust cell barcodes to match with other assays in the Seurat object
  colnames(seurat_subset@assays$RNA@counts.1) <- rownames(seurat_subset@meta.data)
  #colnames(seurat_subset@assays$RNA@data) <- rownames(seurat_subset@meta.data)
  seurat_subset[["integrated_sct_umap"]] <- RenameCells(
    object = seurat_subset[["integratedSCTumap"]],
    new.names = colnames(seurat_subset@assays$RNA$counts.1))
  
  saveRDS(seurat_subset, file = paste0("/scratch/hnatri/CART/infercnv/seurat_data_infercnv_HMMi3_subclusters_refgroup_seurat.rds"))
  
  message("Finnished running")
  
}

# Running for all cells
make_data(seurat_data)
run_infercnv()

