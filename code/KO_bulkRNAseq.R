#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 07/23/2023
# Description: TGFbRa2KO, CAR/TGFb bulk-RNAseq data analysis: DEGs, GSEA
#==============================================================================#

#==============================================================================#
# Loading libraries
#==============================================================================#

library(edgeR)
library(DESeq2)
library(limma)
library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyr)
library(tidyverse)
library(googlesheets4)
library(readxl)
library(biomaRt)
library(CaSpER)
library(org.Hs.eg.db)
library(clusterProfiler)

#==============================================================================#
# Environment variables
#==============================================================================#

set.seed(1234)

#==============================================================================#
# Import data
#==============================================================================#

# https://docs.google.com/spreadsheets/d/1T3oB6dS_rbrijBWKAlPUCOWNAFc2Zz1yHDsUJRSdpOg/edit?usp=sharing
gs4_deauth()
TGFbRa2_KO_metadata  <- gs4_get("https://docs.google.com/spreadsheets/d/1T3oB6dS_rbrijBWKAlPUCOWNAFc2Zz1yHDsUJRSdpOg/edit?usp=sharing")
sheet_names(TGFbRa2_KO_metadata)
metadata_21446 <- read_sheet(TGFbRa2_KO_metadata, sheet = "21446")
metadata_21547 <- read_sheet(TGFbRa2_KO_metadata, sheet = "21547")

counts_21446 <- read_excel("/tgen_labs/banovich/BCTCSF/TGFbR2_KO/21446/Analysis/counts_cpm.xlsx",
                         sheet = "raw counts")
counts_21547 <- read_excel("/tgen_labs/banovich/BCTCSF/TGFbR2_KO/21547/Analysis/counts_cpm_fpkm.xlsx",
                        sheet = "raw counts")

#colnames(counts_21547) <- c(c("drop", "Geneid"), colnames(cpm_21547)[-1:-2])

# Gene lengths
gene_info <- counts_21446[1:2]
colnames(gene_info) <- c("hgnc_symbol", "length")

# Adding gene positions
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
biomart_info <- getBM(attributes=c('chromosome_name',
                                   'start_position',
                                   'end_position',
                                   'strand',
                                   'ensembl_gene_id',
                                   'hgnc_symbol'),
                      filters=c('hgnc_symbol'),
                      values=list(gene_info$hgnc_symbol),
                      mart=ensembl)
gene_info <- merge(gene_info, biomart_info, by = "hgnc_symbol")

counts_21446 <- counts_21446  %>%
  dplyr::select(!Length) %>%
  rename_with(~ paste0(.x, "-21446")) %>%
  column_to_rownames(var = "Geneid-21446")

counts_21547[c(45352, 45353),]

counts_21547 <- counts_21547  %>%
  dplyr::select(!Length) %>%
  rename_with(~ paste0(.x, "-21547")) %>%
  distinct(`Geneid-21547`, .keep_all = TRUE) %>%
  column_to_rownames(var = "Geneid-21547")

counts_all_samples <- merge(counts_21446, counts_21547, by = 0) %>%
  column_to_rownames(var = "Row.names")

# Metadata for differential expression
metadata <- as.data.frame(matrix(nrow = length(counts_all_samples),
                                 ncol = 6))
colnames(metadata) <- c("Sample_Name", "CellLine", "KO_WT", "CAR", "TGFb", "Batch")
metadata$Sample_Name <- c(colnames(counts_21446), colnames(counts_21547))
metadata$CellLine <- substr(metadata$Sample_Name, start = 1, stop = 5)
metadata$KO_WT <- substr(metadata$Sample_Name, start = 6, stop = 7)
metadata$CAR <- sapply(strsplit(metadata$Sample_Name, split='-', fixed=TRUE), `[`, 2)
metadata$TGFb <- sapply(strsplit(metadata$Sample_Name, split='-', fixed=TRUE), `[`, 3)
metadata$Batch <- sapply(strsplit(metadata$Sample_Name, split='-', fixed=TRUE), `[`, 4)

# Creating the DGEList object
gene_info <- gene_info %>%
  dplyr::filter(hgnc_symbol %in% rownames(counts_all_samples))
gene_info <- gene_info[match(rownames(counts_all_samples), gene_info$hgnc_symbol),]

dge <- DGEList(counts = counts_all_samples, genes = gene_info, group = metadata$KO_WT)
colnames(dge) <- colnames(counts_all_samples)
dge$samples$CellLine <- metadata$CellLine
dge$samples$KO_WT <- metadata$KO_WT
dge$samples$TGFb <- metadata$TGFb
dge$samples$Batch <- metadata$Batch

#==============================================================================#
# Dimensionality reduction and differential expression for KO and WT
#==============================================================================#

# Filtering
keep <- filterByExpr(dge, group = "group")
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Creating a design model matri
KO_WT <- dge$samples$KO_WT
KO_WT <- factor(KO_WT, levels = c("KO", "WT"))
design <- model.matrix(~0+KO_WT)
colnames(design) <- levels(KO_WT)

design

# Running limma voom. Normalizes expression intensities so that the log-ratios
# have similar distributions across a set of arrays.
v <- voomWithQualityWeights(dge, design, plot=TRUE, normalize.method="quantile")

#===============================
# Clustering analysis: multi-dimensional scaling plot
#===============================

# If gene.selection is "common", then the top genes are those with the largest standard deviations between samples.
# If gene.selection is "pairwise", then a different set of top genes is selected for each pair of samples.
plotMDS(v, top = 100, plot=TRUE, pch = ifelse(v$targets$KO_WT=="KO",19,15),
        col = ifelse(v$targets$TGFb=="TGFbm","deepskyblue4","firebrick")) #, gene.selection = "pairwise", labels=tissue,
legend("bottomright", pch = c(19,19,15,15),
       col = c("deepskyblue4","firebrick","deepskyblue4","firebrick"),
       legend = c("KO TGFbm", "KO TGFbp", "WT TGFbm", "WT TGFbp"),
       ncol = 2)

plotMDS(v, top = 100, plot=TRUE, pch = ifelse(v$targets$KO_WT=="KO",19,15),
        col = ifelse(v$targets$CellLine=="HD562","springgreen3",
                     ifelse(v$targets$CellLine=="HD624","yellow3","cyan4"))) #, gene.selection = "pairwise", labels=tissue,
legend("bottomright", pch = c(19,19,19,15,15,15),
       col = c("springgreen3","yellow3","cyan4", "springgreen3","yellow3","cyan4"),
       legend = c("KO HD562", "KO HD624", "KO HD680", "WT HD562", "WT HD624", "WT HD680"),
       ncol = 2)

#===============================
# Differential expression
#===============================

# Block design for cell line This is used in KO-WT comparisons with paired samples.
corfit <- duplicateCorrelation(v, design, block = dge$samples$CellLine)
# This should give a positive correlation value. It represents the 
# correlation between measurements made on the same person.
corfit$consensus

# Fitting the linear model
# If using paired damples, add the correlation and block design for sibship
fit <- lmFit(v, design, block = dge$samples$CellLine, correlation = corfit$consensus)

# Contrast design for differential expression
contrast_design <- makeContrasts(KO-WT, levels = design)
contrast_design

# Running contrast analysis
vfit <- contrasts.fit(fit, contrasts = contrast_design)
summary(decideTests(vfit))

# Getting genes with logFC >2 and adjusted P-value <0.01. Robust = should the estimation 
# of the empirical Bayes prior parameters be robustified against outlier sample variances?
# Log2fc of 1 is equivalent to linear fold change of 2.
vebayesfit <- eBayes(vfit, robust=TRUE)
plotSA(vebayesfit, main="Final model: Mean−variance trend")
vtoptable <- topTable(vebayesfit, coef=1, n=Inf, p.value=0.01, lfc=2)
vtoptable

# Results for all genes
vallresults <- topTable(vebayesfit, coef=1, n=Inf)

# Number of selected significant genes
nrow(vtoptable)
table(vtoptable$chromosome_name)

# # DEGs by chromosome
vtoptable <- vtoptable %>% filter(chromosome_name %in% c(seq(1, 22), "X"))
prop_table = as.data.frame(table(vtoptable$chromosome_name))
colnames(prop_table) = c("chr", "DEGs")
prop_table = spread(prop_table, chr, DEGs)
# Converting to percentange
prop_table[,1:length(prop_table)] = (prop_table[,1:length(prop_table)]/rowSums(prop_table[,1:length(prop_table)]))*100
prop_table = gather(prop_table, chr, DEGs, names(prop_table)[1:length(names(prop_table))], factor_key=TRUE)

prop_table$chr <- factor(prop_table$chr, levels = c(seq(1, 22), "X"))

ggplot(prop_table, aes(x=chr, y = DEGs, fill = chr)) +
  geom_bar(stat="identity", position='stack', width = 0.8) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  NoLegend()

# Chr 3 genes
vtoptable %>% filter(chromosome_name == 3)

#==============================================================================#
# Differential expression by treatment
#==============================================================================#

# Creating a DGEList with TGFb treatment as group
dge_TGFb <- DGEList(counts = counts_all_samples, genes = gene_info, group = metadata$TGFb)
colnames(dge_TGFb) <- colnames(counts_all_samples)
dge_TGFb$samples$CellLine <- metadata$CellLine
dge_TGFb$samples$KO_WT <- metadata$KO_WT
dge_TGFb$samples$TGFb <- metadata$TGFb
dge_TGFb$samples$Batch <- metadata$Batch

# Creating a design model matrix
# TGFb- means that it was absent in the coculture media, and TGFb+ means
# TGFb was added in the coculture media
TGFb <- dge_TGFb$samples$TGFb
TGFb <- factor(TGFb, levels = c("TGFbm", "TGFbp"))
design_TGFb <- model.matrix(~0+TGFb)
colnames(design_TGFb) <- levels(TGFb)

design_TGFb

# Running limma voom
v_TGFb <- voomWithQualityWeights(dge_TGFb, design, plot=TRUE, normalize.method="quantile")

# Block design
corfit <- duplicateCorrelation(v_TGFb, design, block = dge_TGFb$samples$CellLine)
corfit$consensus

# Fitting the linear model
fit_TGFb <- lmFit(v_TGFb, design_TGFb, block = dge_TGFb$samples$CellLine,
                  correlation = corfit$consensus)

# Contrast design
contrast_design_TGFb <- makeContrasts(TGFbm-TGFbp, levels = design_TGFb)
contrast_design_TGFb

# Running contrast analysis
vfit_TGFb <- contrasts.fit(fit_TGFb, contrasts = contrast_design_TGFb)
summary(decideTests(vfit_TGFb))

# Getting genes with logFC >2 and adjusted P-value <0.01
vebayesfit_TGF <- eBayes(vfit_TGFb, robust=TRUE)
plotSA(vebayesfit_TGF, main="Final model: Mean−variance trend")
vtoptable_TGF <- topTable(vebayesfit_TGF, coef=1, n=Inf, p.value=0.01, lfc=2)
dim(vtoptable_TGF)

# Results for all genes
vallresults_TGF <- topTable(vebayesfit_TGF, coef=1, n=Inf)

# Number of selected significant genes
nrow(vtoptable_TGF)
table(vtoptable$chromosome_name_TGF)

#==============================================================================#
# Differential expression by TGFb treatment, stratified by KO/WT
#==============================================================================#

dge_TGFb_KO <- dge_TGFb[,which(dge_TGFb$samples$KO_WT == "KO")]
dge_TGFb_WT <- dge_TGFb[,which(dge_TGFb$samples$KO_WT == "WT")]

get_degs <- function(dgelist){
  # Creating a design model matrix
  TGFb <- dgelist$samples$TGFb
  TGFb <- factor(TGFb, levels = c("TGFbm", "TGFbp"))
  design_TGFb <- model.matrix(~0+TGFb)
  colnames(design_TGFb) <- levels(TGFb)
  
  design_TGFb
  
  # Running limma voom
  v_TGFb <- voomWithQualityWeights(dgelist, design_TGFb, plot=TRUE, normalize.method="quantile")
  
  # Block design
  corfit <- duplicateCorrelation(v_TGFb, design_TGFb, block = dgelist$samples$CellLine)
  
  # Fitting the linear model
  fit_TGFb <- lmFit(v_TGFb, design_TGFb, block = dgelist$samples$CellLine,
                    correlation = corfit$consensus)
  
  # Contrast design
  contrast_design_TGFb <- makeContrasts(TGFbm-TGFbp, levels = design_TGFb)
  
  # Running contrast analysis
  vfit_TGFb <- contrasts.fit(fit_TGFb, contrasts = contrast_design_TGFb)
  vebayesfit_TGF <- eBayes(vfit_TGFb, robust=TRUE)
  
  # Results for all genes
  vallresults_TGF <- topTable(vebayesfit_TGF, coef=1, n=Inf)
  
  vallresults_TGF
}

degs_KO <- get_degs(dge_TGFb_KO)
degs_WT <- get_degs(dge_TGFb_WT)

# Significant DEGs
degs_KO_sig <- degs_KO %>% filter(abs(logFC)>2,
                                  adj.P.Val < 0.01)
degs_WT_sig <- degs_WT %>% filter(abs(logFC)>2,
                                  adj.P.Val < 0.01)

length(intersect(rownames(degs_KO_sig), rownames(degs_WT_sig)))
length(setdiff(rownames(degs_KO_sig), rownames(degs_WT_sig)))
length(setdiff(rownames(degs_WT_sig), rownames(degs_KO_sig)))

#==============================================================================#
# Differential expression by KO/WT, stratified by TGFb treatment
#==============================================================================#

dge_TGFbm <- dge[,which(dge$samples$TGFb == "TGFbm")]
dge_TGFbp <- dge[,which(dge$samples$TGFb == "TGFbp")]

get_degs <- function(dgelist){
  # Creating a design model matrix
  KOWT <- dgelist$samples$KO_WT
  KOWT <- factor(KOWT, levels = c("KO", "WT"))
  design_KOWT <- model.matrix(~0+KOWT)
  colnames(design_KOWT) <- levels(KOWT)
  
  # Running limma voom
  v_KOWT <- voomWithQualityWeights(dgelist, design_KOWT, plot=TRUE, normalize.method="quantile")
  
  # Block design
  corfit <- duplicateCorrelation(v_KOWT, design_KOWT, block = dgelist$samples$CellLine)
  
  # Fitting the linear model
  fit_KOWT <- lmFit(v_KOWT, design_KOWT, block = dgelist$samples$CellLine,
                    correlation = corfit$consensus)
  
  # Contrast design
  contrast_design_KOWT <- makeContrasts(KO-WT, levels = design_KOWT)
  
  # Running contrast analysis
  vfit_KOWT <- contrasts.fit(fit_KOWT, contrasts = contrast_design_KOWT)
  vebayesfit_KOWT <- eBayes(vfit_KOWT, robust=TRUE)
  
  # Results for all genes
  vallresults_KOWT <- topTable(vebayesfit_KOWT, coef=1, n=Inf)
  
  vallresults_KOWT
}

degs_TGFbm <- get_degs(dge_TGFbm)
degs_TGFbp <- get_degs(dge_TGFbp)

# Significant DEGs
degs_TGFbm_sig <- degs_TGFbm %>% filter(abs(logFC)>2,
                                  adj.P.Val < 0.01)
degs_TGFbp_sig <- degs_TGFbp %>% filter(abs(logFC)>2,
                                  adj.P.Val < 0.01)

length(intersect(rownames(degs_TGFbm_sig), rownames(degs_TGFbp_sig)))
length(setdiff(rownames(degs_TGFbm_sig), rownames(degs_TGFbp_sig)))
length(setdiff(rownames(degs_TGFbp_sig), rownames(degs_TGFbm_sig)))

#==============================================================================#
# Pathways analysis with clusterProfiler
#==============================================================================#

# Nontreated
TGFbm_gene_list <- degs_TGFbm_sig$logFC
names(TGFbm_gene_list) <- degs_TGFbm_sig$ensembl_gene_id

# omit any NA values 
TGFbm_gene_list <- na.omit(TGFbm_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
TGFbm_gene_list <- sort(TGFbm_gene_list, decreasing = TRUE)

TGFbm_gse <- gseGO(geneList = TGFbm_gene_list, 
             ont = "ALL", 
             keyType = "ENSEMBL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db", 
             pAdjustMethod = "none")

gsp1 <- dotplot(TGFbm_gse, showCategory=10, split=".sign") + facet_grid(.~.sign)

ridgeplot(TGFbm_gse) + labs(x = "enrichment distribution")

# TGFb treated
TGFbp_gene_list <- degs_TGFbp_sig$logFC
names(TGFbp_gene_list) <- degs_TGFbp_sig$ensembl_gene_id

# omit any NA values 
TGFbp_gene_list <- na.omit(TGFbp_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
TGFbp_gene_list <- sort(TGFbp_gene_list, decreasing = TRUE)

TGFbp_gse <- gseGO(geneList = TGFbp_gene_list, 
                   ont = "ALL", 
                   keyType = "ENSEMBL", 
                   nPerm = 10000, 
                   minGSSize = 3, 
                   maxGSSize = 800, 
                   pvalueCutoff = 0.05, 
                   verbose = TRUE, 
                   OrgDb = "org.Hs.eg.db", 
                   pAdjustMethod = "none")

gsp2 <- dotplot(TGFbp_gse, showCategory=10, split=".sign") + facet_grid(.~.sign)

ridgeplot(TGFbp_gse) + labs(x = "enrichment distribution")

gsp1 + gsp2

#==============================================================================#
# Over-representation analysis of shared and unique DEGs
#==============================================================================#

length(intersect(rownames(degs_TGFbm_sig), rownames(degs_TGFbp_sig)))
length(setdiff(rownames(degs_TGFbm_sig), rownames(degs_TGFbp_sig)))
length(setdiff(rownames(degs_TGFbp_sig), rownames(degs_TGFbm_sig)))

# Gene groups
TGFbp_unique <- setdiff(degs_TGFbp_sig$ensembl_gene_id, degs_TGFbm_sig$ensembl_gene_id)
TGFbm_unique <- setdiff(degs_TGFbm_sig$ensembl_gene_id, degs_TGFbp_sig$ensembl_gene_id)
shared <- intersect(degs_TGFbm_sig$ensembl_gene_id, degs_TGFbp_sig$ensembl_gene_id)

go_enrich_TGFbp <- enrichGO(gene = na.omit(TGFbp_unique),
                            universe = degs_TGFbm$ensembl_gene_id,
                            OrgDb = "org.Hs.eg.db", 
                            keyType = 'ENSEMBL',
                            readable = T,
                            ont = "BP",
                            pvalueCutoff = 0.05, 
                            qvalueCutoff = 0.10)

go_enrich_TGFbm <- enrichGO(gene = na.omit(TGFbm_unique),
                            universe = degs_TGFbm$ensembl_gene_id,
                            OrgDb = "org.Hs.eg.db", 
                            keyType = 'ENSEMBL',
                            readable = T,
                            ont = "BP",
                            pvalueCutoff = 0.05, 
                            qvalueCutoff = 0.10)

go_enrich_shared <- enrichGO(gene = na.omit(shared),
                             universe = degs_TGFbm$ensembl_gene_id,
                             OrgDb = "org.Hs.eg.db", 
                             keyType = 'ENSEMBL',
                             readable = T,
                             ont = "BP",
                             pvalueCutoff = 0.05, 
                             qvalueCutoff = 0.10)

bp1 <- barplot(go_enrich_TGFbp, 
        drop = TRUE, 
        showCategory = 10, 
        title = "DE between KO and WT in TGFb+",
        font.size = 8)

bp2 <- barplot(go_enrich_TGFbm, 
        drop = TRUE, 
        showCategory = 10, 
        title = "DE between KO and WT in TGFb-",
        font.size = 8)

bp3 <- barplot(go_enrich_shared, 
        drop = TRUE, 
        showCategory = 10, 
        title = "DE between KO and WT in TGFb+ and TGFb-",
        font.size = 8)

bp1 + bp2 + bp3

#==============================================================================#
# CNVs
#==============================================================================#

# Cytoband information can be downloaded from UCSC.
# http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz
genome_info <- fread("/scratch/hnatri/cytoBand.txt.gz", header = F)
genome_info <- genome_info %>% 
  dplyr::filter(V1 %in% paste0("chr", c(seq(1, 22), "X"))) %>%
  as.data.frame()

cytoband <- data.frame(V1=gsub("chr", "", genome_info[,1]),
                       V2=genome_info[,2],
                       V3=genome_info[,3],
                       V4=substring(genome_info$V4, 1, 1),
                       stringsAsFactors = F)
start <- do.call(rbind, lapply(split(cytoband$V2, paste0(cytoband$V1, cytoband$V4)), min))
end <- do.call(rbind, lapply(split(cytoband$V3, paste0(cytoband$V1, cytoband$V4)), max))
cytoband <- data.frame(V1=gsub("p", "", gsub("q", "", rownames(start))), V2=start, V3=end, V4=rownames(start), stringsAsFactors=F)
cytoband <- cytoband[as.vector(unlist(sapply(c(1:22, "X"), function(x) which(cytoband$V1 %in% x)))), ]
cytoband$V4[grep("q", cytoband$V4)] <- "q"
cytoband$V4[grep("p", cytoband$V4)] <- "p"
rownames(cytoband) <- NULL

# Centromere information can be downloaded from UCSC.
#curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/centromeres.txt.gz" | gunzip -c | grep acen | head
centromere <- fread("/scratch/hnatri/centromeres.txt.gz", header = F)
centromere <- centromere %>% 
  dplyr::filter(V2 %in% paste0("chr", c(seq(1, 22), "X"))) %>%
  as.data.frame()

annotation <- generateAnnotation(id_type = "hgnc_symbol",
                                 genes = rownames(v$genes),
                                 id_type = "hgnc_symbol",
                                 ishg19 = F,
                                 centromere)
  
