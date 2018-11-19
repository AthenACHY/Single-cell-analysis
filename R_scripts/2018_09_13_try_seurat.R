lapply(c("dplyr", "ggplot2", "cowplot", "tidyr", "latex2exp", "gridExtra", "ggthemes", "statmod", "Matrix"), require, character.only = TRUE)
library("Seurat")

pbmc.data <- Read10X(data.dir = "/home/a/Documents/Single_cell_tutorial/pbmc3k_filtered_gene_bc_matrices/hg19/")

pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 3, min.genes = 200, 
                           project = "10X_PBMC")
mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)

# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)













d10x <- Read10X_h5("/home/a/Documents/Cancer_CRISPR_lib/Maggie_chow/apop.sc/apopsc/outs/filtered_gene_bc_matrices_h5.h5")
apsoc <- CreateSeuratObject(raw.data = d10x, min.cells = 3, min.genes = 200, 
                           project = "sc_apop")

d10x<-Read10X_h5("/home/a/Documents/Cancer_CRISPR_lib/Maggie_chow/apop.sc/ApopSc180910/apopsImatapop/outs/raw_gene_bc_matrices_h5.h5")
rawapop<-CreateSeuratObject(raw.data = d10x, min.cells = 3, min.genes = 200, 
                            project = "sc_apop")

#### how to aggregate counts? ###

d10x <- Read10X_h5("/home/a/Documents/Cancer_CRISPR_lib/Maggie_chow/apop.sc/ApopSc180910/apopAggr/outs/filtered_gene_bc_matrices_h5.h5")



d10x_seurat <- CreateSeuratObject(
  d10x,
  project = "scapop_agg", 
  min.cells = 3,
  min.genes = 200)


ids <- read.csv("/home/a/Documents/Cancer_CRISPR_lib/Maggie_chow/apop.sc/ApopSc180910/apopAggr/outs/cell_barcodes.csv", stringsAsFactors=F, sep=",", header=F)
ids$samples<-sapply(ids[, 2], function(x) substring(x, nchar(x), nchar(x)))
samples=data.frame(ids$samples)
rownames(samples)<-ids[, 2]
test <- AddMetaData(object = d10x_seurat, metadata = samples)
