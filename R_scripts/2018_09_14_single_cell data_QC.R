#### lookinto Maggies spaperated data ###
### Aggregated data ###
d10x <- Read10X_h5("/Users/athchu/Documents/apop.sc/ApopSc180910/apopsImatapop/outs/raw_gene_bc_matrices_h5.h5")
apsoc <- CreateSeuratObject(raw.data = d10x, min.cells = 3, min.genes = 200, 
                            project = "sc_apop")
### Contlive ###
d10x<-Read10X_h5("/Users/ac/Documents/apop.sc/ApopSc180910/apopscContlive/outs/raw_gene_bc_matrices_h5.h5")
rawCont<-CreateSeuratObject(raw.data = d10x, min.cells = 3, min.genes = 200, 
                            project = "raw_Contlive")

d10x<-Read10X_h5("/Users/ac/Documents/apop.sc/ApopSc180910/apopscContlive/outs/filtered_gene_bc_matrices_h5.h5")
filCont<-CreateSeuratObject(raw.data = d10x, min.cells = 3, min.genes = 200, 
                            project = "filter_Contlive")


### Imatlive ###
d10x<-Read10X_h5("/Users/ac/Documents/apop.sc/ApopSc180910/apopscImatlive/outs/raw_gene_bc_matrices_h5.h5")
rawImat<-CreateSeuratObject(raw.data = d10x, min.cells = 3, min.genes = 200, 
                            project = "raw_imatlive")

d10x<-Read10X_h5("/Users/ac/Documents/apop.sc/ApopSc180910/apopscImatlive/outs/filtered_gene_bc_matrices_h5.h5")
filImat<-CreateSeuratObject(raw.data = d10x, min.cells = 3, min.genes = 200, 
                            project = "filter_imatlive")

### imatapop ###
d10x<-Read10X_h5("/Users/ac/Documents/apop.sc/ApopSc180910/apopsImatapop/outs/raw_gene_bc_matrices_h5.h5")
rawapop<-CreateSeuratObject(raw.data = d10x, min.cells = 3, min.genes = 200, 
                            project = "raw_apop")

d10x<-Read10X_h5("/Users/ac/Documents/apop.sc/ApopSc180910/apopsImatapop/outs/filtered_gene_bc_matrices_h5.h5")
filapop<-CreateSeuratObject(raw.data = d10x, min.cells = 3, min.genes = 200, 
                            project = "imiat_apop")



### compare UMI curve ###
bind_rows(rawCont@meta.data %>% dplyr::mutate(samplename=rownames(rawCont@meta.data), UMIrank=rank(-nUMI), names="raw", project="Contlive"),
          filCont@meta.data %>% dplyr::mutate(samplename=rownames(filCont@meta.data), UMIrank=rank(-nUMI), names="filtered", project="Contlive"),
          rawImat@meta.data %>% dplyr::mutate(samplename=rownames(rawImat@meta.data), UMIrank=rank(-nUMI), names="raw", project="imatlive"),
          filImat@meta.data %>% dplyr::mutate(samplename=rownames(filImat@meta.data), UMIrank=rank(-nUMI), names="filtered", project="imatlive"),
          rawapop@meta.data %>% dplyr::mutate(samplename=rownames(rawapop@meta.data), UMIrank=rank(-nUMI), names="raw", project="imatapop"),
          filapop@meta.data %>% dplyr::mutate(samplename=rownames(filapop@meta.data), UMIrank=rank(-nUMI), names="filtered", project="imatapop")) %>%
  ggplot(aes(x=UMIrank, y=nUMI, colour=names, shape=project)) +
  geom_point()+
  facet_grid(.~project)

bind_rows(rawCont@meta.data %>% dplyr::mutate(samplename=rownames(rawCont@meta.data), UMIrank=rank(-nUMI), names="raw", project="Contlive"),
          filCont@meta.data %>% dplyr::mutate(samplename=rownames(filCont@meta.data), UMIrank=rank(-nUMI), names="filtered", project="Contlive"),
          rawImat@meta.data %>% dplyr::mutate(samplename=rownames(rawImat@meta.data), UMIrank=rank(-nUMI), names="raw", project="imatlive"),
          filImat@meta.data %>% dplyr::mutate(samplename=rownames(filImat@meta.data), UMIrank=rank(-nUMI), names="filtered", project="imatlive"),
          rawapop@meta.data %>% dplyr::mutate(samplename=rownames(rawapop@meta.data), UMIrank=rank(-nUMI), names="raw", project="imatapop"),
          filapop@meta.data %>% dplyr::mutate(samplename=rownames(filapop@meta.data), UMIrank=rank(-nUMI), names="filtered", project="imatapop")) %>%
  ggplot(aes(x=interaction(project, names), y=nUMI, colour=names, shape=project)) +
  geom_boxplot(position = "dodge")
#### ok so, Imatapop generally have cells with a lot of UMi counts, thus the truncation of the cell tail is dramatic? ### 
#### for the filtered Apop, min nGene = 1173 and min nUMI = 3968 ###
#### for the filtered Imat, min nGene = 838 and min nUMI = 1743 ###
#### for the filtered Cont, min nGene = 618 and min nUMI = 1292 ###

median(rawapop@meta.data$nUMI)
#[1] 866.5

median(rawCont@meta.data$nUMI)
#[1] 4948

median(rawImat@meta.data$nUMI)
#[1] 7131
### use UMI threshold for subsetting all samples ###
rawapop_filtered<- FilterCells(
  object = rawapop, 
  subset.names = c("nUMI"), 
  high.thresholds = c(30000),
  low.thresholds=c(173))

rawCont_filtered<- FilterCells(
  object = rawCont, 
  subset.names = c("nUMI"), 
  high.thresholds = c(30000),
  low.thresholds=c(989))

rawImat_filtered<- FilterCells(
  object = rawImat, 
  subset.names = c("nUMI"), 
  high.thresholds = c(30000),
  low.thresholds=c(1426))

####Merge the Seurat objects ###

Mag.combined <- MergeSeurat(object1 = rawCont_filtered, 
                             object2 = rawImat_filtered, 
                             add.cell.id1 = "Cont", 
                             add.cell.id2 = "Imat", project = "10xAggr",
                             do.normalize=F)
table(Mag.combined@meta.data$orig.ident)

Mag.combined <- AddSamples(object = Mag.combined, new.data=rawapop_filtered@data,
                           add.cell.id = "Apop", do.normalize =F)
### orig.ident were messed up what to do? ###
Mag.combined@meta.data$orig.ident<-Mag.combined@ident

### Convert to singlecell object (Seurat Convert did not work....)###
Mag_sce <- SingleCellExperiment(assay=list(counts=as.matrix(Mag.combined@raw.data)), colData=Mag.combined@meta.data)
### find genes that express in some cells ###
keep_feature <-rowSums(counts(Mag_sce) > 0) > 1000
Mag_sce_keep <- umi[keep_feature, ]
### try to do QC matrix by count itself###
Mag_sce_QC<-calculateQCMetrics(Mag_sce)
#Note that the names of some metrics have changed, see 'Renamed metrics' in ?calculateQCMetrics.
#Old names are currently maintained for back-compatibility, but may be removed in future releases.
plotHighestExprs(Mag_sce_QC)
### R crash...#
plotPCA(object = Mag_sce_QC, colour_by = "ident")
### subset a list of cells for PCA? ###


### PCA ###
pbmc <- NormalizeData(object = Mag.combined, normalization.method = "LogNormalize", 
                      scale.factor = 10000)

pbmc<- FindVariableGenes(object =  Mag.combined, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
pbmc <- ScaleData(object = Mag.combined, vars.to.regress = c("nUMI"))

pbmc <- RunPCA(object = Mag.combined, pc.genes =  Mag.combined@var.genes, do.print = TRUE, pcs.print = 1:5, 
               genes.print = 5)

PCAPlot(object = pbmc, dim.1 = 1, dim.2 = 2)
#> rawapop
#An object of class seurat in project sc_apop 
#13482 genes across 1518 samples.
# So 10x pipeline cut out 1400 cells? #
umi_plot<-rawapop@meta.data %>% dplyr::arrange(desc(nUMI)) %>%
  dplyr::mutate(samplename=rownames(rawapop@meta.data), UMIrank=rank(-nUMI)) %>%
  ggplot(aes(x=UMIrank, y=nUMI)) +
  geom_point()

