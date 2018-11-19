### try Seurat PCA with different group of cells ###
### Contlive ###
d10x<-Read10X_h5("/Users/ac/Documents/apop.sc/ApopSc180910/apopscContlive/outs/raw_gene_bc_matrices_h5.h5")
rawCont<-CreateSeuratObject(raw.data = d10x, min.cells = 3, min.genes = 200, 
                            project = "raw_Contlive")
### Imatlive ###
d10x<-Read10X_h5("/Users/ac/Documents/apop.sc/ApopSc180910/apopscImatlive/outs/raw_gene_bc_matrices_h5.h5")
rawImat<-CreateSeuratObject(raw.data = d10x, min.cells = 3, min.genes = 200, 
                            project = "raw_imatlive")
### imatapop ###
d10x<-Read10X_h5("/Users/ac/Documents/apop.sc/ApopSc180910/apopsImatapop/outs/raw_gene_bc_matrices_h5.h5")
rawapop<-CreateSeuratObject(raw.data = d10x, min.cells = 3, min.genes = 200, 
                            project = "raw_apop")

### basic stats ###

bind_rows(rawapop@meta.data, rawCont@meta.data, rawImat@meta.data) %>%
  ggplot(aes(x=as.numeric(nUMI), y=as.numeric(nGene)))+
  geom_bin2d(binwidth = c(1000, 1000)) +
  facet_grid(orig.ident ~., scales="free")+
  scale_fill_gradient(low="blue", high="red")+
  geom_vline(xintercept=30000)+
  ylab("nGene")+
  xlab("nUMI")+
  theme_classic()



### get cells have nUMI >10000 ###
rawapop_filtered<- FilterCells(object = rawapop, subset.names = c("nUMI"), 
low.thresholds=c(10000))
rawCont_filtered<- FilterCells(object = rawCont, subset.names = c("nUMI"), 
low.thresholds=c(10000))
rawImat_filtered<- FilterCells(object = rawImat, subset.names = c("nUMI"), 
low.thresholds=c(10000))
### merge Seurat object ###

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
### perform log2notmalisation to total lib size ###
Mag.combined<-NormalizeData(object = Mag.combined, normalization.method = "LogNormalize", 
              scale.factor = 10000)
Mag.combined<- FindVariableGenes(object =  Mag.combined, mean.function = ExpMean, dispersion.function = LogVMR, 
                         x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
### center gene expression ###
Mag.combined <- ScaleData(object = Mag.combined, vars.to.regress = c("nUMI"))

Mag.combined <- RunPCA(object = Mag.combined, pc.genes =  Mag.combined@var.genes, do.print = TRUE, pcs.print = 1:5, 
               genes.print = 5)

w<-PCAPlot(object = Mag.combined, dim.1 = 1, dim.2 = 2)
ggsave(w, file="/Users/athchu/Documents/apop.sc/Apop_10000uMI_pca.png", width=6)
w<-PCElbowPlot(object = Mag.combined) + scale_y_continuous(breaks=seq(1, 9, 1))
ggsave(w, file="/Users/athchu/Documents/apop.sc/Apop_10000uMI_pcasd.png", width=6)

VizPCA(object =  Mag.combined, pcs.use = 1:2)
### jackstraw ###
Mag.combined <- JackStraw(object = Mag.combined, num.replicate = 100, display.progress = FALSE)
JackStrawPlot(object = Mag.combined, PCs = 1:12)
### jackstraw showed that PC 1-12 are still significant, use cut -off as PC10 as it is 1st PC that is less significant ###
### perofrom t-sne ###
Mag.combined <- RunTSNE(object = Mag.combined, dims.use = 1:10, do.fast = TRUE)
w<-TSNEPlot(object = Mag.combined)
ggsave(w, file="/Users/athchu/Documents/apop.sc/Apop_10000uMI_tsne.png", width=6)

### find key clustering markers compared to all remaining cells, report ###
apop.markers <- FindAllMarkers(object =  Mag.combined, only.pos = TRUE, min.pct = 0.25, 
                               thresh.use = 0.25)
apop.markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
w<-FeaturePlot(object = Mag.combined, features.plot = c("NME1", "SRM", "NHP2", "PPP1R14B", "MYC"), cols.use = c("grey", "blue"), 
            reduction.use = "pca")
ggsave(w, file="/Users/athchu/Documents/apop.sc/Apop_10000uMI_Contlive_feat.png", width=6)


### apply other threshold? ###
### get cells have nUMI >200 & <3000 ###
rawapop_filtered<- FilterCells(object = rawapop, subset.names = c("nUMI"), 
                               high.thresholds=c(3000), low.thresholds=c(200))
rawCont_filtered<- FilterCells(object = rawCont, subset.names = c("nUMI"), 
                               high.thresholds=c(3000), low.thresholds=c(200))
rawImat_filtered<- FilterCells(object = rawImat, subset.names = c("nUMI"), 
                               high.thresholds=c(3000), low.thresholds=c(200))
table(bind_rows(rawapop@meta.data, rawCont@meta.data, rawImat@meta.data)$orig.ident)
### merge Seurat object ###

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
table(Mag.combined@meta.data$orig.ident)
### perform log2notmalisation to total lib size ###
Mag.combined<-NormalizeData(object = Mag.combined, normalization.method = "LogNormalize", 
                            scale.factor = 10000)

### x: expression level?, y: dispersion level? ###
Mag.combined<- FindVariableGenes(object =  Mag.combined, mean.function = ExpMean, dispersion.function = LogVMR, 
                                 x.low.cutoff = 1, x.high.cutoff = 3, y.cutoff = 0.5)
length(Mag.combined@var.genes)
### center gene expression ###
Mag.combined <- ScaleData(object = Mag.combined, vars.to.regress = c("nUMI"))
length(Mag.combined@var.genes)

Mag.combined <- RunPCA(object = Mag.combined, pc.genes =  Mag.combined@var.genes, do.print = TRUE, pcs.print = 1:5, 
                       genes.print = 5)

w<-PCAPlot(object = Mag.combined, dim.1 = 1, dim.2 = 2)
ggsave(w, file="/Users/athchu/Documents/apop.sc/Apop_10000uMI_pca.png", width=6)
w<-PCElbowPlot(object = Mag.combined) + scale_y_continuous(breaks=seq(1, 9, 1))
ggsave(w, file="/Users/athchu/Documents/apop.sc/Apop_10000uMI_pcasd.png", width=6)
VizPCA(object =  Mag.combined, pcs.use = 1:2)
### try turkey HSD test using PC1 top5 genes ###
gfi<-c("SLC25A37", "HBA1", "MT-ND4L", "HBA2", "AKAP9", "VIM", "PRDX1", "TUBA1B", "LDHA", "ACTG1" )
f_df<-data.frame(group=Mag.combined@meta.data$orig.ident, exp=Mag.combined@scale.data[gfi[1],])
table(f_df[f_df$exp>0,]$group)
exp.lm <- lm(exp ~ group, data = f_df)
exp.av <- aov(exp.lm)
summary(exp.av)
TukeyHSD(exp.av)


### jackstraw ###
Mag.combined <- JackStraw(object = Mag.combined, num.replicate = 100, display.progress = FALSE)
JackStrawPlot(object = Mag.combined, PCs = 1:10)
### prop.freq for variable genes too low ###


### perofrom t-sne ###
Mag.combined <- RunTSNE(object = Mag.combined, dims.use = 1:10, do.fast = TRUE)
w<-TSNEPlot(object = Mag.combined)
ggsave(w, file="/Users/athchu/Documents/apop.sc/Apop_10000uMI_tsne.png", width=6)

### find key clustering markers compared to all remaining cells, report ###
apop.markers <- FindAllMarkers(object =  Mag.combined, only.pos = TRUE, min.pct = 0.25, 
                               thresh.use = 0.25)
apop.markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
w<-FeaturePlot(object = Mag.combined, features.plot = c("NME1", "SRM", "NHP2", "PPP1R14B", "MYC"), cols.use = c("grey", "blue"), 
               reduction.use = "pca")
ggsave(w, file="/Users/athchu/Documents/apop.sc/Apop_10000uMI_Contlive_feat.png", width=6)


### try the better separation from cell-ranger ###
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
table(Mag.combined@meta.data$orig.ident)
### perform log2notmalisation to total lib size ###
Mag.combined<-NormalizeData(object = Mag.combined, normalization.method = "LogNormalize", 
                            scale.factor = 10000)

### x: expression level?, y: dispersion level? ###
Mag.combined<- FindVariableGenes(object =  Mag.combined, mean.function = ExpMean, dispersion.function = LogVMR, 
                                 x.low.cutoff = 0.5, x.high.cutoff = 3, y.cutoff = 1.5)
length(Mag.combined@var.genes)
### center gene expression based on nUMI and sample origin###
Mag.combined <- ScaleData(object = Mag.combined, vars.to.regress = c("nUMI"))
### Mag.combined <- ScaleData(object = Mag.combined, vars.to.regress = c("nUMI"))
length(Mag.combined@var.genes)
Mag.combined <- RunPCA(object = Mag.combined, pc.genes =  Mag.combined@var.genes, do.print = TRUE, pcs.print = 1:5, 
                       genes.print = 5)
PCAPlot(object = Mag.combined, dim.1 = 1, dim.2 = 2)
PCElbowPlot(object = Mag.combined) + scale_y_continuous(breaks=seq(1, 9, 1))
VizPCA(object =  Mag.combined, pcs.use = 1:2)
Mag.combined <- JackStraw(object = Mag.combined, num.replicate = 100, display.progress = FALSE)
JackStrawPlot(object = Mag.combined, PCs = 1:12)
Mag.combined <- RunTSNE(object = Mag.combined, dims.use = 1:10, do.fast = TRUE)
TSNEPlot(object = Mag.combined)
Mag.combined <- FindClusters(object = Mag.combined, reduction.type = "pca", dims.use = 1:10, 
                     resolution = 0.6, print.output = 0, save.SNN = TRUE)

TSNEPlot(object = Mag.combined, color.use=Mag.combined@meta.data$orig.ident)

### find key markers ###
Mag.markers <- FindAllMarkers(object =  Mag.combined, only.pos = TRUE, min.diff.pct = 0.1, 
                              logfc.threshold = 0.25 ,test.use = "roc")
Mag.markers %>% dplyr::filter(avg_logFC>1)
### find the top genes for each cluster, power >0.5 ###
Mag.markers %>% dplyr::filter(power>0.5)
Marker.apop<-FindMarkers(object = Mag.combined, ident.1 = c(4, 5), ident.2 = c(0, 1, 2, 3), 
                                min.pct = 0.5, only.pos = TRUE)

FindMarkers(object = Mag.combined, ident.1 = c(2), ident.2 = c(1), 
            min.pct = 0.5, only.pos = TRUE)
Marker.imat<-FindMarkers(object = Mag.combined, ident.1 = c(1, 3, 4), ident.2 = c(0, 2, 4, 5), 
                      min.pct = 0.5, only.pos = TRUE)
Marker.cont<-FindMarkers(object = Mag.combined, ident.1 = c(0), ident.2 = c(1, 2, 3, 4, 5), 
                         min.pct = 0.5, only.pos = TRUE)

#### look for constant markers - low dispersion, respectable expression level ###
Mag_genes_dispersion<-data.frame(gene=row.names(Mag.combined@raw.data),
           gene.mean= Mag.combined@hvg.info$gene.mean,
           gene.dispersion=Mag.combined@hvg.info$gene.dispersion.scaled,
           non_zero=apply(Mag.combined@raw.data, 1, function(x) sum(x>0)),
           total_cell=apply(Mag.combined@raw.data, 1, function(x) sum(x>=0)) )
#### look at dispersion of traditionally regarded as house keeping genes ###
### E Eisenberg - 2013
HK_genes_2013<-c( "CHMP2A", "EMC7", "GPI", "PSMB2", "PSMB4", "RAB7A", "REEP5", "SNRPD3", "VCP", "VPS29")
### VLT Hoang - 2017
HK_genes_2017<-c("RPS13", "RPL7A", "EEF1B2", "RPS27A", "RPLP0", "RPL38", "EEF1A1", "RPL11", "RPL9", "GAPDH", "RPL23", "HPRT1", "ACTB")
###Tatiana M. Tilli BMC genomics
HK_genes_2016<-c("DHX9", "MZT2B", "UBXN4", "LARP1", "TAF2", "CCSER2", "STX5", "SYMPK", "TMEM11")
Mag_genes_dispersion[c(HK_genes_2013, HK_genes_2017, HK_genes_2016),]
Mag_genes_dispersion[c(HK_genes_2013, HK_genes_2017, HK_genes_2016),] %>%
  dplyr::filter(gene.dispersion <0.1, non_zero >4000)
### RPS13, RPL7A, RPLP0, RPL38, GAPDH, RPL23, ###plot.new()
par(mfrow=c(2,3)) 
for (i in HK_genes_2017[1:6]){
  hist(Mag.combined@raw.data[i,], main=i, fre=TRUE)}
dev.off()
Mag_genes_dispersion %>%
  dplyr::filter(gene.mean<1, gene.dispersion<0.05) %>%
  dplyr::filter(gene.mean>0.3, non_zero>5000) %>% dplyr::arrange(abs(gene.dispersion), desc(non_zero))
### 23 genes ###
### there are only 9 genes fulfill this criteria ###
Mag.constant<-Mag_genes_dispersion %>%
  dplyr::filter(gene.mean<1, gene.dispersion<0.05) %>%
  dplyr::filter(gene.mean>0.3, non_zero>5000) %>% dplyr::arrange(abs(gene.dispersion), desc(non_zero)) %>%
  dplyr::select(gene) %>% ungroup()
Mag.constant<-lapply(Mag.constant, as.character)
Mag.constant<-unlist(Mag.constant)
Mag.constant<-c(Mag.constant, "RPS13", "RPL7A", "RPLP0", "RPL38", "GAPDH", "RPL23")
### double check constant genes, use it to make PCA###

constant_pca <- RunPCA(object = Mag.combined, pc.genes =  Mag.constant, do.print = TRUE, pcs.print = 1:5, 
                       genes.print = 5)
PCAPlot(object = constant_pca, dim.1 = 1, dim.2 = 2)
PCElbowPlot(object = constant_pca) + scale_y_continuous(breaks=seq(1, 9, 1))
JackStraw(object = constant_pca, num.replicate = 100, display.progress = FALSE)
JackStrawPlot(object = constant_pca, PCs = 1:12)
###  fail PCA, so I guess these are useful genes ###
par(mfrow=c(2,2)) 
plot(density(Mag.combined@raw.data["TPM3",]) , main="TPM3")
plot(density(Mag.combined@raw.data["RPS13",]), main="RPS13")
hist(Mag.combined@raw.data["TPM3",], main="TPM3", freq=TRUE)
hist(Mag.combined@raw.data["RPS13",], main="RPS13", freq=TRUE)

### find useful marker ###
Apop_marker<-row.names(Marker.apop[Marker.apop$p_val_adj<0.01 & Marker.apop$avg_logFC >0.5, ])
### 20 genes ###
imat_marker<-row.names(Marker.imat[Marker.imat$p_val_adj<0.01 & Marker.imat$avg_logFC >0.5, ])
### 18 genes ###
cont_marker<-row.names(Marker.cont[Marker.cont$p_val_adj<0.01 & Marker.cont$avg_logFC >0.5 &
                        Marker.cont$pct.1 - Marker.cont$pct.2 >0.3, ])
### 17 genes ###

