### Repeat find clusster with higher stringency ###
Mag.combined<- FindVariableGenes(object =  Mag.combined, mean.function = ExpMean, dispersion.function = LogVMR, 
                                 x.low.cutoff = 0.5, x.high.cutoff = 3, y.cutoff = 1.5)
length(Mag.combined@var.genes)
Mag.combined <- ScaleData(object = Mag.combined, vars.to.regress = c("nUMI"))
Mag.combined <- RunPCA(object = Mag.combined, pc.genes =  Mag.combined@var.genes, do.print = TRUE, pcs.print = 1:5, 
                       genes.print = 5)
PCAPlot(object = Mag.combined, dim.1 = 1, dim.2 = 2)
PCElbowPlot(object = Mag.combined) 
Mag.combined <- JackStraw(object = Mag.combined, num.replicate = 100, display.progress = FALSE)
JackStrawPlot(object = Mag.combined, PCs = 1:12)
Mag.combined <- RunTSNE(object = Mag.combined, dims.use = 1:11, do.fast = TRUE)
TSNEPlot(object = Mag.combined)
Mag_ident<-Mag.combined@ident
### change clusterinf res to 0.7 ###
Mag.combined <- FindClusters(object = Mag.combined, reduction.type = "pca", dims.use = 1:11,
                             force.recalc = TRUE,
                             resolution = 0.7, print.output = 0, save.SNN = TRUE)
TSNEPlot(object = Mag.combined)
### change clusterinf res to 0.8 ###
###Mag.combined <- FindClusters(object = Mag.combined, reduction.type = "pca", dims.use = 1:11,
###                             force.recalc = TRUE,
###                             resolution = 0.8, print.output = 0, save.SNN = TRUE)
TSNEPlot(object = Mag.combined)
TSNEPlot(object = Mag.combined, color.use=Mag.combined@meta.data$orig.ident)
### imat1_imat2_markers
test_M<-FindMarkers(object = Mag.combined, ident.1 = c(1), ident.2 = c(2), 
                         min.pct = 0.7, only.pos = TRUE)
###see if these genes in variable genes ###
imat1_imat2<-intersect(row.names(test_M), Mag.combined@var.genes)
test_M %>% dplyr::mutate(gene=row.names(test_M)) %>%
  dplyr::filter(gene %in% Mag.combined@var.genes, avg_logFC>0.9) 

paste(imat1_imat2, collapse = ", ")

### cont_imat1 markers ###
test_M<-FindMarkers(object = Mag.combined, ident.1 = c(1), ident.2 = c(0), 
                    min.pct = 0.7, only.pos = TRUE)
###see if these genes in variable genes ###
cont_imat1<-intersect(row.names(test_M), Mag.combined@var.genes)
VlnPlot(Mag.combined, cont_imat1)
test_M[cont_imat1,] %>% dplyr::filter(avg_logFC>0.9)
paste(cont_imat1, collapse = ", ")


### cont0_imat2 markers ###
test_M<-FindMarkers(object = Mag.combined, ident.1 = c(2), ident.2 = c(0), 
                    min.pct = 0.7, only.pos = TRUE)
###see if these genes in variable genes ###
test_M %>% dplyr::mutate(gene=row.names(test_M)) %>%
  dplyr::filter(gene %in% Mag.combined@var.genes, avg_logFC>0.9) 

cont_imat2<-intersect(row.names(test_M), Mag.combined@var.genes)
VlnPlot(Mag.combined, cont_imat2)

paste(cont_imat2, collapse = ", ")

###Apop4 vs Cont0+3+Imat1###
test_M<-FindMarkers(object = Mag.combined, ident.1 = c(4), ident.2 = c(0, 3, 1), 
                    min.pct = 0.7, only.pos = TRUE)
###see if these genes in variable genes ###
test_M %>% dplyr::mutate(gene=row.names(test_M)) %>%
  dplyr::filter(gene %in% Mag.combined@var.genes, avg_logFC>0.9) 

Apop4_ContImat<-intersect(row.names(test_M), Mag.combined@var.genes)
VlnPlot(Mag.combined, Apop4_ContImat)

paste(Apop4_ContImat, collapse = ", ")

###Apop4 vs Imat2 ###
test_M<-FindMarkers(object = Mag.combined, ident.1 = c(4), ident.2 = c(2), 
                    min.pct = 0.7, only.pos = TRUE)
###see if these genes in variable genes ###
test_M %>% dplyr::mutate(gene=row.names(test_M)) %>%
  dplyr::filter(gene %in% Mag.combined@var.genes, avg_logFC>0.9) 

###Imat2 vs Apop4 ###
test_M<-FindMarkers(object = Mag.combined, ident.1 = c(2), ident.2 = c(4), 
                    min.pct = 0.7, only.pos = TRUE)
test_M %>% dplyr::mutate(gene=row.names(test_M)) %>%
  dplyr::filter(gene %in% Mag.combined@var.genes, avg_logFC>0.9) 


###Apop4 vs Apop6 ###
test_M<-FindMarkers(object = Mag.combined, ident.1 = c(4), ident.2 = c(6), 
                     min.pct = 0.7, only.pos = TRUE)

###see if these genes in variable genes ###
test_M %>% dplyr::mutate(gene=row.names(test_M)) %>%
  dplyr::filter(gene %in% Mag.combined@var.genes, avg_logFC>0.9) 
### recipical testing ###
test_M<-FindMarkers(object = Mag.combined, ident.1 = c(6), ident.2 = c(4), 
                     min.pct = 0.7, only.pos = TRUE)


test_M %>% dplyr::mutate(gene=row.names(test_M)) %>%
  dplyr::filter(gene %in% Mag.combined@var.genes, avg_logFC>0.9) 
### no genes in return ###
Apop4_ContImat<-intersect(row.names(test_M), Mag.combined@var.genes)
VlnPlot(Mag.combined, Apop4_ContImat)
paste(Apop4_ContImat, collapse = ", ")

### Cont 0 vs Cont3 ###
test_M<-FindMarkers(object = Mag.combined, ident.1 = c(0), ident.2 = c(3), 
                    min.pct = 0.3 , only.pos = TRUE)
### no genes differ more than avg_log_FC of 0.5 ###

### Cont 0, 3 vs the rest ###
test_M<-FindMarkers(object = Mag.combined, ident.1 = c(0, 3), ident.2 = c(1, 2, 4, 6), 
                    min.pct = 0.3 , only.pos = TRUE)
test_M %>% dplyr::mutate(gene=row.names(test_M)) %>%
  dplyr::filter(gene %in% Mag.combined@var.genes) 

### Apop7 vs the rest? ###
test_M<-FindMarkers(object = Mag.combined, ident.1 = c(7), ident.2 = c(0, 1, 2, 3, 4, 5, 6), 
                    min.pct = 0.3 , only.pos = TRUE)

test_M %>% dplyr::mutate(gene=row.names(test_M)) %>%
  dplyr::filter(gene %in% Mag.combined@var.genes, avg_logFC>0.9) 

### are these outliners? extreme high expression for a lot of genes... ###
clus7cells<-attributes(Mag.combined@ident[Mag.combined@ident==7])$names
Mag.combined@meta.data[Mag.combined@meta.data$res.0.7==c("7"),]
Mag.combined@meta.data %>% ggplot(aes(x=nGene, y=nUMI, color=orig.ident)) +
  geom_point()+
  facet_wrap(~res.0.7, ncol=3)+
  scale_x_continuous(breaks=c(1000, 3000, 5000), labels = c("1k", "3k", "5k"))


Mag.combined@meta.data %>% ggplot(aes(x=res.0.7, fill=orig.ident)) +
  geom_bar(stat = "count", color="black")+
  geom_text(data=Mag.combined@meta.data %>% group_by(res.0.7, orig.ident) %>% summarise(count=n()), aes(y=count,label=count), 
  vjust=-0.2)+
  ylim(0, 2300)+
  facet_grid(orig.ident~.)+
  ylab("Number of cells")+
  xlab("cluster number at res.0.7")+
  theme(legend.position = "none")

### the Cont03 markers have avg_logFC, < 0.5 ###
### pick markers not in var.genes instead??
#### deeper levels of finding DE genes ###
### Cont cells in cluster1 and Imat cells in cluster 1 ###
Pop_clures07<-tidyr::unite(Mag.combined@meta.data, "p_c", c("orig.ident", "res.0.7"), sep="_")$p_c
Pop_clures07<-factor(Pop_clures07, levels=unique(Pop_clures07))
attributes(Pop_clures07)$names<-attributes(Mag.combined@ident)$names
Mag.combined@ident<-Pop_clures07
### compare cells within cluster 1###
test_M<-FindMarkers(object = Mag.combined, ident.1 = c("Imat_1"), ident.2 = c("Cont_1"), 
            min.pct = 0.1 , only.pos = TRUE)
### Only HSPA1A is over avg_Log_FC >0.5 ###
VlnPlot(Mag.combined, c("HSPA1A"))

VlnPlot(Mag.combined, c("PKM", "ACTG1", "TMSB10", "SLC25A37", "ALAS2", "HBA1", "MT-ND5", "MT-ND3", "UBC", "SAT1"))

hist(Mag.combined@data["HSPA1A", Mag.combined@ident==1])
hist(Mag.combined@raw.data["HSPA1A", Mag.combined@ident==1])
### find cluster enriched genes ###
### 0+3 vs the rest 
test_M<-FindMarkers(object = Mag.combined, ident.1 = c(0,3), ident.2=c(1, 2, 4, 5, 6, 7),
                    min.pct = 0.2)
### no cluster markers in var.genes ###
VlnPlot(Mag.combined, c("NME1", "HIST1H4C"))
### 1 aginst the rest ###
test_M<-FindMarkers(object = Mag.combined, ident.1 = c(1), ident.2=c(0, 2, 3, 4, 5, 6, 7),
                    min.pct = 0.7)

### 2 aginst the rest ###
test_M<-FindMarkers(object = Mag.combined, ident.1 = c(2), ident.2=c(0, 1, 3, 4, 5, 6, 7),
                    min.pct = 0.7)
test_M %>% dplyr::mutate(gene=row.names(test_M)) %>%
  dplyr::filter(gene %in% Mag.combined@var.genes, avg_logFC>0.9) 

### 4 aginst the rest ###
test_M<-FindMarkers(object = Mag.combined, ident.1 = c(4), ident.2=c(0, 1, 3, 2, 5, 6, 7),
                    min.pct = 0.7)
test_M %>% dplyr::mutate(gene=row.names(test_M)) %>%
  dplyr::filter(gene %in% Mag.combined@var.genes, avg_logFC>0.9) 

