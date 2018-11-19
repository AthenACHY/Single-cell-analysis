### summarize all four different way to look into the expression data ###
### raw count/HK genes, then log2 transform###
marker_SC_tab<-bind_rows(
  data.frame(cell=colnames(apsoc@raw.data), group=rep(c("Apop"), length(colnames(apsoc@raw.data))), as.matrix(t(apsoc@raw.data[marker_genes,]/apply(apsoc@raw.data[HK_gene, ], 2, mean)))),
  data.frame(cell=colnames(rawImat@raw.data), group=rep(c("Imat"), length(colnames(rawImat@raw.data))), as.matrix(t(rawImat@raw.data[marker_genes,]/apply(rawImat@raw.data[HK_gene, ], 2, mean)))),
  data.frame(cell=colnames(rawCont@raw.data), group=rep(c("Cont"), length(colnames(rawCont@raw.data))), as.matrix(t(rawCont@raw.data[marker_genes,]/apply(rawCont@raw.data[HK_gene, ], 2, mean))))
)
marker_SC_tab[, 3:188]<-sapply(marker_SC_tab[, 3:188], function(x) log2(x))

### log2 transformed count, scaled with lib size, then - HK genes ###
marker_SC_tab<-bind_rows(
  data.frame(cell=colnames(apsoc@data), group=rep(c("Apop"), length(colnames(apsoc@data))), as.matrix(t(apsoc@data[marker_genes,]-apply(apsoc@data[HK_gene, ], 2, mean)))),
  data.frame(cell=colnames(rawImat@data), group=rep(c("Imat"), length(colnames(rawImat@data))), as.matrix(t(rawImat@data[marker_genes,]-apply(rawImat@data[HK_gene, ], 2, mean)))),
  data.frame(cell=colnames(rawCont@data), group=rep(c("Cont"), length(colnames(rawCont@data))), as.matrix(t(rawCont@data[marker_genes,]-apply(rawCont@data[HK_gene, ], 2, mean))))
)

### log2 transformed count, scaled with lib size, mean(markers) - mean (HK genes) ###
marker_SC_mean_tab<-bind_rows(
  data.frame(cell=colnames(apsoc@data), group=rep(c("Apop"), length(colnames(apsoc@data))), SC=apply(apsoc@data[interested_genes1$GeneName,], 2, mean)-apply(apsoc@data[HK_gene, ], 2, mean)),
  data.frame(cell=colnames(rawImat@data), group=rep(c("Imat"), length(colnames(rawImat@data))), SC=apply(rawImat@data[interested_genes1$GeneName,], 2, mean)-apply(rawImat@data[HK_gene, ], 2, mean)),
  data.frame(cell=colnames(rawCont@data), group=rep(c("Cont"), length(colnames(rawCont@data))), SC=apply(rawCont@data[interested_genes1$GeneName,], 2, mean)-apply(rawCont@data[HK_gene, ], 2, mean))
)

### log2 transformed count ###
marker_log_tab<-bind_rows(
  data.frame(cell=colnames(apsoc@data), group=rep(c("Apop"), length(colnames(apsoc@data))), as.matrix(t(apsoc@data[marker_genes,]))),
  data.frame(cell=colnames(rawImat@data), group=rep(c("Imat"), length(colnames(rawImat@data))), as.matrix(t(rawImat@data[marker_genes,]))),
  data.frame(cell=colnames(rawCont@data), group=rep(c("Cont"), length(colnames(rawCont@data))), as.matrix(t(rawCont@data[marker_genes,])))
)

### sort out cells identity? ###
table(sapply(rownames(Aggr@meta.data), function(x) substring(x, nchar(x)-1, nchar(x))))
Aggr@meta.data$orig.ident<-sapply(rownames(Aggr@meta.data), function(x) substring(x, nchar(x), nchar(x)))
### extract cellnames ###
Apop_cells<-rownames(Aggr@meta.data)[Aggr@meta.data$orig.ident==3]
Imat_cells<-rownames(Aggr@meta.data)[Aggr@meta.data$orig.ident==2]
Cont_cells<-rownames(Aggr@meta.data)[Aggr@meta.data$orig.ident==1]

### perform cell downsampling + 80/20 (360/80) split ###
Apop_training_cells<-sample(Apop_cells, 360, replace = F)
Imat_training_cells<-sample(Imat_cells, 360, replace = F)
Cont_training_cells<-sample(Cont_cells, 360, replace = F)
training_cells<-c(Apop_training_cells, Imat_training_cells, Cont_training_cells)
testing_cells<-c(sample(setdiff(Apop_cells, training_cells), 80, replace=F), sample(setdiff(Imat_cells, training_cells), 80, replace=F), sample(setdiff(Cont_cells, training_cells), 80, replace=F))

### build cdf for each genes ###
make_D_score<-function(genes, HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells){
training_cells<-c(Apop_training_cells, Imat_training_cells, Cont_training_cells)
### get mean expression of housekeeping genes, then log2 transform ###
HK_mean_counts<-log2(apply(Aggr@raw.data[HK_gene, training_cells], 2, mean))
test_mean_counts<-log2(apply(Aggr@raw.data[HK_gene, testing_cells], 2, mean))
### get mean expression of marker genes, then log2 transform ###



relative_gene_count<-Aggr@raw.data[genes, training_cells]
relative_gene_count<-as.matrix(relative_gene_count)
for (i in (1:nrow(relative_gene_count))){relative_gene_count[i, ]<-relative_gene_count[i, ]-HK_mean_counts}
relative_test_count<-log2(as.matrix(Aggr@raw.data[genes, testing_cells]))
for (i in (1:nrow(relative_test_count))){relative_test_count[i, ]<-relative_test_count[i, ]-test_mean_counts}



testing_apop_p<-data.frame(matrix(ncol = length(genes), nrow = 0))
testing_Imat_p<-data.frame(cells=testing_cells)
testing_Cont_p<-data.frame(cells=testing_cells)

for (i in (1:length(genes))){
### -Inf caanot be used for density ###
Apop_relative_count<-relative_gene_count[genes[i], Apop_training_cells]
Imat_relative_count<-relative_gene_count[genes[i], Imat_training_cells]
Cont_relative_count<-relative_gene_count[genes[i], Cont_training_cells]
### build ecdf, the p-value of X and more extreme = 1 - ecdf(X) ###
Apop_fit<-ecdf(Apop_relative_count[is.finite(Apop_relative_count)])
Imat_fit<-ecdf(Imat_relative_count[is.finite(Imat_relative_count)])
Cont_fit<-ecdf(Imat_relative_count[is.finite(Cont_relative_count)])
### get P values for testing data ###
### not setting a predictive range make computing probability invalid at the extremes? ###
testing_apop_p<-data.frame(testing_apop_p, 1-Apop_fit(relative_test_count[genes[i],]))
testing_Imat_p<-data.frame(testing_Imat_p, 1-Imat_fit(relative_test_count[genes[i],]))
testing_Cont_p<-data.frame(testing_Cont_p, 1-Cont_fit(relative_test_count[genes[i],]))
}
### get D score per cell per condition:
rownames(testing_apop_p)<-testing_cells
colnames(testing_apop_p)<-c("cells", genes)
row.names(testing_Imat_p)<-testing_cells
colnames(testing_Imat_p)<-c("cells", genes)
row.names(testing_Cont_p)<-testing_cells
colnames(testing_Cont_p)<-c("cells", genes)
### what to deal with missing values???? -Inf fall outside the range of prediction but D calculation withi missing values will lead to incomplete representation? ###
testing_apop_p<--1*log10(testing_apop_p)
testing_apop_D<-apply(testing_apop_p, 1, sum)
testing_imat_p<--1*log10(testing_imat_p)
testing_imat_D<-apply(testing_imat_p, 1, sum)
testing_cont_p<--1*log10(testing_cont_p)
testing_cont_D<-apply(testing_cont_p, 1, sum)
### determine cell identity ###
i=1
cell_ident_list<-c()
for (i in (1:length(testing_cells))) {
  Apop<-testing_apop_D[i]
  Imat<-testing_imat_D[i]
  Cont<-testing_cont_D[i]
  cell_ident<-c("NA")
  if(Apop>Imat & Apop >Cont){cell_ident<-c("Apop")}
  if(Imat>Apop & Imat >Cont){cell_ident<-c("Imat")}
  if(Cont>Imat & Cont > Apop){cell_ident<-c("Cont")}
  cell_ident_list<-c(cell_ident_list, cell_ident)
}
}
