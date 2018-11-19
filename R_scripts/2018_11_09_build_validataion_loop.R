### put function here ###
### build density function for each genes relative to HK genes###
make_D_score<-function(Aggr, genes, HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells){
  training_cells<-c(Apop_training_cells, Imat_training_cells, Cont_training_cells)
  ### get mean expression of housekeeping genes, then log2 transform ###
  train_HK_counts<-log2(apply(Aggr@raw.data[HK_gene, training_cells], 2, mean)+1)
  test_HK_counts<-log2(apply(Aggr@raw.data[HK_gene, testing_cells], 2, mean)+1)
  ### get mean expression of marker genes, then log2 transform ###
  ### log2p transform each gene ###
  relative_train_count<-as.matrix(Aggr@raw.data[genes, training_cells])
  relative_train_count<-relative_train_count+1
  relative_train_count<-log2(relative_train_count)
  ### ###
  relative_test_count<-as.matrix(Aggr@raw.data[genes, testing_cells])
  relative_test_count<-relative_test_count + 1
  relative_test_count<-log2(relative_test_count)
  ### normalise with HK genes 
  for (i in (1:nrow(relative_train_count))){relative_train_count[i, ]<-relative_train_count[i, ]-train_HK_counts}
  for (i in (1:nrow(relative_test_count))){relative_test_count[i, ]<-relative_test_count[i, ]-test_HK_counts}
  ### make a plot to look at the distribution ###
  d<-data.frame(t(relative_train_count)) %>% dplyr::mutate(names=row.names(t(relative_train_count))) %>%
    dplyr::mutate(group=substring(names, nchar(names), nchar(names))) %>%
    gather(gene, exp, 1:length(genes)) %>%
    ggplot(aes(x=group, y=exp, colour=group)) +
    geom_boxplot()+
    facet_wrap(~gene, ncol=5)
  ### make empty df for probability index ###
  testing_apop_p<-data.frame(matrix(ncol = length(genes), nrow = length(testing_cells)))
  testing_Imat_p<-data.frame(matrix(ncol = length(genes), nrow = length(testing_cells)))
  testing_Cont_p<-data.frame(matrix(ncol = length(genes), nrow = length(testing_cells)))
  ### do the modelling and get Pr(X>Xg) ###
  for (i in (1:length(genes))){
    Apop_relative_count<-relative_train_count[genes[i], Apop_training_cells]
    Imat_relative_count<-relative_train_count[genes[i], Imat_training_cells]
    Cont_relative_count<-relative_train_count[genes[i], Cont_training_cells]
    ### build empirical pdf, the p-value of X and more extreme = 1 - ecdf(X) ###
    Apop_fit<-approxfun(density(Apop_relative_count[is.finite(Apop_relative_count)]))
    Imat_fit<-approxfun(density(Imat_relative_count[is.finite(Imat_relative_count)]))
    Cont_fit<-approxfun(density(Cont_relative_count[is.finite(Cont_relative_count)]))
    ### get P values for testing data ###
    ### not setting a predictive range make computing probability invalid at the extremes? ###
    testing_apop_p[, i]<-Apop_fit(relative_test_count[genes[i],])
    testing_Imat_p[, i]<-Imat_fit(relative_test_count[genes[i],])
    testing_Cont_p[, i]<-Cont_fit(relative_test_count[genes[i],])
  }
  ### get D score per cell per condition:
  rownames(testing_apop_p)<-testing_cells
  colnames(testing_apop_p)<-c(genes)
  row.names(testing_Imat_p)<-testing_cells
  colnames(testing_Imat_p)<-c(genes)
  row.names(testing_Cont_p)<-testing_cells
  colnames(testing_Cont_p)<-c(genes)
  ### what to deal with missing values???? -Inf fall outside the range of prediction but D calculation withi missing values will lead to incomplete representation? ###
  ##testing_apop_p<--1*log10(testing_apop_p)
  testing_apop_D<-apply(testing_apop_p, 1, sum)
  ###testing_Imat_p<--1*log10(testing_Imat_p)
  testing_imat_D<-apply(testing_Imat_p, 1, sum)
  ###testing_Cont_p<--1*log10(testing_Cont_p)
  testing_cont_D<-apply(testing_Cont_p, 1, sum)
  ### determine cell identity .... hold on lets look at the D distribution ###
  D_tab<-bind_rows(data.frame(D=testing_apop_D, cells=testing_cells, group=c("Apop")),
                   data.frame(D=testing_imat_D, cells=testing_cells, group=c("Imat")),
                   data.frame(D=testing_cont_D, cells=testing_cells, group=c("Cont")))
  D_tab <- D_tab %>% dplyr::mutate(orig=substring(as.character(cells), nchar(as.character(cells)), nchar(as.character(cells))))
  D_tab$orig[D_tab$orig==c("3")]<-c("Apop")
  D_tab$orig[D_tab$orig==c("2")]<-c("Imat")
  D_tab$orig[D_tab$orig==c("1")]<-c("Cont")
  g<- D_tab %>%  ggplot(aes(x=orig, y=D, colour=orig)) +
    geom_boxplot()+
    facet_wrap(~group, ncol=3)
  
  gridExtra::grid.arrange(d, newpage = T)
  gridExtra::grid.arrange(g, newpage = T)
  return(D_tab)
}
### build density function for each genes relative to library size? use Seurat normalised @data ###
make_DLB_score<-function(Aggr, genes, HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells){
  training_cells<-c(Apop_training_cells, Imat_training_cells, Cont_training_cells)
  ### get mean expression of marker genes, then log2 transform ###
  ### log2p transform each gene ###
  relative_train_count<-as.matrix(Aggr@data[genes, training_cells])
  ### ###
  relative_test_count<-as.matrix(Aggr@data[genes, testing_cells])
  ### make a plot to look at the distribution ###
  d<-data.frame(t(relative_train_count)) %>% dplyr::mutate(names=row.names(t(relative_train_count))) %>%
    dplyr::mutate(group=substring(names, nchar(names), nchar(names))) %>%
    gather(gene, exp, 1:length(genes)) %>%
    ggplot(aes(x=group, y=exp, colour=group)) +
    geom_boxplot()+
    facet_wrap(~gene, ncol=5)
  ### make empty df for probability index ###
  testing_apop_p<-data.frame(matrix(ncol = length(genes), nrow = length(testing_cells)))
  testing_Imat_p<-data.frame(matrix(ncol = length(genes), nrow = length(testing_cells)))
  testing_Cont_p<-data.frame(matrix(ncol = length(genes), nrow = length(testing_cells)))
  ### do the modelling and get Pr(X>Xg) ###
  for (i in (1:length(genes))){
    Apop_relative_count<-relative_train_count[genes[i], Apop_training_cells]
    Imat_relative_count<-relative_train_count[genes[i], Imat_training_cells]
    Cont_relative_count<-relative_train_count[genes[i], Cont_training_cells]
    ### build empirical distribution ###
    Apop_fit<-approxfun(density(Apop_relative_count[is.finite(Apop_relative_count)]))
    Imat_fit<-approxfun(density(Imat_relative_count[is.finite(Imat_relative_count)]))
    Cont_fit<-approxfun(density(Cont_relative_count[is.finite(Cont_relative_count)]))
    ### get P values for testing data ###
    ### not setting a predictive range make computing probability invalid at the extremes? ###
    testing_apop_p[, i]<-Apop_fit(relative_test_count[genes[i],])
    testing_Imat_p[, i]<-Imat_fit(relative_test_count[genes[i],])
    testing_Cont_p[, i]<-Cont_fit(relative_test_count[genes[i],])
  }
  ### get D score per cell per condition:
  rownames(testing_apop_p)<-testing_cells
  colnames(testing_apop_p)<-c(genes)
  row.names(testing_Imat_p)<-testing_cells
  colnames(testing_Imat_p)<-c(genes)
  row.names(testing_Cont_p)<-testing_cells
  colnames(testing_Cont_p)<-c(genes)
  ### what to deal with missing values???? -Inf fall outside the range of prediction but D calculation withi missing values will lead to incomplete representation? ###
  ##testing_apop_p<--1*log10(testing_apop_p)
  testing_apop_D<-apply(testing_apop_p, 1, sum)
  ###testing_Imat_p<--1*log10(testing_Imat_p)
  testing_imat_D<-apply(testing_Imat_p, 1, sum)
  ###testing_Cont_p<--1*log10(testing_Cont_p)
  testing_cont_D<-apply(testing_Cont_p, 1, sum)
  ### determine cell identity .... hold on lets look at the D distribution ###
  D_tab<-bind_rows(data.frame(D=testing_apop_D, cells=testing_cells, group=c("Apop")),
                   data.frame(D=testing_imat_D, cells=testing_cells, group=c("Imat")),
                   data.frame(D=testing_cont_D, cells=testing_cells, group=c("Cont")))
  D_tab <- D_tab %>% dplyr::mutate(orig=substring(as.character(cells), nchar(as.character(cells)), nchar(as.character(cells))))
  D_tab$orig[D_tab$orig==c("3")]<-c("Apop")
  D_tab$orig[D_tab$orig==c("2")]<-c("Imat")
  D_tab$orig[D_tab$orig==c("1")]<-c("Cont")
  g<- D_tab %>%  ggplot(aes(x=orig, y=D, colour=orig)) +
    geom_boxplot()+
    facet_wrap(~group, ncol=3)
  
  gridExtra::grid.arrange(d, newpage = T)
  gridExtra::grid.arrange(g, newpage = T)
  return(D_tab)
}
### build density function for each gene relative to library size ##
make_D1p_score<-function(Aggr, genes, HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells){
  training_cells<-c(Apop_training_cells, Imat_training_cells, Cont_training_cells)
  ### get total library count per cells ###
  train_lib_sum<-colSums(Aggr@raw.data[, training_cells])
  test_lib_sum<-colSums(Aggr@raw.data[, testing_cells])
  ### get mean expression of marker genes, then log2 transform ###
  ### log1p transform each gene/total library size ###
  relative_train_count<-as.matrix(Aggr@raw.data[genes, training_cells])
  ### ###
  relative_test_count<-as.matrix(Aggr@raw.data[genes, testing_cells])
  ### normalise with HK genes 
  for (i in (1:nrow(relative_train_count))){
    relative_train_count[i, ]<-relative_train_count[i, ]+1
    relative_train_count[i, ]<-relative_train_count[i, ]/train_lib_sum
    relative_train_count[i, ]<-log10(relative_train_count[i, ])}
  for (i in (1:nrow(relative_test_count))){
    relative_test_count[i, ]<-relative_test_count[i, ]+1
    relative_test_count[i, ]<-relative_test_count[i, ]/test_lib_sum
    relative_test_count[i, ]<-log10(relative_test_count[i, ])}
    ### make a plot to look at the distribution ###
  d<-data.frame(t(relative_train_count)) %>% dplyr::mutate(names=row.names(t(relative_train_count))) %>%
    dplyr::mutate(group=substring(names, nchar(names), nchar(names))) %>%
    gather(gene, exp, 1:length(genes)) %>%
    ggplot(aes(x=group, y=exp, colour=group)) +
    geom_boxplot()+
    facet_wrap(~gene, ncol=5)
  ### make empty df for probability index ###
  testing_apop_p<-data.frame(matrix(ncol = length(genes), nrow = length(testing_cells)))
  testing_Imat_p<-data.frame(matrix(ncol = length(genes), nrow = length(testing_cells)))
  testing_Cont_p<-data.frame(matrix(ncol = length(genes), nrow = length(testing_cells)))
  ### do the modelling and get Pr(X>Xg) ###
  for (i in (1:length(genes))){
    Apop_relative_count<-relative_train_count[genes[i], Apop_training_cells]
    Imat_relative_count<-relative_train_count[genes[i], Imat_training_cells]
    Cont_relative_count<-relative_train_count[genes[i], Cont_training_cells]
    ### build empirical pdf, the p-value of X and more extreme = 1 - ecdf(X) ###
    Apop_fit<-approxfun(density(Apop_relative_count[is.finite(Apop_relative_count)]))
    Imat_fit<-approxfun(density(Imat_relative_count[is.finite(Imat_relative_count)]))
    Cont_fit<-approxfun(density(Cont_relative_count[is.finite(Cont_relative_count)]))
    ### get P values for testing data ###
    ### not setting a predictive range make computing probability invalid at the extremes? ###
    testing_apop_p[, i]<-Apop_fit(relative_test_count[genes[i],])
    testing_Imat_p[, i]<-Imat_fit(relative_test_count[genes[i],])
    testing_Cont_p[, i]<-Cont_fit(relative_test_count[genes[i],])
  }
  ### get D score per cell per condition:
  rownames(testing_apop_p)<-testing_cells
  colnames(testing_apop_p)<-c(genes)
  row.names(testing_Imat_p)<-testing_cells
  colnames(testing_Imat_p)<-c(genes)
  row.names(testing_Cont_p)<-testing_cells
  colnames(testing_Cont_p)<-c(genes)
  ### what to deal with missing values???? -Inf fall outside the range of prediction but D calculation withi missing values will lead to incomplete representation? ###
  ##testing_apop_p<--1*log10(testing_apop_p)
  testing_apop_D<-apply(testing_apop_p, 1, sum)
  ###testing_Imat_p<--1*log10(testing_Imat_p)
  testing_imat_D<-apply(testing_Imat_p, 1, sum)
  ###testing_Cont_p<--1*log10(testing_Cont_p)
  testing_cont_D<-apply(testing_Cont_p, 1, sum)
  ### determine cell identity .... hold on lets look at the D distribution ###
  D_tab<-bind_rows(data.frame(D=testing_apop_D, cells=testing_cells, group=c("Apop")),
                   data.frame(D=testing_imat_D, cells=testing_cells, group=c("Imat")),
                   data.frame(D=testing_cont_D, cells=testing_cells, group=c("Cont")))
  D_tab <- D_tab %>% dplyr::mutate(orig=substring(as.character(cells), nchar(as.character(cells)), nchar(as.character(cells))))
  D_tab$orig[D_tab$orig==c("3")]<-c("Apop")
  D_tab$orig[D_tab$orig==c("2")]<-c("Imat")
  D_tab$orig[D_tab$orig==c("1")]<-c("Cont")
  g<- D_tab %>%  ggplot(aes(x=orig, y=D, colour=orig)) +
    geom_boxplot()+
    facet_wrap(~group, ncol=3)
  
  gridExtra::grid.arrange(d, newpage = T)
  gridExtra::grid.arrange(g, newpage = T)
  return(D_tab)
}



### decision making step - how to pool all likelihood together? ###
### calculate confusion matrix ###
create_confusion_matrix<-function(D_tab){
  ### NA is created when one predictor gene yield undefined results - leading later calulation as NA ###
  ## We manintain the NA and decided that if one gene is missing from the overall prediction, we define that prediction as invalid - assign 0 in the log likelihood in the end ###
  D_tab$D[is.na(D_tab$D)]<-0
  D_tab<- D_tab %>% spread(group, D)
  D_tab$est<-colnames(D_tab[3:5])[max.col(D_tab[3:5],ties.method="random")]
  ### make confusion matrix ###
  confusion_m<-D_tab %>% group_by(est, orig) %>% summarize(count=n()) %>% spread(orig, count)
  confusion_m[is.na(confusion_m)]<-0
  Apop_sens<-confusion_m[1, "Apop"]/sum(confusion_m$Apop)
  Cont_sens<-confusion_m[2, "Cont"]/sum(confusion_m$Cont)
  Imat_sens<-confusion_m[3, "Imat"]/sum(confusion_m$Imat)
  Apop_spec<-(confusion_m[2, "Cont"]+ confusion_m[3, "Imat"] + confusion_m[3, "Cont"]+ confusion_m[2, "Imat"])/sum(colSums(confusion_m[, c("Cont", "Imat")]))
  Cont_spec<-(confusion_m[1, "Apop"] + confusion_m[3, "Apop"] + confusion_m[1, "Imat"] + confusion_m[3, "Imat"])/(sum(confusion_m$Apop) + sum(confusion_m$Imat))
  Imat_spec<-(confusion_m[1, "Apop"] + confusion_m[2, "Apop"] + confusion_m[1, "Cont"] + confusion_m[2, "Cont"])/(sum(confusion_m$Apop) + sum(confusion_m$Cont))
  Acc<-(confusion_m[1, "Apop"] + confusion_m[2, "Cont"] +  confusion_m[3, "Imat"])/sum(colSums(confusion_m[, c("Apop", "Imat", "Cont")]))
  print(confusion_m)
  message(paste(c("Sensitivity Apop="), round(Apop_sens *100, 2), c("%, "), c("Sepcificity Apop="), round(Apop_spec *100, 2), c("%")))
  message(paste(c("Sensitivity Cont="), round(Cont_sens *100, 2), c("%, "), c("Sepcificity Cont="), round(Cont_spec *100, 2), c("%")))
  message(paste(c("Sensitivity Imat="), round(Imat_sens *100, 2), c("%, "), c("Sepcificity Imat="), round(Imat_spec *100, 2), c("%")))
  message(paste(c("Accuracy:"), round(Acc*100, 2), c("%")))
  ### output a table row for later stats? ###
  ### modelling group, sensitivity, specificity, Acc fo the predictor markers ###
  sum_tab<-data.frame(model=c("Apop", "Cont", "Imat"))
  sum_tab$Sensitivity=c(Apop_sens, Cont_sens, Imat_sens)
  sum_tab$Specificity=c(Apop_spec, Cont_spec, Imat_spec)
  sum_tab$Accuracy=rep(Acc, nrow(sum_tab))
  return(sum_tab)
}


### sort out cells identity? ###
table(sapply(rownames(Aggr@meta.data), function(x) substring(x, nchar(x)-1, nchar(x))))
Aggr@meta.data$orig.ident<-sapply(rownames(Aggr@meta.data), function(x) substring(x, nchar(x), nchar(x)))
### extract cellnames ###
Apop_cells<-rownames(Aggr@meta.data)[Aggr@meta.data$orig.ident==3]
Imat_cells<-rownames(Aggr@meta.data)[Aggr@meta.data$orig.ident==2]
Cont_cells<-rownames(Aggr@meta.data)[Aggr@meta.data$orig.ident==1]
#### test a list of markers ###
Apop_marker<-c("PLCG2", "MT-ND5", "MT-ND3")
imat1_marker<-c("PKM", "ACTG1", "VIM")
imat2_marker<-c("ALAS2", "HBZ")
cont_marker<-c("NME1", "HIST1H4C")
csvdir=c("/Users/athchu/Documents/apop.sc/")
interested_genes1<-read.csv(paste(csvdir, c("4 Genes (1).csv"), sep = ""), stringsAsFactors = F, header = T)$GeneName
interested_genes2<-read.csv(paste(csvdir, c("4 Genes (2).csv"), sep = ""), stringsAsFactors = F, header = T)$GeneName
interested_genes3<-read.csv(paste(csvdir, c("4 Genes (3).csv"), sep = ""), stringsAsFactors = F, header = T)$GeneName
interested_genes4<-read.csv(paste(csvdir, c("4 Genes (4).csv"), sep = ""), stringsAsFactors = F, header = T)$GeneName


### perform cell downsampling + 80/20 (360/80) split ###
performance_df<-data.frame(matrix(nrow=0, ncol=6))
for (iter in (1:10)) {
Apop_training_cells<-sample(Apop_cells, 360, replace = F)
Imat_training_cells<-sample(Imat_cells, 360, replace = F)
Cont_training_cells<-sample(Cont_cells, 360, replace = F)
training_cells<-c(Apop_training_cells, Imat_training_cells, Cont_training_cells)
testing_cells<-c(sample(setdiff(Apop_cells, training_cells), 80, replace=F), sample(setdiff(Imat_cells, training_cells), 80, replace=F), sample(setdiff(Cont_cells, training_cells), 80, replace=F))
### calculate D ###
Apopmarker_seurat_D<-make_D_score(Aggr, Apop_marker, HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells)
Imat2marker_seurat_D<-make_D_score(Aggr, imat2_marker, HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells)
interested_geen1_D<-make_D_score(Aggr, interested_genes1, HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells)
interested_geen2_D<-make_D_score(Aggr, interested_genes2, HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells)
interested_geen3_D<-make_D_score(Aggr, interested_genes3, HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells)
interested_geen4_D<-make_D_score(Aggr, interested_genes4, HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells)
### stack df ###
cm<-create_confusion_matrix(Apopmarker_seurat_D)
performance_df<-rbind(performance_df, data.frame(predictor_gene=c("Apopmarkers"), iteration=iter, cm))
cm<-create_confusion_matrix(Imat2marker_seurat_D) 
performance_df<-rbind(performance_df, data.frame(predictor_gene=c("Imat2markers"), iteration=iter, cm))
cm<-create_confusion_matrix(interested_geen1_D) 
performance_df<-rbind(performance_df, data.frame(predictor_gene=c("Clust1"), iteration=iter, cm))
cm<-create_confusion_matrix(interested_geen2_D) 
performance_df<-rbind(performance_df, data.frame(predictor_gene=c("Clust2"), iteration=iter, cm))
cm<-create_confusion_matrix(interested_geen3_D) 
performance_df<-rbind(performance_df, data.frame(predictor_gene=c("Clust3"), iteration=iter, cm))
cm<-create_confusion_matrix(interested_geen4_D) 
performance_df<-rbind(performance_df, data.frame(predictor_gene=c("Clust4"), iteration=iter, cm))
rm(cm)
}
### using HK genes for normalisation is a very invalid ideas? cuz all difference are gone. ###
#### make a loop for iteration? ###
performance_df %>% 
  group_by(model, predictor_gene) %>%
  summarize(Avg_Sens=mean(as.numeric(Sensitivity)), Avg_Spec=mean(as.numeric(Specificity)), Sd_sen=sd(as.numeric(Sensitivity)), Sd_spec=sd(as.numeric(Specificity))) %>%
  ggplot(aes(y=Avg_Sens, x=Avg_Spec, color=predictor_gene)) +
  geom_point(size=3)+
  geom_crossbar(aes(xmin=Avg_Spec-Sd_spec, xmax=Avg_Spec+Sd_spec, ymin=Avg_Sens-Sd_sen, ymax=Avg_Sens+Sd_sen))+  
  facet_grid(.~model)+
  ylab("Sensitivity")+
  xlab("Specificity")+
  coord_fixed(xlim=c(0, 1), ylim=c(0, 1))+
  theme(panel.background = element_rect(fill= NA, colour= "black"),
      panel.grid.major = element_line(colour = "grey"))


###repeat with D_score_1D ###
performance_1p<-data.frame(matrix(nrow=0, ncol=6))
for (iter in (1:10)) {
  Apop_training_cells<-sample(Apop_cells, 360, replace = F)
  Imat_training_cells<-sample(Imat_cells, 360, replace = F)
  Cont_training_cells<-sample(Cont_cells, 360, replace = F)
  training_cells<-c(Apop_training_cells, Imat_training_cells, Cont_training_cells)
  testing_cells<-c(sample(setdiff(Apop_cells, training_cells), 80, replace=F), sample(setdiff(Imat_cells, training_cells), 80, replace=F), sample(setdiff(Cont_cells, training_cells), 80, replace=F))
  ### calculate D ###
  Apopmarker_seurat_D<-make_D1p_score(Aggr, Apop_marker, HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells)
  Imat2marker_seurat_D<-make_D1p_score(Aggr, imat2_marker, HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells)
  interested_geen1_D<-make_D1p_score(Aggr, interested_genes1, HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells)
  interested_geen2_D<-make_D1p_score(Aggr, interested_genes2, HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells)
  interested_geen3_D<-make_D1p_score(Aggr, interested_genes3, HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells)
  interested_geen4_D<-make_D1p_score(Aggr, interested_genes4, HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells)
  ### stack df ###
  cm<-create_confusion_matrix(Apopmarker_seurat_D)
  performance_1p<-rbind(performance_df, data.frame(predictor_gene=c("Apopmarkers"), iteration=iter, cm))
  cm<-create_confusion_matrix(Imat2marker_seurat_D) 
  performance_1p<-rbind(performance_df, data.frame(predictor_gene=c("Imat2markers"), iteration=iter, cm))
  cm<-create_confusion_matrix(interested_geen1_D) 
  performance_1p<-rbind(performance_df, data.frame(predictor_gene=c("Clust1"), iteration=iter, cm))
  cm<-create_confusion_matrix(interested_geen2_D) 
  performance_1p<-rbind(performance_df, data.frame(predictor_gene=c("Clust2"), iteration=iter, cm))
  cm<-create_confusion_matrix(interested_geen3_D) 
  performance_1p<-rbind(performance_df, data.frame(predictor_gene=c("Clust3"), iteration=iter, cm))
  cm<-create_confusion_matrix(interested_geen4_D) 
  performance_1p<-rbind(performance_df, data.frame(predictor_gene=c("Clust4"), iteration=iter, cm))
  rm(cm)
}

### look into genes set, can subsetting the gene set yield higher sensitivity? ###
### use top 25 marker genes for each cluster ###

performance_25<-data.frame(matrix(nrow=0, ncol=6))
for (iter in (1:10)) {
  Apop_training_cells<-sample(Apop_cells, 360, replace = F)
  Imat_training_cells<-sample(Imat_cells, 360, replace = F)
  Cont_training_cells<-sample(Cont_cells, 360, replace = F)
  training_cells<-c(Apop_training_cells, Imat_training_cells, Cont_training_cells)
  testing_cells<-c(sample(setdiff(Apop_cells, training_cells), 80, replace=F), sample(setdiff(Imat_cells, training_cells), 80, replace=F), sample(setdiff(Cont_cells, training_cells), 80, replace=F))
  ### calculate D ###
  interested_geen1_D<-make_D1p_score(Aggr, interested_genes1[1:25], HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells)
  interested_geen2_D<-make_D1p_score(Aggr, interested_genes2[1:25], HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells)
  interested_geen3_D<-make_D1p_score(Aggr, interested_genes3[1:25], HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells)
  interested_geen4_D<-make_D1p_score(Aggr, interested_genes4[1:25], HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells)
  ### stack df ###
  cm<-create_confusion_matrix(interested_geen1_D) 
  performance_25<-rbind(performance_25, data.frame(predictor_gene=c("Clust1"), iteration=iter, cm))
  cm<-create_confusion_matrix(interested_geen2_D) 
  performance_25<-rbind(performance_25, data.frame(predictor_gene=c("Clust2"), iteration=iter, cm))
  cm<-create_confusion_matrix(interested_geen3_D) 
  performance_25<-rbind(performance_25, data.frame(predictor_gene=c("Clust3"), iteration=iter, cm))
  cm<-create_confusion_matrix(interested_geen4_D) 
  performance_25<-rbind(performance_25, data.frame(predictor_gene=c("Clust4"), iteration=iter, cm))
  rm(cm)
}
performance_10<-data.frame(matrix(nrow=0, ncol=6))
for (iter in (1:10)) {
  Apop_training_cells<-sample(Apop_cells, 360, replace = F)
  Imat_training_cells<-sample(Imat_cells, 360, replace = F)
  Cont_training_cells<-sample(Cont_cells, 360, replace = F)
  training_cells<-c(Apop_training_cells, Imat_training_cells, Cont_training_cells)
  testing_cells<-c(sample(setdiff(Apop_cells, training_cells), 80, replace=F), sample(setdiff(Imat_cells, training_cells), 80, replace=F), sample(setdiff(Cont_cells, training_cells), 80, replace=F))
  ### calculate D ###
  interested_geen1_D<-make_D1p_score(Aggr, interested_genes1[1:10], HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells)
  interested_geen2_D<-make_D1p_score(Aggr, interested_genes2[1:10], HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells)
  interested_geen3_D<-make_D1p_score(Aggr, interested_genes3[1:10], HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells)
  interested_geen4_D<-make_D1p_score(Aggr, interested_genes4[1:10], HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells)
  ### stack df ###
  cm<-create_confusion_matrix(interested_geen1_D) 
  performance_10<-rbind(performance_10, data.frame(predictor_gene=c("Clust1"), iteration=iter, cm))
  cm<-create_confusion_matrix(interested_geen2_D) 
  performance_10<-rbind(performance_10, data.frame(predictor_gene=c("Clust2"), iteration=iter, cm))
  cm<-create_confusion_matrix(interested_geen3_D) 
  performance_10<-rbind(performance_10, data.frame(predictor_gene=c("Clust3"), iteration=iter, cm))
  cm<-create_confusion_matrix(interested_geen4_D) 
  performance_10<-rbind(performance_10, data.frame(predictor_gene=c("Clust4"), iteration=iter, cm))
  rm(cm)
}


performance_10%>% ggplot(aes(y=as.numeric(Sensitivity), x=as.numeric(Specificity), colour=predictor_gene)) +
  geom_point()+
  facet_grid(.~model)+
  ylab("Sensitivity")+
  xlab("Specificity")+
  coord_fixed(xlim=c(0, 1), ylim=c(0, 1))+
  theme(panel.background = element_rect(fill= NA, colour= "black"),
        panel.grid.major = element_line(colour = "grey"))


bind_rows(performance_10 %>% dplyr::mutate(gene_used=10), 
          performance_25 %>% dplyr::mutate(gene_used=25),
          performance_1p %>% dplyr::mutate(gene_used=50)) %>%
  dplyr::filter(grepl("Clust", predictor_gene))%>% 
  group_by(model, predictor_gene, gene_used) %>%
  summarize(Avg_Sens=mean(as.numeric(Sensitivity)), Avg_Spec=mean(as.numeric(Specificity)), Sd_sen=sd(as.numeric(Sensitivity)), Sd_spec=sd(as.numeric(Specificity))) %>%
  ggplot(aes(y=Avg_Sens, x=Avg_Spec, shape=predictor_gene, color=gene_used)) +
  geom_point(size=4)+
  geom_crossbar(aes(xmin=Avg_Spec-Sd_spec, xmax=Avg_Spec+Sd_spec, ymin=Avg_Sens-Sd_sen, ymax=Avg_Sens+Sd_sen))+  facet_grid(.~model)+
  ylab("Sensitivity")+
  xlab("Specificity")+
  coord_fixed(xlim=c(0, 1), ylim=c(0, 1))+
  scale_color_gradient(low="blue", high="red")+
  theme(panel.background = element_rect(fill= NA, colour= "black"),
        panel.grid.major = element_line(colour = "grey"))

### now use house keeping genes to repeat model ###
performance_25R<-data.frame(matrix(nrow=0, ncol=6))
for (iter in (1:10)) {
  Apop_training_cells<-sample(Apop_cells, 360, replace = F)
  Imat_training_cells<-sample(Imat_cells, 360, replace = F)
  Cont_training_cells<-sample(Cont_cells, 360, replace = F)
  training_cells<-c(Apop_training_cells, Imat_training_cells, Cont_training_cells)
  testing_cells<-c(sample(setdiff(Apop_cells, training_cells), 80, replace=F), sample(setdiff(Imat_cells, training_cells), 80, replace=F), sample(setdiff(Cont_cells, training_cells), 80, replace=F))
  ### calculate D ###
  interested_geen1_D<-make_D_score(Aggr, interested_genes1[1:25], HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells)
  interested_geen2_D<-make_D_score(Aggr, interested_genes2[1:25], HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells)
  interested_geen3_D<-make_D_score(Aggr, interested_genes3[1:25], HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells)
  interested_geen4_D<-make_D_score(Aggr, interested_genes4[1:25], HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells)
  ### stack df ###
  cm<-create_confusion_matrix(interested_geen1_D) 
  performance_25R<-rbind(performance_25R, data.frame(predictor_gene=c("Clust1"), iteration=iter, cm))
  cm<-create_confusion_matrix(interested_geen2_D) 
  performance_25R<-rbind(performance_25R, data.frame(predictor_gene=c("Clust2"), iteration=iter, cm))
  cm<-create_confusion_matrix(interested_geen3_D) 
  performance_25R<-rbind(performance_25R, data.frame(predictor_gene=c("Clust3"), iteration=iter, cm))
  cm<-create_confusion_matrix(interested_geen4_D) 
  performance_25R<-rbind(performance_25R, data.frame(predictor_gene=c("Clust4"), iteration=iter, cm))
  rm(cm)
}

performance_10R<-data.frame(matrix(nrow=0, ncol=6))
for (iter in (1:10)) {
  Apop_training_cells<-sample(Apop_cells, 360, replace = F)
  Imat_training_cells<-sample(Imat_cells, 360, replace = F)
  Cont_training_cells<-sample(Cont_cells, 360, replace = F)
  training_cells<-c(Apop_training_cells, Imat_training_cells, Cont_training_cells)
  testing_cells<-c(sample(setdiff(Apop_cells, training_cells), 80, replace=F), sample(setdiff(Imat_cells, training_cells), 80, replace=F), sample(setdiff(Cont_cells, training_cells), 80, replace=F))
  ### calculate D ###
  interested_geen1_D<-make_D_score(Aggr, interested_genes1[1:10], HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells)
  interested_geen2_D<-make_D_score(Aggr, interested_genes2[1:10], HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells)
  interested_geen3_D<-make_D_score(Aggr, interested_genes3[1:10], HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells)
  interested_geen4_D<-make_D_score(Aggr, interested_genes4[1:10], HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells)
  ### stack df ###
  cm<-create_confusion_matrix(interested_geen1_D) 
  performance_10R<-rbind(performance_10R, data.frame(predictor_gene=c("Clust1"), iteration=iter, cm))
  cm<-create_confusion_matrix(interested_geen2_D) 
  performance_10R<-rbind(performance_10R, data.frame(predictor_gene=c("Clust2"), iteration=iter, cm))
  cm<-create_confusion_matrix(interested_geen3_D) 
  performance_10R<-rbind(performance_10R, data.frame(predictor_gene=c("Clust3"), iteration=iter, cm))
  cm<-create_confusion_matrix(interested_geen4_D) 
  performance_10R<-rbind(performance_10R, data.frame(predictor_gene=c("Clust4"), iteration=iter, cm))
  rm(cm)
}

bind_rows(performance_10R %>% dplyr::mutate(gene_used=10), 
          performance_25R %>% dplyr::mutate(gene_used=25),
          performance_df %>% dplyr::mutate(gene_used=50)) %>%
  dplyr::filter(grepl("Clust", predictor_gene))%>% 
  group_by(model, predictor_gene, gene_used) %>%
  summarize(Avg_Sens=mean(as.numeric(Sensitivity)), Avg_Spec=mean(as.numeric(Specificity)), Sd_sen=sd(as.numeric(Sensitivity)), Sd_spec=sd(as.numeric(Specificity))) %>%
  ggplot(aes(y=Avg_Sens, x=Avg_Spec, shape=predictor_gene, color=gene_used)) +
  geom_point(size=4)+
  geom_crossbar(aes(xmin=Avg_Spec-Sd_spec, xmax=Avg_Spec+Sd_spec, ymin=Avg_Sens-Sd_sen, ymax=Avg_Sens+Sd_sen))+  facet_grid(.~model)+
  ylab("Sensitivity")+
  xlab("Specificity")+
  coord_fixed(xlim=c(0, 1), ylim=c(0, 1))+
  scale_color_gradient(low="blue", high="red")+
  theme(panel.background = element_rect(fill= NA, colour= "black"),
        panel.grid.major = element_line(colour = "grey"))
