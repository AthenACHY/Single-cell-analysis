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
  testing_apop_D<-apply(testing_apop_p, 1, sum)
  testing_imat_D<-apply(testing_Imat_p, 1, sum)
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
find_identity<-function(D_tab){
D_tab$D[is.na(D_tab$D)]<-0
D_tab<- D_tab %>% spread(group, D)
D_tab$est<-colnames(D_tab[3:5])[max.col(D_tab[3:5],ties.method="random")]
return(D_tab)
}
### wrap into the function to compute consensus identity ###
creat_consensus1p<-function(Aggr, ...){
  Apopmarker_seurat_D<-make_D1p_score(Aggr, Apop_marker, HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells)
  interested_geen1_D<-make_D1p_score(Aggr, interested_genes1, HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells)
  interested_geen2_D<-make_D1p_score(Aggr, interested_genes2, HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells)
  interested_geen3_D<-make_D1p_score(Aggr, interested_genes3, HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells)
  interested_geen4_D<-make_D1p_score(Aggr, interested_genes4, HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells)
  Apop_ident<-find_identity(Apopmarker_seurat_D)
  int1_ident<-find_identity(interested_geen1_D)
  int2_ident<-find_identity(interested_geen2_D)
  int3_ident<-find_identity(interested_geen3_D)
  int4_ident<-find_identity(interested_geen4_D)
  consensus_tab<-left_join(Apop_ident[, c("cells", "est")], int1_ident[, c("cells", "est")], by="cells")
  consensus_tab<-left_join(consensus_tab, int2_ident[, c("cells", "est")], by="cells")
  consensus_tab<-left_join(consensus_tab, int3_ident[, c("cells", "est")], by="cells")
  consensus_tab<-left_join(consensus_tab, int4_ident[, c("cells", "est")], by="cells")
#### vote for consensus for each cell ###
  row.names(consensus_tab)<-consensus_tab$cells
### convert each row to a vector, then look for the most frequent estimation ###  
  consensus_tab$ultimate_est<-apply(consensus_tab[, 2:ncol(consensus_tab)], 1, function(x) names(which.max(table(unlist(x)))))
### return final orig_ultimat est_tab for evlauation ###
  consensus_tab<-right_join(data.frame(unique(as.matrix(interested_geen1_D[, c("cells", "orig")]))), 
                          consensus_tab[, c("cells", "ultimate_est")], by="cells")
  return(consensus_tab)
}

creat_consensus<-function(Aggr, ...){
  Apopmarker_seurat_D<-make_D_score(Aggr, Apop_marker, HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells)
  interested_geen1_D<-make_D_score(Aggr, interested_genes1, HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells)
  interested_geen2_D<-make_D_score(Aggr, interested_genes2, HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells)
  interested_geen3_D<-make_D_score(Aggr, interested_genes3, HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells)
  interested_geen4_D<-make_D_score(Aggr, interested_genes4, HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells)
  Apop_ident<-find_identity(Apopmarker_seurat_D)
  int1_ident<-find_identity(interested_geen1_D)
  int2_ident<-find_identity(interested_geen2_D)
  int3_ident<-find_identity(interested_geen3_D)
  int4_ident<-find_identity(interested_geen4_D)
  consensus_tab<-left_join(Apop_ident[, c("cells", "est")], int1_ident[, c("cells", "est")], by="cells")
  consensus_tab<-left_join(consensus_tab, int2_ident[, c("cells", "est")], by="cells")
  consensus_tab<-left_join(consensus_tab, int3_ident[, c("cells", "est")], by="cells")
  consensus_tab<-left_join(consensus_tab, int4_ident[, c("cells", "est")], by="cells")
  #### vote for consensus for each cell ###
  row.names(consensus_tab)<-consensus_tab$cells
  ### convert each row to a vector, then look for the most frequent estimation ###  
  consensus_tab$ultimate_est<-apply(consensus_tab[, 2:ncol(consensus_tab)], 1, function(x) names(which.max(table(unlist(x)))))
  ### return final orig_ultimat est_tab for evlauation ###
  consensus_tab<-right_join(data.frame(unique(as.matrix(interested_geen1_D[, c("cells", "orig")]))), 
                            consensus_tab[, c("cells", "ultimate_est")], by="cells")
  return(consensus_tab)
}
### make confusion_matrix from consensus_tab
confusion_matrix_consesus_methods<-function(consensus_tab, ...){
  confusion_m<-consensus_tab %>% group_by(ultimate_est, orig) %>% summarize(count=n()) %>% spread(orig, count)
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
### run iteration ###
performance_con<-data.frame(matrix(nrow=0, ncol=6))
for (iter in (1:10)) {
  Apop_training_cells<-sample(Apop_cells, 360, replace = F)
  Imat_training_cells<-sample(Imat_cells, 360, replace = F)
  Cont_training_cells<-sample(Cont_cells, 360, replace = F)
  training_cells<-c(Apop_training_cells, Imat_training_cells, Cont_training_cells)
  testing_cells<-c(sample(setdiff(Apop_cells, training_cells), 80, replace=F), sample(setdiff(Imat_cells, training_cells), 80, replace=F), sample(setdiff(Cont_cells, training_cells), 80, replace=F))
  consensus_tab<-creat_consensus(Aggr, Apop_marker, interested_genes1, interested_genes2, interested_genes3, interested_genes4, HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells)
  cm<-confusion_matrix_consesus_methods(consensus_tab)
  performance_con<-rbind(performance_con, data.frame(predictor_gene=c("Consensus"), iteration=iter, cm))
}
  
bind_rows(performance_df, performance_con) %>% 
  group_by(model, predictor_gene) %>%
  summarize(Avg_Sens=mean(as.numeric(Sensitivity)), Avg_Spec=mean(as.numeric(Specificity)), Sd_sen=sd(as.numeric(Sensitivity)), Sd_spec=sd(as.numeric(Specificity))) %>%
  ggplot(aes(y=Avg_Sens, x=Avg_Spec, color=predictor_gene)) +
  geom_point(size=2)+
  geom_crossbar(aes(ymin=Avg_Sens-Sd_sen, ymax=Avg_Sens+Sd_sen))+  
  geom_errorbarh(aes(xmin=Avg_Spec-Sd_spec, xmax=Avg_Spec+Sd_spec))+
  facet_grid(.~model)+
  ylab("Sensitivity")+
  xlab("Specificity")+
  coord_fixed(xlim=c(0, 1), ylim=c(0, 1))+
  theme(panel.background = element_rect(fill= NA, colour= "black"),
        panel.grid.major = element_line(colour = "grey"))

performance_con25<-data.frame(matrix(nrow=0, ncol=6))
for (iter in (1:10)) {
  Apop_training_cells<-sample(Apop_cells, 360, replace = F)
  Imat_training_cells<-sample(Imat_cells, 360, replace = F)
  Cont_training_cells<-sample(Cont_cells, 360, replace = F)
  training_cells<-c(Apop_training_cells, Imat_training_cells, Cont_training_cells)
  testing_cells<-c(sample(setdiff(Apop_cells, training_cells), 80, replace=F), sample(setdiff(Imat_cells, training_cells), 80, replace=F), sample(setdiff(Cont_cells, training_cells), 80, replace=F))
  consensus_tab<-creat_consensus(Aggr, Apop_marker, interested_genes1[1:25], interested_genes2[1:25], interested_genes3[1:25], interested_genes4[1:25], HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells)
  cm<-confusion_matrix_consesus_methods(consensus_tab)
  performance_con25<-rbind(performance_con25, data.frame(predictor_gene=c("Consensus"), iteration=iter, cm))
}

performance_con10<-data.frame(matrix(nrow=0, ncol=6))
for (iter in (1:10)) {
  Apop_training_cells<-sample(Apop_cells, 360, replace = F)
  Imat_training_cells<-sample(Imat_cells, 360, replace = F)
  Cont_training_cells<-sample(Cont_cells, 360, replace = F)
  training_cells<-c(Apop_training_cells, Imat_training_cells, Cont_training_cells)
  testing_cells<-c(sample(setdiff(Apop_cells, training_cells), 80, replace=F), sample(setdiff(Imat_cells, training_cells), 80, replace=F), sample(setdiff(Cont_cells, training_cells), 80, replace=F))
  consensus_tab<-creat_consensus(Aggr, Apop_marker, interested_genes1[1:10], interested_genes2[1:10], interested_genes3[1:10], interested_genes4[1:10], HK_gene, Apop_training_cells, Imat_training_cells, Cont_training_cells, testing_cells)
  cm<-confusion_matrix_consesus_methods(consensus_tab)
  performance_con10<-rbind(performance_con10, data.frame(predictor_gene=c("Consensus"), iteration=iter, cm))
}
bind_rows(bind_rows(performance_10R %>% dplyr::mutate(gene_used=10), 
                    performance_25R %>% dplyr::mutate(gene_used=25),
                    performance_df %>% dplyr::mutate(gene_used=50)) %>%
            dplyr::filter(grepl("Clust", predictor_gene)),
bind_rows(performance_con %>% dplyr::mutate(gene_used=50),
          performance_con25 %>% dplyr::mutate(gene_used=25),
          performance_con10 %>% dplyr::mutate(gene_used=10)) )%>%
  group_by(model, predictor_gene, gene_used) %>%
  summarize(Avg_Sens=mean(as.numeric(Sensitivity)), Avg_Spec=mean(as.numeric(Specificity)), Sd_sen=sd(as.numeric(Sensitivity)), Sd_spec=sd(as.numeric(Specificity))) %>%
  ggplot(aes(y=Avg_Sens, x=Avg_Spec, shape=predictor_gene, color=gene_used)) +
  geom_point(size=4)+
  geom_crossbar(aes(ymin=Avg_Sens-Sd_sen, ymax=Avg_Sens+Sd_sen))+  facet_grid(.~model)+
  geom_errorbarh(aes(xmin=Avg_Spec-Sd_spec, xmax=Avg_Spec+Sd_spec))+
    ylab("Sensitivity")+
  xlab("Specificity")+
  coord_fixed(xlim=c(0, 1), ylim=c(0, 1))+
  scale_color_gradient(low="blue", high="red")+
  theme(panel.background = element_rect(fill= NA, colour= "black"),
        panel.grid.major = element_line(colour = "grey"))


bind_rows(
data.frame(as.matrix(t(Aggr@raw.data[interested_genes1[1:10], Apop_cells]))) %>%
  gather(gene, exp, 1:10) %>% dplyr::mutate(group=c("Apop")),
data.frame(as.matrix(t(Aggr@raw.data[interested_genes1[1:10], Imat_cells]))) %>%
  gather(gene, exp, 1:10) %>% dplyr::mutate(group=c("Imat")),
data.frame(as.matrix(t(Aggr@raw.data[interested_genes1[1:10], Cont_cells]))) %>%
  gather(gene, exp, 1:10) %>% dplyr::mutate(group=c("Cont"))
) %>%
  ggplot(aes(x=gene, y=exp)) +
  geom_boxplot()+
  facet_grid(group~.)+
  scale_y_log10()


bind_rows(
  data.frame(as.matrix(t(Aggr@raw.data[interested_genes2[1:10], Apop_cells]))) %>%
    gather(gene, exp, 1:10) %>% dplyr::mutate(group=c("Apop")),
  data.frame(as.matrix(t(Aggr@raw.data[interested_genes2[1:10], Imat_cells]))) %>%
    gather(gene, exp, 1:10) %>% dplyr::mutate(group=c("Imat")),
  data.frame(as.matrix(t(Aggr@raw.data[interested_genes2[1:10], Cont_cells]))) %>%
    gather(gene, exp, 1:10) %>% dplyr::mutate(group=c("Cont"))
) %>%
  ggplot(aes(x=gene, y=exp)) +
  geom_boxplot()+
  facet_grid(group~.)+
  scale_y_log10()


bind_rows(
  data.frame(as.matrix(t(Aggr@raw.data[interested_genes3[1:10], Apop_cells]))) %>%
    gather(gene, exp, 1:10) %>% dplyr::mutate(group=c("Apop")),
  data.frame(as.matrix(t(Aggr@raw.data[interested_genes3[1:10], Imat_cells]))) %>%
    gather(gene, exp, 1:10) %>% dplyr::mutate(group=c("Imat")),
  data.frame(as.matrix(t(Aggr@raw.data[interested_genes3[1:10], Cont_cells]))) %>%
    gather(gene, exp, 1:10) %>% dplyr::mutate(group=c("Cont"))
) %>%
  ggplot(aes(x=gene, y=exp)) +
  geom_boxplot()+
  facet_grid(group~.)+
  scale_y_log10()


bind_rows(
  data.frame(as.matrix(t(Aggr@raw.data[c(interested_genes4[1:10], "MT-ND3"), Apop_cells]))) %>%
    gather(gene, exp, 1:11) %>% dplyr::mutate(group=c("Apop")),
  data.frame(as.matrix(t(Aggr@raw.data[c(interested_genes4[1:10], "MT-ND3"), Imat_cells]))) %>%
    gather(gene, exp, 1:11) %>% dplyr::mutate(group=c("Imat")),
  data.frame(as.matrix(t(Aggr@raw.data[c(interested_genes4[1:10], "MT-ND3"), Cont_cells]))) %>%
    gather(gene, exp, 1:11) %>% dplyr::mutate(group=c("Cont"))
) %>%
  ggplot(aes(x=gene, y=exp)) +
  geom_boxplot()+
  facet_grid(group~.)+
  scale_y_log10()


bind_rows(
  data.frame(as.matrix(t(Aggr@raw.data[HK_gene, Apop_cells]))) %>%
    gather(gene, exp, 1:15) %>% dplyr::mutate(group=c("Apop")),
  data.frame(as.matrix(t(Aggr@raw.data[HK_gene, Imat_cells]))) %>%
    gather(gene, exp, 1:15) %>% dplyr::mutate(group=c("Imat")),
  data.frame(as.matrix(t(Aggr@raw.data[HK_gene, Cont_cells]))) %>%
    gather(gene, exp, 1:15) %>% dplyr::mutate(group=c("Cont"))
) %>%
  ggplot(aes(x=gene, y=exp)) +
  geom_boxplot()+
  facet_grid(group~.)+
  scale_y_log10()+
  theme(axis.text.x=element_text(angle = 45, vjust=0.5))
