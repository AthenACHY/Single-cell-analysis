### look at raw count SC for markers identified by 10x genomics ###
csvdir=c("/Users/athchu/Documents/apop.sc/")
interested_genes1<-read.csv(paste(csvdir, c("4 Genes (1).csv"), sep = ""), stringsAsFactors = F, header = T)
interested_genes2<-read.csv(paste(csvdir, c("4 Genes (2).csv"), sep = ""), stringsAsFactors = F, header = T)
interested_genes3<-read.csv(paste(csvdir, c("4 Genes (3).csv"), sep = ""), stringsAsFactors = F, header = T)
interested_genes4<-read.csv(paste(csvdir, c("4 Genes (4).csv"), sep = ""), stringsAsFactors = F, header = T)
marker_genes<-unique(c(interested_genes1$GeneName, interested_genes2$GeneName, interested_genes3$GeneName, interested_genes4$GeneName))
### calculate SC_scores for each gene ###
### draw a table : cellbarcode, orig.ident, gene, genespecific SC score ###
HK_gene<-c("ATP5E","COX7C","EEF1D","FAU","RPL10A","RPL24","RPL30","RPL35A","RPS10","RPS20", "RPS21","RPS29","SERF2","TOMM7","UQCRB")
gene_specific_SC<-c()
marker_SC_tab<-bind_rows(
  data.frame(cell=colnames(apsoc@raw.data), group=rep(c("Apop"), length(colnames(apsoc@raw.data))), as.matrix(t(apsoc@raw.data[marker_genes,]/apply(apsoc@raw.data[HK_gene, ], 2, mean)))),
  data.frame(cell=colnames(rawImat@raw.data), group=rep(c("Imat"), length(colnames(rawImat@raw.data))), as.matrix(t(rawImat@raw.data[marker_genes,]/apply(rawImat@raw.data[HK_gene, ], 2, mean)))),
  data.frame(cell=colnames(rawCont@raw.data), group=rep(c("Cont"), length(colnames(rawCont@raw.data))), as.matrix(t(rawCont@raw.data[marker_genes,]/apply(rawCont@raw.data[HK_gene, ], 2, mean))))
  )
marker_SC_tab[, 3:188]<-sapply(marker_SC_tab[, 3:188], function(x) log2(x))
marker_SC_tab<-marker_SC_tab %>% gather(gene, SC, 3:188) %>%
  dplyr::filter(SC > -10000) %>% dplyr::filter(SC < 10000)
### plot 20 plots each time ###
start<-seq(1, 188, 20)
end<-c(seq(20, 188, 20), 188)
plot_density<-function(i)
{
g<-marker_SC_tab %>% dplyr::filter(gene %in% marker_genes[start[i]:end[i]]) %>%
  ggplot(aes(x=SC, group=group, colour=group))+
  geom_density(position="identity", fill="NA")+
  facet_wrap(~gene, ncol=4)
text_df<-marker_SC_tab %>% dplyr::filter(gene %in% marker_genes[start[i]:end[i]]) %>% group_by(gene, group) %>% summarize(label=as.character(n()))
text_df$y<-0
text_df$y[text_df$group==c("Apop")]<-0
text_df$y[text_df$group==c("Imat")]<-0.05
text_df$y[text_df$group==c("Cont")]<-0.1
g+ geom_text(
  data=text_df,
  mapping=aes(x = -Inf, y = y, label = label),
  hjust   = -0.1,
  vjust   = -1
  )
}

d<-marker_SC_tab %>% dplyr::filter(gene %in% interested_genes1$GeneName) %>%
  group_by(group, cell) %>% summarise(mean_SC=mean(SC)) %>%
  ggplot(aes(x=group, y=mean_SC, colour=group))+
  ggtitle("Top 50 genes cluster 1")+
  geom_boxplot()

e<-marker_SC_tab %>% dplyr::filter(gene %in% interested_genes2$GeneName) %>%
  group_by(group, cell) %>% summarise(mean_SC=mean(SC)) %>%
  ggplot(aes(x=group, y=mean_SC, colour=group))+
  ggtitle("Top 50 genes cluster 2")+
  geom_boxplot()

f<-marker_SC_tab %>% dplyr::filter(gene %in% interested_genes3$GeneName) %>%
  group_by(group, cell) %>% summarise(mean_SC=mean(SC)) %>%
  ggplot(aes(x=group, y=mean_SC, colour=group))+
  ggtitle("Top 50 genes cluster 3")+
  geom_boxplot()

g<-marker_SC_tab %>% dplyr::filter(gene %in% interested_genes4$GeneName) %>%
  group_by(group, cell) %>% summarise(mean_SC=mean(SC)) %>%
  ggplot(aes(x=group, y=mean_SC, colour=group))+
  ggtitle("Top 50 genes cluster 4")+
  geom_boxplot()

gridExtra::grid.arrange(d,e,f,g, nrow=2)

### why Apop shift right all the time??? ###
### check if raw count of HK genes are even amongst all population ###
HKmarker_SC_tab<-bind_rows(
  data.frame(cell=colnames(apsoc@raw.data), group=rep(c("Apop"), length(colnames(apsoc@raw.data))), as.matrix(t(apsoc@raw.data[HK_gene,]))),
  data.frame(cell=colnames(rawImat@raw.data), group=rep(c("Imat"), length(colnames(rawImat@raw.data))), as.matrix(t(rawImat@raw.data[HK_gene,]))),
  data.frame(cell=colnames(rawCont@raw.data), group=rep(c("Cont"), length(colnames(rawCont@raw.data))), as.matrix(t(rawCont@raw.data[HK_gene,])))
)
HKmarker_SC_tab[, 3:17]<-sapply(HKmarker_SC_tab[, 3:17], function(x) log2(x))
HKmarker_SC_tab<-HKmarker_SC_tab %>% gather(gene, SC, 3:17) %>%
  dplyr::filter(SC > -10000) %>% dplyr::filter(SC < 10000)
g<-HKmarker_SC_tab %>%
  ggplot(aes(x=SC, group=group, colour=group))+
  geom_density(position="identity", fill="NA")+
  facet_wrap(~gene, ncol=4)

#### use normalise data for library size instead###
HK_gene<-c("ATP5E","COX7C","EEF1D","FAU","RPL10A","RPL24","RPL30","RPL35A","RPS10","RPS20", "RPS21","RPS29","SERF2","TOMM7","UQCRB")
marker_SC_tab<-bind_rows(
  data.frame(cell=colnames(apsoc@data), group=rep(c("Apop"), length(colnames(apsoc@data))), as.matrix(t(apsoc@data[marker_genes,]-apply(apsoc@data[HK_gene, ], 2, mean)))),
  data.frame(cell=colnames(rawImat@data), group=rep(c("Imat"), length(colnames(rawImat@data))), as.matrix(t(rawImat@data[marker_genes,]-apply(rawImat@data[HK_gene, ], 2, mean)))),
  data.frame(cell=colnames(rawCont@data), group=rep(c("Cont"), length(colnames(rawCont@data))), as.matrix(t(rawCont@data[marker_genes,]-apply(rawCont@data[HK_gene, ], 2, mean))))
)
marker_SC_tab<-marker_SC_tab %>% gather(gene, SC, 3:188) %>%
  dplyr::filter(SC > -10000) %>% dplyr::filter(SC < 10000)

marker_log_tab<-bind_rows(
  data.frame(cell=colnames(apsoc@data), group=rep(c("Apop"), length(colnames(apsoc@data))), as.matrix(t(apsoc@data[marker_genes,]))),
  data.frame(cell=colnames(rawImat@data), group=rep(c("Imat"), length(colnames(rawImat@data))), as.matrix(t(rawImat@data[marker_genes,]))),
  data.frame(cell=colnames(rawCont@data), group=rep(c("Cont"), length(colnames(rawCont@data))), as.matrix(t(rawCont@data[marker_genes,])))
)

marker_SC_mean_tab<-bind_rows(
  data.frame(cell=colnames(apsoc@data), group=rep(c("Apop"), length(colnames(apsoc@data))), SC=apply(apsoc@data[interested_genes1$GeneName,], 2, mean)-apply(apsoc@data[HK_gene, ], 2, mean)),
  data.frame(cell=colnames(rawImat@data), group=rep(c("Imat"), length(colnames(rawImat@data))), SC=apply(rawImat@data[interested_genes1$GeneName,], 2, mean)-apply(rawImat@data[HK_gene, ], 2, mean)),
  data.frame(cell=colnames(rawCont@data), group=rep(c("Cont"), length(colnames(rawCont@data))), SC=apply(rawCont@data[interested_genes1$GeneName,], 2, mean)-apply(rawCont@data[HK_gene, ], 2, mean))
)

d<-marker_SC_mean_tab%>% 
  group_by(group, cell) %>% summarise(mean_SC=mean(SC)) %>%
  ggplot(aes(x=group, y=mean_SC, colour=group))+
  ggtitle("Top 50 genes cluster 1")+
  geom_boxplot()

marker_SC_mean_tab<-bind_rows(
  data.frame(cell=colnames(apsoc@data), group=rep(c("Apop"), length(colnames(apsoc@data))), SC=apply(apsoc@data[interested_genes2$GeneName,], 2, mean)-apply(apsoc@data[HK_gene, ], 2, mean)),
  data.frame(cell=colnames(rawImat@data), group=rep(c("Imat"), length(colnames(rawImat@data))), SC=apply(rawImat@data[interested_genes2$GeneName,], 2, mean)-apply(rawImat@data[HK_gene, ], 2, mean)),
  data.frame(cell=colnames(rawCont@data), group=rep(c("Cont"), length(colnames(rawCont@data))), SC=apply(rawCont@data[interested_genes2$GeneName,], 2, mean)-apply(rawCont@data[HK_gene, ], 2, mean))
)

e<-marker_SC_mean_tab %>% 
  group_by(group, cell) %>% summarise(mean_SC=mean(SC)) %>%
  ggplot(aes(x=group, y=mean_SC, colour=group))+
  ggtitle("Top 50 genes cluster 2")+
  geom_boxplot()

marker_SC_mean_tab<-bind_rows(
  data.frame(cell=colnames(apsoc@data), group=rep(c("Apop"), length(colnames(apsoc@data))), SC=apply(apsoc@data[interested_genes3$GeneName,], 2, mean)-apply(apsoc@data[HK_gene, ], 2, mean)),
  data.frame(cell=colnames(rawImat@data), group=rep(c("Imat"), length(colnames(rawImat@data))), SC=apply(rawImat@data[interested_genes3$GeneName,], 2, mean)-apply(rawImat@data[HK_gene, ], 2, mean)),
  data.frame(cell=colnames(rawCont@data), group=rep(c("Cont"), length(colnames(rawCont@data))), SC=apply(rawCont@data[interested_genes3$GeneName,], 2, mean)-apply(rawCont@data[HK_gene, ], 2, mean))
)

f<-marker_SC_mean_tab %>% 
  group_by(group, cell) %>% summarise(mean_SC=mean(SC)) %>%
  ggplot(aes(x=group, y=mean_SC, colour=group))+
  ggtitle("Top 50 genes cluster 3")+
  geom_boxplot()
marker_SC_mean_tab<-bind_rows(
  data.frame(cell=colnames(apsoc@data), group=rep(c("Apop"), length(colnames(apsoc@data))), SC=apply(apsoc@data[interested_genes4$GeneName,], 2, mean)-apply(apsoc@data[HK_gene, ], 2, mean)),
  data.frame(cell=colnames(rawImat@data), group=rep(c("Imat"), length(colnames(rawImat@data))), SC=apply(rawImat@data[interested_genes4$GeneName,], 2, mean)-apply(rawImat@data[HK_gene, ], 2, mean)),
  data.frame(cell=colnames(rawCont@data), group=rep(c("Cont"), length(colnames(rawCont@data))), SC=apply(rawCont@data[interested_genes4$GeneName,], 2, mean)-apply(rawCont@data[HK_gene, ], 2, mean))
)


g<-marker_SC_mean_tab%>% 
  group_by(group, cell) %>% summarise(mean_SC=mean(SC)) %>%
  ggplot(aes(x=group, y=mean_SC, colour=group))+
  ggtitle("Top 50 genes cluster 4")+
  geom_boxplot()




gridExtra::grid.arrange(d,e,f,g, nrow=2)

marker_log_tab %>% dplyr::filter(gene %in% interested_genes1$GeneName[1:10]) %>%
  dplyr::filter(lgFC>0)%>%
  ggplot(aes(x=lgFC, group=group, colour=group))+
  geom_density(position="identity", fill="NA")+
  facet_wrap(~gene, ncol=4)


a<-marker_log_tab%>% dplyr::filter(gene %in% interested_genes1$GeneName[1:5]) %>%
  group_by(group, cell)%>% summarize(mean_exp=mean(lgFC)) %>%
  ggplot(aes(y=mean_exp, x=group, colour=group))+
  ggtitle("Top 5 genes cluster 1")+
  geom_boxplot()

b<-marker_log_tab%>% dplyr::filter(gene %in% interested_genes1$GeneName[1:10]) %>%
  group_by(group, cell)%>% summarize(mean_exp=mean(lgFC)) %>%
  ggplot(aes(y=mean_exp, x=group, colour=group))+
  ggtitle("Top 10 genes cluster 1")+
  geom_boxplot()

c<-marker_log_tab%>% dplyr::filter(gene %in% interested_genes1$GeneName[1:30]) %>%
  group_by(group, cell)%>% summarize(mean_exp=mean(lgFC)) %>%
  ggplot(aes(y=mean_exp, x=group, colour=group))+
  ggtitle("Top 30 genes cluster 1")+
  geom_boxplot()

d<-marker_log_tab%>% dplyr::filter(gene %in% interested_genes1$GeneName) %>%
  group_by(group, cell)%>% summarize(mean_exp=mean(lgFC)) %>%
  ggplot(aes(y=mean_exp, x=group, colour=group))+
  ggtitle("Top 50 genes cluster 1")+
  geom_boxplot()

e<-marker_log_tab%>% dplyr::filter(gene %in% interested_genes2$GeneName) %>%
  group_by(group, cell)%>% summarize(mean_exp=mean(lgFC)) %>%
  ggplot(aes(y=mean_exp, x=group, colour=group))+
  ggtitle("Top 50 genes cluster 2")+
  geom_boxplot()

### use markers to compute probability ###
cluster1_apop_SC

f<-marker_log_tab%>% dplyr::filter(gene %in% interested_genes3$GeneName) %>%
  group_by(group, cell)%>% summarize(mean_exp=mean(lgFC)) %>%
  ggplot(aes(y=mean_exp, x=group, colour=group))+
  ggtitle("Top 50 genes cluster 3")+
  geom_boxplot()

g<-marker_log_tab%>% dplyr::filter(gene %in% interested_genes4$GeneName) %>%
  group_by(group, cell)%>% summarize(mean_exp=mean(lgFC)) %>%
  ggplot(aes(y=mean_exp, x=group, colour=group))+
  ggtitle("Top 50 genes cluster 4")+
  geom_boxplot()
gridExtra::grid.arrange(d, e, f, g, nrow=2)

marker_log_tab %>% dplyr::filter(gene %in% interested_genes1$GeneName[1:10]) %>%
  group_by(group, gene)%>% summarize(mean_exp=mean(lgFC)) %>%
  ggplot(aes(y=mean_exp, x=group, colour=group))+
  geom_bar(stat="identity")+
  facet_wrap(~gene, ncol=4)


marker_log_tab %>% dplyr::filter(gene %in% interested_genes1$GeneName[1:20]) %>%
  ggplot(aes(y=lgFC, x=group, colour=group))+
  geom_boxplot()+
  facet_wrap(~gene, ncol=5)

  
### plot 20 plots each time ###
start<-seq(1, 188, 20)
end<-c(seq(20, 188, 20), 188)
plot_density<-function(i)
{
  g<-marker_SC_tab %>% dplyr::filter(gene %in% marker_genes[start[i]:end[i]]) %>%
    ggplot(aes(x=SC, group=group, colour=group))+
    geom_density(position="identity", fill="NA")+
    facet_wrap(~gene, ncol=4)
  text_df<-marker_SC_tab %>% dplyr::filter(gene %in% marker_genes[start[i]:end[i]]) %>% group_by(gene, group) %>% summarize(label=as.character(n()))
  text_df$y<-0
  text_df$y[text_df$group==c("Apop")]<-0
  text_df$y[text_df$group==c("Imat")]<-0.05
  text_df$y[text_df$group==c("Cont")]<-0.1
  g+ geom_text(
    data=text_df,
    mapping=aes(x = -Inf, y = y, label = label),
    hjust   = -0.1,
    vjust   = -1
  )
}

plot_density(1)
#### look at the SC seperations from 10x cluster genes ###
d<-marker_SC_tab %>% dplyr::filter(gene %in% interested_genes1$GeneName) %>%
  group_by(group, cell) %>% summarise(mean_SC=mean(SC)) %>%
  ggplot(aes(x=group, y=mean_SC, colour=group))+
  ggtitle("Top 50 genes cluster 1")+
  geom_boxplot()

e<-marker_SC_tab %>% dplyr::filter(gene %in% interested_genes2$GeneName) %>%
  group_by(group, cell) %>% summarise(mean_SC=mean(SC)) %>%
  ggplot(aes(x=group, y=mean_SC, colour=group))+
  ggtitle("Top 50 genes cluster 2")+
  geom_boxplot()

f<-marker_SC_tab %>% dplyr::filter(gene %in% interested_genes3$GeneName) %>%
  group_by(group, cell) %>% summarise(mean_SC=mean(SC)) %>%
  ggplot(aes(x=group, y=mean_SC, colour=group))+
  ggtitle("Top 50 genes cluster 3")+
  geom_boxplot()

g<-marker_SC_tab %>% dplyr::filter(gene %in% interested_genes4$GeneName) %>%
  group_by(group, cell) %>% summarise(mean_SC=mean(SC)) %>%
  ggplot(aes(x=group, y=mean_SC, colour=group))+
  ggtitle("Top 50 genes cluster 4")+
  geom_boxplot()

gridExtra::grid.arrange(d,e,f,g, nrow=2)

###build functions for density cluster 1, 2, 3, 4 ###
### first look at log2(normalised count) ###
clust1_apop_df<-approxfun(density(apply(as.matrix(apsoc@data[interested_genes1$GeneName,]), 2, mean)))
clust1_imat_df<-approxfun(density(apply(as.matrix(rawImat@data[interested_genes1$GeneName,]), 2, mean)))
clust1_cont_df<-approxfun(density(apply(as.matrix(rawCont@data[interested_genes1$GeneName,]), 2, mean)))
clust1_apop_LN<-clust1_apop_df(apply(as.matrix(Aggr@data[interested_genes1$GeneName,]), 2, mean))
clust1_imat_LN<-clust1_imat_df(apply(as.matrix(Aggr@data[interested_genes1$GeneName,]), 2, mean))
clust1_cont_LN<-clust1_cont_df(apply(as.matrix(Aggr@data[interested_genes1$GeneName,]), 2, mean))
clust1_apop_LN[is.na(clust1_apop_LN)]<-0
clust1_imat_LN[is.na(clust1_imat_LN)]<-0
clust1_cont_LN[is.na(clust1_cont_LN)]<-0
clust1_ident<-c()
i=1
for (i in 1:length(clust1_apop_LN)){
  Apop_LN<-clust1_apop_LN[i]
  Imat_LN<-clust1_imat_LN[i]
  Cont_LN<-clust1_cont_LN[i]
  Ident<-c("undetermined")
  if (Apop_LN > Imat_LN & Apop_LN >Cont_LN){Ident<-c("Apop")}
  if (Cont_LN > Apop_LN & Cont_LN > Imat_LN){Ident<-c("Cont")}
  if (Imat_LN >Apop_LN & Imat_LN > Cont_LN){Ident<-c("Imat")}
  clust1_ident<-c(clust1_ident, Ident)
}
clust1_ident_table<-data.frame(cell=colnames(Aggr@data), ident=clust1_ident)
clust1_ident_table<- clust1_ident_table %>% separate(cell, c("cell", "index"))
compare_clust1<-right_join(SC_tab[, c("cell", "group")], clust1_ident_table, by=c("cell"))
table(compare_clust1$group)
table(compare_clust1$group[compare_clust1$group==compare_clust1$ident])

### cluster 2 ###
clust2_apop_df<-approxfun(density(apply(as.matrix(apsoc@data[interested_genes2$GeneName,]), 2, mean)))
clust2_imat_df<-approxfun(density(apply(as.matrix(rawImat@data[interested_genes2$GeneName,]), 2, mean)))
clust2_cont_df<-approxfun(density(apply(as.matrix(rawCont@data[interested_genes2$GeneName,]), 2, mean)))
clust2_apop_LN<-clust2_apop_df(apply(as.matrix(Aggr@data[interested_genes2$GeneName,]), 2, mean))
clust2_imat_LN<-clust2_imat_df(apply(as.matrix(Aggr@data[interested_genes2$GeneName,]), 2, mean))
clust2_cont_LN<-clust2_cont_df(apply(as.matrix(Aggr@data[interested_genes2$GeneName,]), 2, mean))
clust2_apop_LN[is.na(clust2_apop_LN)]<-0
clust2_imat_LN[is.na(clust2_imat_LN)]<-0
clust2_cont_LN[is.na(clust2_cont_LN)]<-0
clust2_ident<-c()
i=1
for (i in 1:length(clust2_apop_LN)){
  Apop_LN<-clust2_apop_LN[i]
  Imat_LN<-clust2_imat_LN[i]
  Cont_LN<-clust2_cont_LN[i]
  Ident<-c("undetermined")
  if (Apop_LN > Imat_LN & Apop_LN >Cont_LN){Ident<-c("Apop")}
  if (Cont_LN > Apop_LN & Cont_LN > Imat_LN){Ident<-c("Cont")}
  if (Imat_LN >Apop_LN & Imat_LN > Cont_LN){Ident<-c("Imat")}
  clust2_ident<-c(clust2_ident, Ident)
}
clust2_ident_table<-data.frame(cell=colnames(Aggr@data), ident=clust2_ident)
clust2_ident_table<- clust2_ident_table %>% separate(cell, c("cell", "index"))
compare_clust2<-right_join(SC_tab[, c("cell", "group")], clust2_ident_table, by=c("cell"))
table(compare_clust2$group)
table(compare_clust2$group[compare_clust2$group==compare_clust2$ident])

### cluster 3 ###
clust3_apop_df<-approxfun(density(apply(as.matrix(apsoc@data[interested_genes3$GeneName,]), 2, mean)))
clust3_imat_df<-approxfun(density(apply(as.matrix(rawImat@data[interested_genes3$GeneName,]), 2, mean)))
clust3_cont_df<-approxfun(density(apply(as.matrix(rawCont@data[interested_genes3$GeneName,]), 2, mean)))
clust3_apop_LN<-clust3_apop_df(apply(as.matrix(Aggr@data[interested_genes3$GeneName,]), 2, mean))
clust3_imat_LN<-clust3_imat_df(apply(as.matrix(Aggr@data[interested_genes3$GeneName,]), 2, mean))
clust3_cont_LN<-clust3_cont_df(apply(as.matrix(Aggr@data[interested_genes3$GeneName,]), 2, mean))
clust3_apop_LN[is.na(clust3_apop_LN)]<-0
clust3_imat_LN[is.na(clust3_imat_LN)]<-0
clust3_cont_LN[is.na(clust3_cont_LN)]<-0
clust3_ident<-c()
i=1
for (i in 1:length(clust3_apop_LN)){
  Apop_LN<-clust3_apop_LN[i]
  Imat_LN<-clust3_imat_LN[i]
  Cont_LN<-clust3_cont_LN[i]
  Ident<-c("undetermined")
  if (Apop_LN > Imat_LN & Apop_LN >Cont_LN){Ident<-c("Apop")}
  if (Cont_LN > Apop_LN & Cont_LN > Imat_LN){Ident<-c("Cont")}
  if (Imat_LN >Apop_LN & Imat_LN > Cont_LN){Ident<-c("Imat")}
  clust3_ident<-c(clust3_ident, Ident)
}
clust3_ident_table<-data.frame(cell=colnames(Aggr@data), ident=clust3_ident)
clust3_ident_table<- clust3_ident_table %>% separate(cell, c("cell", "index"))
compare_clust3<-right_join(SC_tab[, c("cell", "group")], clust3_ident_table, by=c("cell"))
table(compare_clust3$group)
table(compare_clust3$group[compare_clust3$group==compare_clust3$ident])

### cluster 4 ###
clust4_apop_df<-approxfun(density(apply(as.matrix(apsoc@data[interested_genes4$GeneName,]), 2, mean)))
clust4_imat_df<-approxfun(density(apply(as.matrix(rawImat@data[interested_genes4$GeneName,]), 2, mean)))
clust4_cont_df<-approxfun(density(apply(as.matrix(rawCont@data[interested_genes4$GeneName,]), 2, mean)))
clust4_apop_LN<-clust4_apop_df(apply(as.matrix(Aggr@data[interested_genes4$GeneName,]), 2, mean))
clust4_imat_LN<-clust4_imat_df(apply(as.matrix(Aggr@data[interested_genes4$GeneName,]), 2, mean))
clust4_cont_LN<-clust4_cont_df(apply(as.matrix(Aggr@data[interested_genes4$GeneName,]), 2, mean))
clust4_apop_LN[is.na(clust4_apop_LN)]<-0
clust4_imat_LN[is.na(clust4_imat_LN)]<-0
clust4_cont_LN[is.na(clust4_cont_LN)]<-0
clust4_ident<-c()
i=1
for (i in 1:length(clust4_apop_LN)){
  Apop_LN<-clust4_apop_LN[i]
  Imat_LN<-clust4_imat_LN[i]
  Cont_LN<-clust4_cont_LN[i]
  Ident<-c("undetermined")
  if (Apop_LN > Imat_LN & Apop_LN >Cont_LN){Ident<-c("Apop")}
  if (Cont_LN > Apop_LN & Cont_LN > Imat_LN){Ident<-c("Cont")}
  if (Imat_LN >Apop_LN & Imat_LN > Cont_LN){Ident<-c("Imat")}
  clust4_ident<-c(clust4_ident, Ident)
}
clust4_ident_table<-data.frame(cell=colnames(Aggr@data), ident=clust4_ident)
clust4_ident_table<- clust4_ident_table %>% separate(cell, c("cell", "index"))
compare_clust4<-right_join(SC_tab[, c("cell", "group")], clust4_ident_table, by=c("cell"))
table(compare_clust4$group)
table(compare_clust4$group[compare_clust4$group==compare_clust4$ident])


#### not great because a few density curve is over 1, need a better way to compute likelihood ###
library(pracma);trapz(dens$x[dens$x < 2], dens$y[dens$x < 2])
