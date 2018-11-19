###randomize 400 cells from Aggregated samples to reveal their identity ####
### use raw count to calculate SC score???###
HK_gene<-c("ATP5E","COX7C","EEF1D","FAU","RPL10A","RPL24","RPL30","RPL35A","RPS10","RPS20", "RPS21","RPS29","SERF2","TOMM7","UQCRB")
Apop_marker<-c("PLCG2", "MT-ND5", "MT-ND3")
imat1_marker<-c("PKM", "ACTG1", "VIM")
imat2_marker<-c("ALAS2", "HBZ")
cont_marker<-c("NME1", "HIST1H4C")

d10x<-Read10X_h5("/Users/athchu/Documents/apop.sc/ApopSc180910/apopAggr/outs/raw_gene_bc_matrices_h5.h5")
Aggr<-CreateSeuratObject(raw.data = d10x, min.cells = 3, min.genes = 200, 
                            project = "Aggr")

Aggr<-NormalizeData(Aggr, scale.factor = 10000)
### check mean HK_genes diff ###
bind_rows(data.frame(SC=log2(apply(apsoc@raw.data[HK_gene, ], 2, mean))) %>% dplyr::mutate(group="Apop"),
          data.frame(SC=log2(apply(rawImat@raw.data[HK_gene, ], 2, mean))) %>% dplyr::mutate(group="Imat"),
          data.frame(SC=log2(apply(rawCont@raw.data[HK_gene, ], 2, mean)))%>% dplyr::mutate(group="Cont")) %>%
  ggplot(aes(x=SC, group=group, colour=group)) + geom_density()



### try to caclulate SC scroe from raw data ###
Apop_SC_apop<-log2(apply(apsoc@raw.data[Apop_marker, ], 2, mean)/apply(apsoc@raw.data[HK_gene, ], 2, mean))
Apop_SC_Imat<-log2(apply(rawImat@raw.data[Apop_marker, ], 2, mean)/apply(rawImat@raw.data[HK_gene, ], 2, mean))
Apop_SC_Cont<-log2(apply(rawCont@raw.data[Apop_marker, ], 2, mean)/apply(rawCont@raw.data[HK_gene, ], 2, mean))


sc_tab<-bind_rows(data.frame(SC=Apop_SC_apop) %>% dplyr::mutate(group=c("Apop")),
          data.frame(SC=Apop_SC_Imat) %>% dplyr::mutate(group=c("Imat")),
          data.frame(SC=Apop_SC_Cont) %>% dplyr::mutate(group=c("cont")))

sc_tab %>% dplyr::filter(SC<100) %>%
          ggplot(aes(x=SC, colour=group, group=group))+
          geom_density(position="identity", fill="NA")
### Imat markers

Imat2_SC_apop<-log2(apply(apsoc@raw.data[imat2_marker, ], 2, mean)/apply(apsoc@raw.data[HK_gene, ], 2, mean))
Imat2_SC_Imat<-log2(apply(rawImat@raw.data[imat2_marker, ], 2, mean)/apply(rawImat@raw.data[HK_gene, ], 2, mean))
Imat2_SC_Cont<-log2(apply(rawCont@raw.data[imat2_marker, ], 2, mean)/apply(rawCont@raw.data[HK_gene, ], 2, mean))


sc_tab<-bind_rows(data.frame(SC=Imat2_SC_apop) %>% dplyr::mutate(group=c("Apop")),
                  data.frame(SC=Imat2_SC_Imat) %>% dplyr::mutate(group=c("Imat")),
                  data.frame(SC=Imat2_SC_Cont) %>% dplyr::mutate(group=c("cont")))

sc_tab %>% dplyr::filter(SC<100) %>%
  ggplot(aes(x=SC, colour=group, group=group))+
  geom_density(position="identity", fill="NA")


### set gates ###
SC_tab<-bind_rows(data.frame(cell=names(Apop_SC_apop), SCapop=Apop_SC_apop, SCImat=Imat2_SC_apop) %>% dplyr::mutate(group=c("Apop")),
          data.frame(cell=names(Apop_SC_Imat), SCapop=Apop_SC_Imat, SCImat=Imat2_SC_Imat) %>% dplyr::mutate(group=c("Imat")),
          data.frame(cell=names(Apop_SC_Cont), SCapop=Apop_SC_Cont, SCImat=Imat2_SC_Cont) %>% dplyr::mutate(group=c("Cont")))

SC_tab %>% ggplot(aes(x=SCapop, y=SCImat, colour=group)) +
  geom_point()

#### calculate SC_raw count from Aggr ###
Imat2_SC_Aggr<-log2(apply(Aggr@raw.data[imat2_marker, ], 2, mean)/apply(Aggr@raw.data[HK_gene, ], 2, mean))
Apop_SC_Aggr<-log2(apply(Aggr@raw.data[Apop_marker, ], 2, mean)/apply(Aggr@raw.data[HK_gene, ], 2, mean))
Aggr_SC_tab<-data.frame(cell=names(Apop_SC_Aggr), SCapop=Apop_SC_Aggr, SCImat=Imat2_SC_Aggr)
###create function ####
Apop_apop_df<-approxfun(density(Apop_SC_apop[!is.na(Apop_SC_apop)]))
Apop_cont_df<-approxfun(density(Apop_SC_Cont[!is.na(Apop_SC_Cont)]))
Imat_apop_df<-approxfun(density(Imat2_SC_apop[!is.na(Imat2_SC_apop)]))
Imat_cont_df<-approxfun(density(Imat2_SC_Cont[!is.na(Imat2_SC_Cont)]))

#### apply to Aggr_SC_tab ###
### cont population has strange behaviour ###
Aggr_SC_tab$Apop_apop_LN<-Apop_apop_df(Aggr_SC_tab$SCapop)
Aggr_SC_tab$Apop_cont_LN<-Apop_cont_df(Aggr_SC_tab$SCapop)
Aggr_SC_tab$Imat_apop_LN<-Imat_apop_df(Aggr_SC_tab$SCImat)
Aggr_SC_tab$Imat_cont_LN<-Imat_cont_df(Aggr_SC_tab$SCImat)

### So compute if there are any agreements ###
Aggr_SC_tab<-Aggr_SC_tab %>% dplyr::mutate(apop_df=ifelse(Apop_apop_LN>Apop_cont_LN, "Apop", "cont"), Imat_df=ifelse(Imat_apop_LN>Imat_cont_LN, "Apop", "cont"))
Aggr_SC_tab<-separate(Aggr_SC_tab, cell, c("cell", "index"), sep=c("-"))
SC_tab<-separate(SC_tab, cell, c("cell", "index"), sep=c("-"))
compare_results<-left_join(SC_tab[, c("cell", "group")], Aggr_SC_tab, by=c("cell"))
compare_results<-compare_results[!is.na(compare_results$SCapop),]
### 7130 guesses in Aggr ###
nrow(compare_results[compare_results$group==compare_results$apop_df,])
nrow(compare_results[compare_results$group==compare_results$Imat_df,])
nrow(compare_results[compare_results$group==compare_results$apop_df & compare_results$group==c("Apop"),])
nrow(compare_results[compare_results$group==compare_results$Imat_df & compare_results$group==c("Apop"),])

