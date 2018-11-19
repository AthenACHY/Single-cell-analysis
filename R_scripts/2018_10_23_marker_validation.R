### only choose the best genes ###
Apop_marker<-c("PLCG2", "MT-ND5", "MT-ND3")
imat1_marker<-c("PKM", "ACTG1", "VIM")
imat2_marker<-c("ALAS2", "HBZ")
cont_marker<-c("NME1", "HIST1H4C")
HK_marker<-c("RPL8" , "RPS13" , "C9orf16" , "RPS6" , "COX8A" , "FTH1" , "PFDN5" , "ATP5G2" , "NACA" , "ERP29" , "HNRNPA1" , "POMP" , "TPT1" , "RPL6" , "RPLP0" , "SRP14" , "B2M" , "KIAA0101" , "RPL4" , "RPLP1" , "LSM7" , "RPS13" , "RPL7A" , "RPLP0" , "RPL38" , "GAPDH" , "RPL23")
### check for dispersion ###
Mag_genes_dispersion[Apop_marker,]
Mag_genes_dispersion[imat1_marker,]
Mag_genes_dispersion[imat2_marker,]
Mag_genes_dispersion[cont_marker,]
Mag_genes_dispersion[HK_marker,]

VlnPlot(Mag.combined, HK_marker[1:10])

### now calculate singe-cell gene signature score for each cell ###
### extract scaled-expression for each gene ###
averger_SC_cont<-apply(Mag.combined@scale.data[HK_marker, ], 2, mean)
Apop_relative<-apply(Mag.combined@scale.data[Apop_marker, ], 2, mean) - averger_SC_cont
Imat1_relative<-apply(Mag.combined@scale.data[imat1_marker, ], 2, mean) - averger_SC_cont
Imat2_relative<-apply(Mag.combined@scale.data[imat2_marker, ], 2, mean) - averger_SC_cont

Cont_relative<-apply(Mag.combined@scale.data[cont_marker, ], 2, mean) - averger_SC_cont

### plot density in synchronixe way ###
plot_tab<-data.frame(samples=attributes(averger_SC_cont), exp=averger_SC_cont) 
plot_tab<-plot_tab %>% mutate(group=substring(names, 1, 4))
plot_tab %>% ggplot(aes(x=exp, group=group)) +
  geom_density()+
  facet_grid(group~.)
### HKgenes distribution looks awful !!! ###
plot_tab<-data.frame(samples=attributes(Apop_relative), ident=Mag.combined@ident, exp=Apop_relative) 
plot_tab<-plot_tab %>% dplyr::filter(ident %in% c(0, 1, 2, 3, 4)) %>%
  mutate(group=substring(names, 1, 4))
plot_tab %>% ggplot(aes(x=exp,  fill=ident)) +
  geom_density()+
  facet_grid(ident~.)
### 
plot_tab<-data.frame(samples=attributes(Imat1_relative), exp=Imat1_relative) 
plot_tab<-plot_tab %>% mutate(group=substring(names, 1, 4))
plot_tab %>% ggplot(aes(x=exp, group=group)) +
  geom_density()+
  facet_grid(group~.)
###
plot_tab<-data.frame(samples=attributes(Imat2_relative), exp=Imat2_relative) 
plot_tab<-plot_tab %>% mutate(group=substring(names, 1, 4))
plot_tab %>% ggplot(aes(x=exp, group=group)) +
  geom_density()+
  facet_grid(group~.)
###
plot_tab<-data.frame(samples=attributes(Cont_relative), exp=Cont_relative) 
plot_tab<-plot_tab %>% mutate(group=substring(names, 1, 4))
plot_tab %>% ggplot(aes(x=exp, group=group)) +
  geom_density()+
  facet_grid(group~.)

### plot overall density ###
bind_rows(data.frame(samples=attributes(Apop_relative), exp=Apop_relative , ident=Mag.combined@ident) %>% dplyr::mutate(marker=c("Apop_marker")), 
          data.frame(samples=attributes(Imat1_relative), exp=Imat1_relative, ident=Mag.combined@ident) %>% dplyr::mutate(marker=c("Imat1_marker")), 
          data.frame(samples=attributes(Imat2_relative), exp=Imat2_relative, ident=Mag.combined@ident) %>% dplyr::mutate(marker=c("Imat2_marker")), 
          data.frame(samples=attributes(Cont_relative), exp=Cont_relative, ident=Mag.combined@ident) %>% dplyr::mutate(marker=c("Cont_marker"))) %>% 
  mutate(group=substring(names, 1, 4)) %>%
  dplyr::filter(ident %in% c(0, 1, 2, 3, 4)) %>%
  ggplot(aes(x=exp, group=ident)) +
  geom_density()+
  facet_grid(ident~marker)


### compute overlapping between ident ###
Mag_cells<-colnames(Mag.combined@raw.data)
Mag_cells2<-Mag_cells[Mag.combined@ident==2]
Mag_cells1<-Mag_cells[Mag.combined@ident==1]
Mag_cellsC<-c(Mag_cells[Mag.combined@ident==0], Mag_cells[Mag.combined@ident==3])
Mag_cells4<-Mag_cells[Mag.combined@ident==4]
### create distribution function for Apop ###
A_p<-approxfun(density(Apop_relative[Mag_cells4]))

I1_p<-ecdf(Apop_relative[Mag_cells1])
I2_p<-ecdf(Apop_relative[Mag_cells2])
C_p<-ecdf(Apop_relative[Mag_cellsC])

### make distribution for Imat2 ###
A_p<-ecdf(Imat2_relative[Mag_cells4])
I1_p<-ecdf(Imat2_relative[Mag_cells1])
I2_p<-ecdf(Imat2_relative[Mag_cells2])
C_p<-ecdf(Imat2_relative[Mag_cellsC])

C_p(0)
I2_p(0)
A_p(0)

