Apop_marker<-c("PLCG2" , "MALAT1" , "MT-ND5" , "APOE" , "SLC25A37" , "HBZ" , "MT-ND2" , "MT-ND3" , "MT-CYB" , "HBA1" , "MT-CO1" , "MT-ATP6" , "MT-ND4" , "MT-ND1" , "APOC1" , "MT-CO2" , "MT-CO3" , "EID1" , "EIF5B" , "WBP5")
imat_marker<-c("EEF1A1" , "NPC2" , "ACTG1" , "S100A10" , "KRT19" , "CTSL" , "KIAA0101" , "ALDOA" , "ANXA2" , "TUBB2B" , "KRT8" , "ALAS2" , "VIM" , "PKM" , "PRDX1" , "KRT18" , "TMSB10" , "TIMP1")
cont_marker<-c("HIST1H4C" , "SRM" , "NME1" , "PPP1R14B" , "GCSH" , "MYC" , "FABP5" , "ODC1" , "CHCHD10" , "C17orf89" , "TMEM147" , "EBNA1BP2" , "BOLA3" , "SAMSN1" , "FAM173A" , "PRMT1" , "PIM1")
HK_marker<-c("ENY2" , "EIF3H" , "RPL8" , "RPS13" , "C9orf16" , "RPS6" , "COX8A" , "FTH1" , "PFDN5" , "ATP5G2" , "NACA" , "ERP29" , "HNRNPA1" , "POMP" , "TPT1" , "RPL6" , "RPLP0" , "SRP14" , "B2M" , "KIAA0101" , "RPL4" , "RPLP1" , "LSM7" , "RPS13" , "RPL7A" , "RPLP0" , "RPL38" , "GAPDH" , "RPL23")

### now calculate singe-cell gene signature score for each cell ###
### extract scaled-expression for each gene ###
averger_SC_cont<-apply(Mag.combined@scale.data[HK_marker, ], 2, mean)
Apop_relative<-apply(Mag.combined@scale.data[Apop_marker, ], 2, mean) - averger_SC_cont
Imat_relative<-apply(Mag.combined@scale.data[imat_marker, ], 2, mean) - averger_SC_cont
Cont_relative<-apply(Mag.combined@scale.data[cont_marker, ], 2, mean) - averger_SC_cont

### check for dispersion ###
Mag_genes_dispersion[Apop_marker,]


par(mfrow=c(3,1))
plot(density((data.frame(averger_SC_cont)[grep("^Cont", row.names(data.frame(averger_SC_cont))), ])), main="Cont")
plot(density((data.frame(averger_SC_cont)[grep("^Apop", row.names(data.frame(averger_SC_cont))), ])), main="Apop")
plot(density((data.frame(averger_SC_cont)[grep("^Imat", row.names(data.frame(averger_SC_cont))), ])), main="Imat")

par(mfrow=c(3,1))
plot(density((data.frame(Apop_relative)[grep("^Cont", row.names(data.frame(Apop_relative))), ])), main="Cont")
plot(density((data.frame(Apop_relative)[grep("^Apop", row.names(data.frame(Apop_relative))), ])), main="Apop")
plot(density((data.frame(Apop_relative)[grep("^Imat", row.names(data.frame(Apop_relative))), ])), main="Imat")
par(mfrow=c(3,1))
plot(density((data.frame(Imat_relative)[grep("^Cont", row.names(data.frame(Imat_relative))), ])), main="Cont")
plot(density((data.frame(Imat_relative)[grep("^Apop", row.names(data.frame(Imat_relative))), ])), main="Apop")
plot(density((data.frame(Imat_relative)[grep("^Imat", row.names(data.frame(Imat_relative))), ])), main="Imat")
par(mfrow=c(3,1))
plot(density((data.frame(Cont_relative)[grep("^Cont", row.names(data.frame(Cont_relative))), ])), main="Cont")
plot(density((data.frame(Cont_relative)[grep("^Apop", row.names(data.frame(Cont_relative))), ])), main="Apop")
plot(density((data.frame(Cont_relative)[grep("^Imat", row.names(data.frame(Cont_relative))), ])), main="Imat")


### plot correlation? ###
data.frame(cells=attributes(Apop_relative), Apop_SC=Apop_relative, Imat_SC=Imat_relative, Cont_Sc=Cont_relative) %>%
  dplyr::mutate()