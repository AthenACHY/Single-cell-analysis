### find cluster enriched genes ###
### 0+3 vs the rest 
test_M<-FindMarkers(object = Mag.combined, ident.1 = c(0,3), ident.2=c(1, 2, 4, 5, 6, 7),
                    min.pct = 0.2)
test_M %>% dplyr::mutate(gene=row.names(test_M)) %>%
  dplyr::filter(avg_logFC>0.9) 
### no cluster markers in var.genes ###
VlnPlot(Mag.combined, c("NME1", "HIST1H4C"))
### 1 aginst the rest ###
test_M<-FindMarkers(object = Mag.combined, ident.1 = c(1), ident.2=c(0, 2, 3, 4, 5, 6, 7),
                    min.pct = 0.7)
VlnPlot(Mag.combined, c("ACTG1", "PKM", "VIM"))

### 2 aginst the rest ###
test_M<-FindMarkers(object = Mag.combined, ident.1 = c(2), ident.2=c(0, 1, 3, 4, 5, 6, 7),
                    min.pct = 0.7)
VlnPlot(Mag.combined, c("ALAS2", "HBZ"))
### no clear separation ###
### 4 aginst the rest ###
test_M<-FindMarkers(object = Mag.combined, ident.1 = c(4), ident.2=c(0, 1, 3, 2, 5, 6, 7),
                    min.pct = 0.7)
VlnPlot(Mag.combined, c("PLCG2", "MT-ND5", "MT-ND3"))

