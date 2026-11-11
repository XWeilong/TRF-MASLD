########
#Fig.S3
########
#clean environmental variables in case of unnecessary troubles
rm(list =ls()) 

#set working path
setwd("C:/Users/xiaoweilong/Desktop/Now!/Code/Figure 3/")
getwd()

#load needed R package
if(!require(Seurat))install.packages("Seurat")
if(!require(qs))install.packages("qs")
#-------------------------------------------------------------------------------
#########
#Fig.S3A
#########

Integration <- qread("../Figure 3/Zonation.qs")

FeaturePlot(Integration,features = c("Hnf4a","Rsad2","percent.mt",
                        "Cyp2e1","Hamp","Cyp2f2"), reduction = "umap", ncol = 3)
#-------------------------------------------------------------------------------
#########
#Fig.S3B
#########

Idents(Integration) <- "Subtype"
Integration$Subtype <- factor(Integration$Subtype, levels = c("IS","MitoHigh",
           "Hep_Pcdh11x+","Hep_Egr1+","Hep_Car3+","Hep_Cidec+","Hep_PP","Hep_MZ","Hep_PC"))

DotPlot(Integration,features = c("Slc1a2","Cps1","Sds","Cidec","Car3","Egr1",
                                 "Pcdh11x","Mup20","Rsad2"),group.by = "Subtype")
#-------------------------------------------------------------------------------
#########
#Fig.S3C
#########

library(clusterProfiler)
library(org.Mm.eg.db)

Idents(Integration) <- "Subtype"
Markers <- FindAllMarkers(Integration, logfc.threshold = 0.5, return.thresh = 0.05, 
                          test.use = "wilcox", min.pct = 0.5, only.pos = T)

Markers <- subset(Markers, Markers$p_val_adj < 0.05)
#-------------------------------------------------------------------------------
Markers$cluster <- factor(Markers$cluster, levels = rev(c("IS","MitoHigh",
       "Hep_Pcdh11x+","Hep_Egr1+","Hep_Car3+","Hep_Cidec+","Hep_PP","Hep_MZ","Hep_PC")))

xx.formula <- compareCluster(gene~cluster, data=Markers,keyType = "SYMBOL",
                             fun='enrichGO', OrgDb='org.Mm.eg.db', ont = "BP")

dotplot(xx.formula,showCategory = 2,font.size = 7)

write.csv(xx.formula@compareClusterResult, "Fig.S3C_GO.csv")
#-------------------------------------------------------------------------------