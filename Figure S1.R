########
#Fig.S1
########
#clean environmental variables in case of unnecessary troubles
rm(list =ls()) 

#set working path
setwd("C:/Users/xiaoweilong/Desktop/Now!/Code/Figure 2/")
getwd()

#load needed R package
if(!require(Seurat))install.packages("Seurat")
if(!require(qs))install.packages("qs")

if(!require(ggplot2))install.packages("ggplot2")
if(!require(scRNAtoolVis))devtools::install_github('junjunlab/scRNAtoolVis')
#-------------------------------------------------------------------------------
#########
#Fig.S1A
#########
Integration <- qread("Integration.qs")

Integration$orig.ident <- factor(Integration$orig.ident,levels = c("ND1","ND2","ND3",
                                                              "HFD1","HFD2","HFD3",
                                                              "TRF1","TRF2","TRF3"))

VlnPlot(Integration,features = c("nCount_RNA","nFeature_RNA","percent.mt"),
        group.by = "orig.ident",ncol=3,pt.size = 0,log = F)
#-------------------------------------------------------------------------------
aggregate(Integration$nCount_RNA, by = list(type = Integration$orig.ident), FUN = mean)
aggregate(Integration$nCount_RNA, by = list(type = Integration$orig.ident), FUN = sd)

aggregate(Integration$nFeature_RNA, by = list(type = Integration$orig.ident), FUN = mean)
aggregate(Integration$nFeature_RNA, by = list(type = Integration$orig.ident), FUN = sd)

aggregate(Integration$percent.mt, by = list(type = Integration$orig.ident), FUN = mean)
aggregate(Integration$percent.mt, by = list(type = Integration$orig.ident), FUN = sd)
#-------------------------------------------------------------------------------
#########
#Fig.S1B
#########
Integration$Group <- factor(Integration$Group,levels = c("ND","HFD","TRF"))
                                                                   
FeatureScatter(Integration,feature1 = "nFeature_RNA",feature2 = "nCount_RNA",
               pt.size = 1,log = T,raster = T,group.by = "Group")
#-------------------------------------------------------------------------------
#########
#Fig.S1C
#########
p1 <- DimPlot(Integration,group.by = "Group",pt.size = 0.5)

p2 <- DimPlot(Integration,split.by = "orig.ident",ncol = 3,
        group.by = "orig.ident",pt.size = 0.5)

p1|p2
#-------------------------------------------------------------------------------
#########
#Fig.S1D
#########
FeaturePlot(Integration,label.size = 1,features = c("Hnf4a","Scd1","Rsad2","Pkhd1",
                                     "Reln","Thsd4","Stab2","Ptprc",
                                     "Cd5l","Skap1","Ccr2","Irf8"),ncol = 4)
#-------------------------------------------------------------------------------