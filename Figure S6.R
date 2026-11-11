#clean environmental variables in case of unnecessary troubles
rm(list =ls()) 
gc()

#set working path
setwd("C:/Users/xiaoweilong/Desktop/Now!/Code/Figure 7/")
getwd()

library("Seurat")
library("qs")
if(!require(scRNAtoolVis))devtools::install_github('junjunlab/scRNAtoolVis')
#-------------------------------------------------------------------------------
###########
#Fig.S6A-B
###########
##Experiment data.
#-------------------------------------------------------------------------------
#########
#Fig.S6C
#########
Integration <- qread("../Figure 3/Zonation.qs")

#Bile acid metabolic process
Genelist <- read.csv("Fig.S6A_Gene_set.csv", header = T)
Genelist <- as.data.frame(Genelist[,1], ncol = 1)
Genelist <- as.list(Genelist)

DefaultAssay(Integration) <- "SCT"
Integration <- AddModuleScore(Integration, assay = "SCT", features = Genelist, name = "BAM")
colnames(Integration@meta.data)[ncol(Integration@meta.data)] <- "Bile acid metabolic process"

cor(Integration$`Bile acid metabolic process`,Integration$Zonation_Score,method = "pearson")
VlnPlot(Integration,features = "Bile acid metabolic process", group.by = "Subtype", sort = T, pt.size = 0) +NoLegend()

p1 <- FeatureScatter(Integration,feature2 = "Bile acid metabolic process",feature1 = "Zonation_Score", raster = T, group.by = "Subtype")
#-------------------------------------------------------------------------------
#########
#Fig.S6D
#########
gene=c("Nr1h4","Nr0b2","Nr1h2","Nr1i2","Rxra","Hnf4a","Nr5a2","Ppara")

options(repr.plot.width=18, repr.plot.height=10)
p=DotPlot(Integration,features=rev(gene), dot.scale=11, scale.min = 0, scale.max = 100, group.by = "Group")+scale_y_discrete("")+scale_x_discrete("")#+coord_flip()
p2=p+theme(panel.grid=element_blank())+theme_bw(base_size = 13)+#RotatedAxis()+
  scale_color_viridis_c()+labs(x=NULL,y=NULL)+theme(aspect.ratio=0.8)+ggtitle("")
plot(p2)
#-------------------------------------------------------------------------------