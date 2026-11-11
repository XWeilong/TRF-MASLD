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
########
#Fig.7A
########
##Images from TEM.
#-------------------------------------------------------------------------------
########
#Fig.7B
########
Integration <- qread("../Figure 3/Zonation.qs")

Idents(Integration) <- "Group"
DEGs <- FindMarkers(Integration, ident.1 = "TRF", ident.2 = "HFD",
                            min.pct = 0.25, logfc.threshold = 0.5)

DEGs <- subset(DEGs, p_val_adj < 0.05)
DEGs$Direction[DEGs$avg_log2FC>0] <- "up"
DEGs$Direction[DEGs$avg_log2FC<0] <- "down"
DEGs$gene <- rownames(DEGs)

write.csv(DEGs,"Fig.7B_DEGs_Hep.csv")
#-------------------------------------------------------------------------------
##Pathway enrichment in Metascape.
#-------------------------------------------------------------------------------
########
#Fig.7C
########
#Bile acid biosynthetic process
Genelist <- read.csv("Fig.7C-D_Gene_sets.csv", header = T)
Genelist <- as.data.frame(Genelist[,1], ncol = 1)
Genelist <- as.list(Genelist)

DefaultAssay(Integration) <- "SCT"
Integration <- AddModuleScore(Integration, features = Genelist, name = "BAB")
colnames(Integration@meta.data)[ncol(Integration@meta.data)] <- "Bile acid biosynthetic process"
VlnPlot(Integration,features = "Bile acid biosynthetic process",group.by = "Group", sort = T, pt.size = 0) +NoLegend()
FeatureScatter(Integration,feature1 = "Zonation_Score",feature2 = "Bile acid biosynthetic process", group.by = "Subtype")
#-------------------------------------------------------------------------------
########
#Fig.7D
########
#Bile acid and bile salt transport
Genelist <- read.csv("Fig.7C-D_Gene_sets.csv", header = T)
Genelist <- as.data.frame(Genelist[,2], ncol = 1)
Genelist <- as.list(Genelist)

DefaultAssay(Integration) <- "SCT"
Integration <- AddModuleScore(Integration, features = Genelist, name = "BAT")
colnames(Integration@meta.data)[ncol(Integration@meta.data)] <- "Bile acid and bile salt transport"
VlnPlot(Integration,features = "Bile acid and bile salt transport",group.by = "Group", sort = T, pt.size = 0) +NoLegend()
FeatureScatter(Integration,feature1 = "Zonation_Score",feature2 = "Bile acid and bile salt transport", group.by = "Subtype")
#-------------------------------------------------------------------------------
library(ggpubr)

p1 <- ggviolin(Integration@meta.data, x="Group", y="Bile acid biosynthetic process", width = 0.6, #x="group"
               color = "black",#轮廓颜色
               fill="Group",#填充 #fill="group"
               palette = "npg",
               add = 'mean_sd',
               xlab = F, #不显示x轴的标签
               bxp.errorbar=T,#显示误差条
               bxp.errorbar.width=0.5, #误差条大小
               size=1, #箱型图边线的粗细
               outlier.shape=NA, #不显示outlier
               legend = "right")+NoLegend()

p2 <- ggviolin(Integration@meta.data, x="Group", y="Bile acid and bile salt transport", width = 0.6, #x="group"
               color = "black",#轮廓颜色
               fill="Group",#填充 #fill="group"
               palette = "npg",
               add = 'mean_sd',
               xlab = F, #不显示x轴的标签
               bxp.errorbar=T,#显示误差条
               bxp.errorbar.width=0.5, #误差条大小
               size=1, #箱型图边线的粗细
               outlier.shape=NA, #不显示outlier
               legend = "right")+NoLegend()
p1|p2
#-------------------------------------------------------------------------------
kruskal.test(Integration$`Bile acid biosynthetic process` ~ Integration$Group)
DescTools::DunnTest(Integration$`Bile acid biosynthetic process` ~ Integration$Group)

aggregate(Integration$`Bile acid biosynthetic process`, by = list(type = Integration$Group), FUN = mean)
aggregate(Integration$`Bile acid biosynthetic process`, by = list(type = Integration$Group), FUN = sd)
#-------------------------------------------------------------------------------
kruskal.test(Integration$`Bile acid and bile salt transport` ~ Integration$Group)
DescTools::DunnTest(Integration$`Bile acid and bile salt transport` ~ Integration$Group)

aggregate(Integration$`Bile acid and bile salt transport`, by = list(type = Integration$Group), FUN = mean)
aggregate(Integration$`Bile acid and bile salt transport`, by = list(type = Integration$Group), FUN = sd)
#-------------------------------------------------------------------------------
########
#Fig.7E
########
#Synthesis
gene=c("Cyp7a1","Cyp8b1","Cyp27a1","Cyp7b1","Cyp2c70")

options(repr.plot.width=18, repr.plot.height=10)
p=DotPlot(Integration,features=gene, dot.scale=11, scale.min = 20, scale.max = 80, group.by = "Group")+scale_y_discrete("")+scale_x_discrete("")#+coord_flip()
p1=p+theme(panel.grid=element_blank())+theme_bw(base_size = 13)+#RotatedAxis()+
  scale_color_viridis_c()+labs(x=NULL,y=NULL)+theme(aspect.ratio=0.3)+ggtitle("")+NoLegend()
plot(p1)
#-------------------------------------------------------------------------------
########
#Fig.7F
########
#Transport
gene=c("Abcb4","Abcc2","Abcc3","Abcb11","Slc10a1","Slco1a1")

options(repr.plot.width=18, repr.plot.height=10)
p=DotPlot(Integration,features=rev(gene), dot.scale=11, scale.min = 20, scale.max = 80, group.by = "Group")+scale_y_discrete("")+scale_x_discrete("")#+coord_flip()
p2=p+theme(panel.grid=element_blank())+theme_bw(base_size = 13)+#RotatedAxis()+
  scale_color_viridis_c()+labs(x=NULL,y=NULL)+theme(aspect.ratio=0.3)+ggtitle("")+NoLegend()
plot(p2)
#-------------------------------------------------------------------------------
cowplot::plot_grid(p1,p2,ncol = 1)
#-------------------------------------------------------------------------------
########
#Fig.7G
########
##Graphpad Prism 10.
#-------------------------------------------------------------------------------