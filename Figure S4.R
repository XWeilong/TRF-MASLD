########
#Fig.S4
########
#clean environmental variables in case of unnecessary troubles
rm(list =ls()) 
gc()

#set working path
setwd("C:/Users/xiaoweilong/Desktop/Now!/Code/Figure 4/")
getwd()

library(Seurat)
library(qs)
library(ggplot2)
#-------------------------------------------------------------------------------
#########
#Fig.S4A
#########

Integration <- qread("../Figure 3/Zonation.qs")

VlnPlot(Integration,features = c("percent.mt","nCount_RNA"), group.by = "Subtype",
        ncol = 1, pt.size = 0, sort = T) +NoLegend()

aggregate(Integration$percent.mt, by = list(type = Integration$Subtype), FUN = mean)
aggregate(Integration$percent.mt, by = list(type = Integration$Subtype), FUN = sd)

aggregate(Integration$nCount_RNA, by = list(type = Integration$Subtype), FUN = mean)
aggregate(Integration$nCount_RNA, by = list(type = Integration$Subtype), FUN = sd)
#-------------------------------------------------------------------------------
#########
#Fig.S4B
#########

gene=rev(c("Alb","C3","Fga","Fgb","Fgg","Ttr","Trf","Apoa1","Apoe"))

options(repr.plot.width=18, repr.plot.height=10)
p=DotPlot(Integration,features=rev(gene), dot.scale=11, scale.min = 0, scale.max = 100, group.by = "Subtype")+scale_y_discrete("")+scale_x_discrete("")#+coord_flip()
p1=p+theme(panel.grid=element_blank())+theme_bw(base_size = 13)+#RotatedAxis()+
  scale_color_viridis_c()+labs(x=NULL,y=NULL)+theme(aspect.ratio=1.3)+ggtitle("")
plot(p1)
#-------------------------------------------------------------------------------
#########
#Fig.S4C
#########
Genelist <- read.csv("Fig.4F_Plasma protein.csv", header = T)
Genelist <- as.data.frame(Genelist[,1], ncol = 1)
Genelist <- as.list(Genelist)
DefaultAssay(Integration) <- "SCT"
Integration <- AddModuleScore(Integration, assay = "SCT", features = Genelist, name = "Secretion_Score")
Integration@meta.data[,9] <- NULL #Secretion score calculated based on uncorrected SCT.
colnames(Integration@meta.data)[ncol(Integration@meta.data)] <- "Secretion_Score"

Integration$MitoHigh <- "0"
Integration$MitoHigh[Integration$Subtype=="MitoHigh"] <- "MitoHigh"
Integration$MitoHigh[Integration$Subtype!="MitoHigh"] <- "Non-MitoHigh"

p1 <- FeatureScatter(Integration,feature2 = "Secretion_Score",feature1 = "percent.mt", 
               group.by = "MitoHigh")
p1
#-------------------------------------------------------------------------------
#########
#Fig.S4D
#########
MitoHigh <- subset(Integration, Subtype=="MitoHigh")
DefaultAssay(MitoHigh) <- "RNA"
MitoHigh@assays$SCT <- NULL
MitoHigh <- JoinLayers(MitoHigh,assay = "RNA")
MitoHigh[["RNA"]] <- split(MitoHigh[["RNA"]],f=MitoHigh$orig.ident)

MitoHigh <- SCTransform(MitoHigh)
MitoHigh <- RunPCA(MitoHigh)

MitoHigh <- IntegrateLayers(object = MitoHigh, method = HarmonyIntegration,
                            orig.reduction = "pca", new.reduction = 'harmony',
                            assay = "SCT", verbose = TRUE)
#-------------------------------------------------------------------------------
MitoHigh <- FindNeighbors(MitoHigh, reduction = "harmony", dims = 1:30)
MitoHigh <- FindClusters(MitoHigh, resolution = 0.5)
MitoHigh <- RunUMAP(MitoHigh, reduction = "harmony", dims = 1:30)
#-------------------------------------------------------------------------------
DimPlot(MitoHigh, reduction = "umap", label = T, repel = T)+NoLegend()
FeaturePlot(MitoHigh,features = c("Slc1a2","Cyp2e1","Hamp","Cyp2f2"), reduction = "umap", ncol = 2)
DotPlot(MitoHigh,features = c("Slc1a2","Cyp2e1","Hamp","Cyp2f2"))
#-------------------------------------------------------------------------------
MitoHigh$Sublevel <- "Sublevel"
MitoHigh$Sublevel[which(MitoHigh$seurat_clusters==0)] <- "MitoHigh_PP"
MitoHigh$Sublevel[which(MitoHigh$seurat_clusters==1)] <- "MitoHigh_PP"
MitoHigh$Sublevel[which(MitoHigh$seurat_clusters==2)] <- "MitoHigh_MZ"
MitoHigh$Sublevel[which(MitoHigh$seurat_clusters==3)] <- "MitoHigh_PP"
MitoHigh$Sublevel[which(MitoHigh$seurat_clusters==4)] <- "MitoHigh_MZ"
MitoHigh$Sublevel[which(MitoHigh$seurat_clusters==5)] <- "MitoHigh_PC"
#-------------------------------------------------------------------------------
p2 <- FeatureScatter(MitoHigh,group.by = "Sublevel",
      feature2 = "Secretion_Score",feature1 = "Zonation_Score")
p2
#-------------------------------------------------------------------------------
p1|p2
#-------------------------------------------------------------------------------
#########
#Fig.S4E
#########
gene=c("Apoa1","Apoa2","Apoa4","Apoa5",
       "Apob","Apoc1","Apoc2","Apoc3","Apoc4","Apoe")
Genelist <- as.data.frame(gene, ncol = 1)
Genelist <- as.list(Genelist)

DefaultAssay(Integration) <- "SCT"
Integration <- AddModuleScore(Integration, assay = "SCT", features = Genelist, name = "Apolipoproteins")
colnames(Integration@meta.data)[ncol(Integration@meta.data)] <- "Apolipoproteins"
p1 <- VlnPlot(Integration,features = "Apolipoproteins", group.by =  "Subtype", sort = T, pt.size = 0, ncol=1) +NoLegend()

aggregate(Integration$Apolipoproteins, by = list(type = Integration$Subtype), FUN = mean)
aggregate(Integration$Apolipoproteins, by = list(type = Integration$Subtype), FUN = sd)
#-------------------------------------------------------------------------------
#########
#Fig.S4F
#########
#Lipid homeostasis
GeneSet <- read.csv("Fig.S4F_Genesets.csv")
Genelist <- as.data.frame(GeneSet[,3])
Genelist <- as.list(Genelist)
DefaultAssay(Integration) <- "SCT"
Integration <- AddModuleScore(Integration, assay = "SCT", features = Genelist, name = "Lipid homeostasis")
colnames(Integration@meta.data)[ncol(Integration@meta.data)] <- "Lipid homeostasis"
p2 <- VlnPlot(Integration,features = "Lipid homeostasis", group.by =  "Subtype", sort = T, pt.size = 0, ncol=1) +NoLegend()

aggregate(Integration$`Lipid homeostasis`, by = list(type = Integration$Subtype), FUN = mean)
aggregate(Integration$`Lipid homeostasis`, by = list(type = Integration$Subtype), FUN = sd)
#-------------------------------------------------------------------------------
p1|p2
#-------------------------------------------------------------------------------
#Triglyceride homeostasis
GeneSet <- read.csv("Fig.S4F_Genesets.csv")
Genelist <- as.data.frame(GeneSet[,1])
Genelist <- as.list(Genelist)
DefaultAssay(Integration) <- "SCT"
Integration <- AddModuleScore(Integration, assay = "SCT", features = Genelist, name = "TG")
colnames(Integration@meta.data)[ncol(Integration@meta.data)] <- "Triglyceride homeostasis"
VlnPlot(Integration[,Integration$Subtype=="MitoHigh"],features = "Triglyceride homeostasis", group.by =  "Group", sort = T, pt.size = 0, ncol=1) +NoLegend()
#-------------------------------------------------------------------------------
#Cholesterol homeostasis
GeneSet <- read.csv("Fig.S4F_Genesets.csv")
Genelist <- as.data.frame(GeneSet[,2])
Genelist <- as.list(Genelist)
DefaultAssay(Integration) <- "SCT"
Integration <- AddModuleScore(Integration, assay = "SCT", features = Genelist, name = "TC")
colnames(Integration@meta.data)[ncol(Integration@meta.data)] <- "Cholesterol homeostasis"
VlnPlot(Integration[,Integration$Subtype=="MitoHigh"],features = "Cholesterol homeostasis", group.by =  "Group", sort = T, pt.size = 0, ncol=1) +NoLegend()
#-------------------------------------------------------------------------------
library(ggpubr)
p3 <- ggpubr::ggviolin(Integration@meta.data[Integration$Subtype=="MitoHigh",],x="Group", y="Triglyceride homeostasis", width = 0.6,
                 color = "black",
                 fill="Group",
                 palette = "npg",
                 add = 'mean_sd',
                 xlab = F,
                 bxp.errorbar=T,
                 bxp.errorbar.width=0.5,
                 size=0.5,
                 outlier.shape=NA,
                 legend = "right",
                 title = "")
#-------------------------------------------------------------------------------
p4 <- ggpubr::ggviolin(Integration@meta.data[Integration$Subtype=="MitoHigh",],x="Group", y="Cholesterol homeostasis", width = 0.6,
                       color = "black",
                       fill="Group",
                       palette = "npg",
                       add = 'mean_sd',
                       xlab = F,
                       bxp.errorbar=T,
                       bxp.errorbar.width=0.5,
                       size=0.5,
                       outlier.shape=NA,
                       legend = "right",
                       title = "")
#-------------------------------------------------------------------------------
kruskal.test(Integration$`Triglyceride homeostasis`[Integration$Subtype=="MitoHigh"] ~ Integration$Group[Integration$Subtype=="MitoHigh"])
DescTools::DunnTest(Integration$`Triglyceride homeostasis`[Integration$Subtype=="MitoHigh"] ~ Integration$Group[Integration$Subtype=="MitoHigh"])

kruskal.test(Integration$`Cholesterol homeostasis`[Integration$Subtype=="MitoHigh"] ~ Integration$Group[Integration$Subtype=="MitoHigh"])
DescTools::DunnTest(Integration$`Cholesterol homeostasis`[Integration$Subtype=="MitoHigh"] ~ Integration$Group[Integration$Subtype=="MitoHigh"])

aggregate(Integration$`Triglyceride homeostasis`[Integration$Subtype=="MitoHigh"], by = list(type = Integration$Group[Integration$Subtype=="MitoHigh"]), FUN = mean)
aggregate(Integration$`Triglyceride homeostasis`[Integration$Subtype=="MitoHigh"], by = list(type = Integration$Group[Integration$Subtype=="MitoHigh"]), FUN = sd)

aggregate(Integration$`Cholesterol homeostasis`[Integration$Subtype=="MitoHigh"], by = list(type = Integration$Group[Integration$Subtype=="MitoHigh"]), FUN = mean)
aggregate(Integration$`Cholesterol homeostasis`[Integration$Subtype=="MitoHigh"], by = list(type = Integration$Group[Integration$Subtype=="MitoHigh"]), FUN = sd)

p3|p4
#-------------------------------------------------------------------------------