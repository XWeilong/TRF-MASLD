#clean environmental variables in case of unnecessary troubles
rm(list =ls()) 
gc()

#set working path
setwd("C:/Users/xiaoweilong/Desktop/Now!/Code/Figure 4/")
getwd()

library(Seurat)
library(qs)
library(ggpubr)
#-------------------------------------------------------------------------------
########
#Fig.4A
########
Integration <- qread("../Figure 3/Zonation.qs")
Idents(Integration) <- "Subtype"
Integration <- PrepSCTFindMarkers(Integration)
Markers <- FindAllMarkers(Integration, only.pos = T, logfc.threshold = 0.5, min.pct = 0.5)
Markers <- subset(Markers, p_val_adj<0.05)

#GO biological process enrichment
Markers <- subset(Markers,Markers$cluster=="Hep_Cidec+")

library(clusterProfiler)
library(org.Mm.eg.db)

GO <- enrichGO(gene=Markers$gene, keyType = "SYMBOL", ont = "BP", 'org.Mm.eg.db')
GO <- simplify(GO)

dotplot(GO, showCategory = 15, font.size = 7)

Result <- as.data.frame(GO@result)
write.csv(Result,"Fig.4A_Hep_Cidec+.csv")
#-------------------------------------------------------------------------------
# Check cell groups:
table(Idents(Integration))

# Run COSG:
marker_cosg <- cosg(
  Integration,
  groups = "all",
  assay = "SCT",
  slot = "data",
  mu = 1,
  remove_lowly_expressed=TRUE,
  expressed_pct=0.5,
  n_genes_user = 20)

# Check the marker genes:
head(marker_cosg$names)
head(marker_cosg$scores)
#-------------------------------------------------------------------------------
########
#Fig.4B
########
#lipid localization
gene=c("Osbpl3","Cidec","Cd36","Acsl1","Ppara","Plin2","Abcb4","Abcd3",
       "Ces1d","Lima1","Ces1f","Abcc3","Sec24a","Slc25a20","Lpcat3")

#fatty acid transport and lipogenesis
gene=c("Srebf1","Cidec","Cd36","Fabp1","Fabp2","Fabp3","Fabp4",
       "Acly","Acaca","Fasn","Scd1","Elovl6","Dgat1","Dgat2","Plin2")

#fatty acid β-oxidation
gene=c("Acox1","Acox2","Cpt1","Cpt2","Acsl1","Acadl","Ehhadh","Acaa1a")

options(repr.plot.width=20, repr.plot.height=10)
p=DotPlot(Integration,features=rev(gene), dot.scale=8,group.by = "Subtype")+
                     scale_y_discrete("")+scale_x_discrete("")+coord_flip()
p1=p+theme(panel.grid=element_blank())+theme_bw(base_size = 13)+RotatedAxis()+
  scale_color_viridis_c()+labs(x=NULL,y=NULL)+theme(aspect.ratio=1)+ggtitle("")
plot(p1)
#-------------------------------------------------------------------------------
########
#Fig.4C
########
#lipid localization
Genelist <- as.data.frame(c("Osbpl3","Cidec","Cd36","Acsl1","Ppara","Plin2","Abcb4","Abcd3",
                            "Ces1d","Lima1","Ces1f","Abcc3","Sec24a","Slc25a20","Lpcat3"), ncol = 1)
Genelist <- as.list(Genelist)

DefaultAssay(Integration) <- "SCT"
Integration <- AddModuleScore(Integration, assay = "SCT", features = Genelist, name = "Lipid localization")
colnames(Integration@meta.data)[ncol(Integration@meta.data)] <- "Lipid localization"

VlnPlot(Integration,features = "Lipid localization", group.by = "Subtype", sort = T, pt.size = 0) +NoLegend()

aggregate(Integration$`Lipid localization`, by = list(type = Integration$Subtype), FUN = mean)
aggregate(Integration$`Lipid localization`, by = list(type = Integration$Subtype), FUN = sd)

pairwise.wilcox.test(Integration$`Lipid localization`,Integration$Subtype)
#-------------------------------------------------------------------------------
########
#Fig.4D
########
Markers <- FindAllMarkers(Integration, only.pos = T, logfc.threshold = 0.5, min.pct = 0.5)
Markers <- subset(Markers, p_val_adj<0.05)

Markers <- subset(Markers,Markers$cluster=="Hep_Cidec+")
Markers <- Markers[order(Markers$avg_log2FC, decreasing = T),]

GeneModule <- Markers$gene[1:50]
GeneModule <- as.data.frame(GeneModule, ncol = 1)
GeneModule <- as.list(GeneModule)
#-------------------------------------------------------------------------------
Hep_Cidec <- subset(Integration, Subtype=="Hep_Cidec+")
DefaultAssay(Hep_Cidec) <- "RNA"
Hep_Cidec@assays[["SCT"]] <- NULL

Hep_Cidec <- SCTransform(Hep_Cidec)

Hep_Cidec <- AddModuleScore(Hep_Cidec, assay = "SCT", features = GeneModule, name = "Cell identity score")
colnames(Hep_Cidec@meta.data)[ncol(Hep_Cidec@meta.data)] <- "Cell identity score"

p1 <- ggviolin(Hep_Cidec@meta.data, x="Group", y="Cell identity score", width = 0.6, #x="group"
         color = "black",#轮廓颜色
         fill="Group",#填充 #fill="group"
         palette = "npg",
         add = 'mean_sd',
         xlab = F, #不显示x轴的标签
         bxp.errorbar=T,#显示误差条
         bxp.errorbar.width=0.5, #误差条大小
         size=0.5, #箱型图边线的粗细
         outlier.shape=NA, #不显示outlier
         legend = "right",
         title = "")+NoLegend()
#-------------------------------------------------------------------------------
kruskal.test(Hep_Cidec$`Cell identity score` ~ Hep_Cidec$Group)
DescTools::DunnTest(Hep_Cidec$`Cell identity score` ~ Hep_Cidec$Group)
#-------------------------------------------------------------------------------
aggregate(Hep_Cidec$`Cell identity score`, by = list(type = Hep_Cidec$Group), FUN = mean)
aggregate(Hep_Cidec$`Cell identity score`, by = list(type = Hep_Cidec$Group), FUN = sd)
#-------------------------------------------------------------------------------
########
#Fig.4E
########
meta <- subset(Integration@meta.data,Subtype=="Hep_Cidec+")
p2 <- ggviolin(meta,x="Group", y="Lipid localization", width = 0.6, #x="group"
         color = "black",#轮廓颜色
         fill="Group",#填充 #fill="group"
         palette = "npg",
         add = 'mean_sd',
         xlab = F, #不显示x轴的标签
         bxp.errorbar=T,#显示误差条
         bxp.errorbar.width=0.5, #误差条大小
         size=0.5, #箱型图边线的粗细
         outlier.shape=NA, #不显示outlier
         legend = "right",
         title = "")+NoLegend()
#-------------------------------------------------------------------------------
kruskal.test(meta$`Lipid localization` ~ meta$Group)
DescTools::DunnTest(meta$`Lipid localization` ~ meta$Group)
#-------------------------------------------------------------------------------
aggregate(meta$`Lipid localization`, by = list(type = meta$Group), FUN = mean)
aggregate(meta$`Lipid localization`, by = list(type = meta$Group), FUN = sd)
#-------------------------------------------------------------------------------
p1|p2 #Fig.4D+4E
#-------------------------------------------------------------------------------
########
#Fig.4F
########
Genelist <- read.csv("Fig.4F_Plasma protein.csv", header = T)
Genelist <- as.data.frame(Genelist[,1], ncol = 1)
Genelist <- as.list(Genelist)
DefaultAssay(Integration) <- "SCT"
Integration <- AddModuleScore(Integration, assay = "SCT", features = Genelist, name = "Secretion_Score")
Integration@meta.data[,9] <- NULL #Secretion score calculated based on uncorrected SCT.
colnames(Integration@meta.data)[ncol(Integration@meta.data)] <- "Secretion_Score"

p1 <- VlnPlot(Integration,features = "Secretion_Score", group.by = "Subtype", sort = T, pt.size = 0) +NoLegend()

aggregate(Integration$Secretion_Score, by = list(type = Integration$Subtype), FUN = mean)
aggregate(Integration$Secretion_Score, by = list(type = Integration$Subtype), FUN = sd)
#-------------------------------------------------------------------------------
########
#Fig.4G
########
Genelist <- read.csv("Fig.4G_ATP synthesis.csv", header = T)
Genelist <- as.data.frame(Genelist[,1], ncol = 1)
Genelist <- as.list(Genelist)
DefaultAssay(Integration) <- "SCT"
Integration <- AddModuleScore(Integration, assay = "SCT", features = Genelist, name = "ATP biosynthesis")
colnames(Integration@meta.data)[ncol(Integration@meta.data)] <- "ATP biosynthesis"

p2 <- VlnPlot(Integration,features = "ATP biosynthesis", group.by = "Subtype", sort = T, pt.size = 0) +NoLegend()

aggregate(Integration$`ATP biosynthesis`, by = list(type = Integration$Subtype), FUN = mean)
aggregate(Integration$`ATP biosynthesis`, by = list(type = Integration$Subtype), FUN = sd)
#-------------------------------------------------------------------------------
p1|p2 #Fig.4F-G
#-------------------------------------------------------------------------------
########
#Fig.4H
########
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
meta <- subset(MitoHigh@meta.data)

meta$Sublevel <- factor(meta$Sublevel,levels = c("MitoHigh_PC","MitoHigh_MZ","MitoHigh_PP"))

p1 <- ggviolin(meta,x="Sublevel", y="ATP biosynthesis", width = 0.6, #x="group"
               color = "black",#轮廓颜色
               fill="Sublevel",#填充 #fill="group"
               palette = "npg",
               add = 'mean_sd',
               xlab = F, #不显示x轴的标签
               bxp.errorbar=T,#显示误差条
               bxp.errorbar.width=0.5, #误差条大小
               size=0.5, #箱型图边线的粗细
               outlier.shape=NA, #不显示outlier
               #legend = "right",
               title = "")+NoLegend()
#-------------------------------------------------------------------------------
kruskal.test(meta$`ATP biosynthesis` ~ meta$Sublevel)
DescTools::DunnTest(meta$`ATP biosynthesis` ~ meta$Sublevel)
#-------------------------------------------------------------------------------
aggregate(meta$`ATP biosynthesis`, by = list(type = meta$Sublevel), FUN = mean)
aggregate(meta$`ATP biosynthesis`, by = list(type = meta$Sublevel), FUN = sd)
#-------------------------------------------------------------------------------
########
#Fig.4I
########
myCor=cor.test(meta$Secretion_Score,meta$`ATP biosynthesis`,
               alternative = c("two.sided"),method="spearman")

myCor$estimate
myCor$p.value

p2 <- ggplot(meta,aes(y=`ATP biosynthesis`,x=Secretion_Score))+
  geom_smooth(method="lm",  # 拟合曲线
              se=T,         # 是否添加置信区间
              color="red")+
  geom_point()+             # 散点图
  annotate(
    "text", x = min(meta$Secretion_Score,na.rm = T), y = max(meta$ATP,na.rm = T), # 坐标位置，左上角
    label = paste0("R^2=",round(myCor$estimate,2),                                        # 注释标记内容相关性+P值
                   ",P=",round(myCor$p.value,3)),   # P值采用科学计数法，保留2位小数
    size = 5,hjust=0,vjust=1)+
  theme_classic()

write.csv(meta,"Fig.4I_Correlation.csv")
#-------------------------------------------------------------------------------
########
#Fig.4J
########
Markers <- FindAllMarkers(Integration, only.pos = T, logfc.threshold = 0.5, min.pct = 0.5)
Markers <- subset(Markers, p_val_adj<0.05)

Markers <- subset(Markers,Markers$cluster=="MitoHigh")
Markers <- Markers[order(Markers$avg_log2FC, decreasing = T),]

GeneModule <- Markers$gene[1:50]
GeneModule <- as.data.frame(GeneModule, ncol = 1)
GeneModule <- as.list(GeneModule)
#-------------------------------------------------------------------------------
MitoHigh <- subset(Integration, Subtype=="MitoHigh")
DefaultAssay(MitoHigh) <- "RNA"
MitoHigh@assays[["SCT"]] <- NULL
MitoHigh <- SCTransform(MitoHigh)

MitoHigh <- AddModuleScore(MitoHigh, assay = "SCT", features = GeneModule, name = "Cell identity score")
colnames(MitoHigh@meta.data)[ncol(MitoHigh@meta.data)] <- "Cell identity score"

p3 <- ggviolin(MitoHigh@meta.data, x="Group", y="Cell identity score", width = 0.6, #x="group"
               color = "black",#轮廓颜色
               fill="Group",#填充 #fill="group"
               palette = "npg",
               add = 'mean_sd',
               xlab = F, #不显示x轴的标签
               bxp.errorbar=T,#显示误差条
               bxp.errorbar.width=0.5, #误差条大小
               size=0.5, #箱型图边线的粗细
               outlier.shape=NA, #不显示outlier
               #legend = "right",
               title = "")+NoLegend()
#-------------------------------------------------------------------------------
kruskal.test(MitoHigh$`Cell identity score` ~ MitoHigh$Group)
DescTools::DunnTest(MitoHigh$`Cell identity score` ~ MitoHigh$Group)
#-------------------------------------------------------------------------------
aggregate(MitoHigh$`Cell identity score`, by = list(type = MitoHigh$Group), FUN = mean)
aggregate(MitoHigh$`Cell identity score`, by = list(type = MitoHigh$Group), FUN = sd)
#-------------------------------------------------------------------------------
########
#Fig.4K
########
meta <- subset(Integration@meta.data,Subtype=="MitoHigh")

p4 <- ggviolin(meta,x="Group", y="Secretion_Score", width = 0.6, #x="group"
               color = "black",#轮廓颜色
               fill="Group",#填充 #fill="group"
               palette = "npg",
               add = 'mean_sd',
               xlab = F, #不显示x轴的标签
               bxp.errorbar=T,#显示误差条
               bxp.errorbar.width=0.5, #误差条大小
               size=0.5, #箱型图边线的粗细
               outlier.shape=NA, #不显示outlier
               #legend = "right",
               title = "")+NoLegend()
#-------------------------------------------------------------------------------
kruskal.test(meta$Secretion_Score ~ meta$Group)
DescTools::DunnTest(meta$Secretion_Score ~ meta$Group)
#-------------------------------------------------------------------------------
aggregate(meta$Secretion_Score, by = list(type = meta$Group), FUN = mean)
aggregate(meta$Secretion_Score, by = list(type = meta$Group), FUN = sd)
#-------------------------------------------------------------------------------
p1|p2|p3|p4 #Fig4H-K