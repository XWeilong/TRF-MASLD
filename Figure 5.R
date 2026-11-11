#clean environmental variables in case of unnecessary troubles
rm(list =ls()) 
gc()

#set working path
setwd("C:/Users/xiaoweilong/Desktop/Now!/Code/Figure 5/")
getwd()

library(Seurat)
library(qs)
#-------------------------------------------------------------------------------
########
#Fig.5A
########
data <- Read10X(data.dir = "Fig.5A_5k_mouse_liver_CNIK_3pv3", gene.column = 2, cell.column = 1, unique.features = T)
data <- CreateSeuratObject(data, min.cells = 1, min.features = 1)
data[["Percent.mt"]] <- PercentageFeatureSet(data, pattern = "^mt-")
data[["Percent.rp"]] <- PercentageFeatureSet(data, pattern = "^Rp[l|s]")

data <- SCTransform(data)
data <- RunPCA(data)
ElbowPlot(data,ndims = 30)
data <- RunUMAP(data,dims = 1:15)
data <- FindNeighbors(data, dims = 1:15)
data <- FindClusters(data, resolution = 0.3)
data$seurat_clusters <- paste0("C",data$seurat_clusters)
DimPlot(data, reduction = "umap", label = T, group.by = "seurat_clusters")

Idents(data) <- "seurat_clusters"
Markers <- FindAllMarkers(data, only.pos = T, logfc.threshold = 1, min.pct = 0.5)
Markers <- subset(Markers, p_val_adj < 0.05)
table(Markers$cluster)

write.csv(Markers,"Fig.5A_Markers.csv")
#-------------------------------------------------------------------------------
#########
#Fig.5B
#########
Reference_10X <- AverageExpression(data, group.by = "seurat_clusters", features = data@assays[["SCT"]]@var.features, assays = "SCT", layer = 'scale.data')
Reference_10X <- as.data.frame(Reference_10X)

Integration <- qread("../Figure 2/Integration.qs")

Query_TRF <- AverageExpression(Integration, group.by = "CellType", features = Integration@assays[["SCT"]]@var.features, assays = "SCT", layer = 'scale.data')
Query_TRF <- as.data.frame(Query_TRF)

overlap <- intersect(rownames(Reference_10X),rownames(Query_TRF))
Reference_10X <- Reference_10X[overlap,]
Query_TRF <- Query_TRF[overlap,]

table(rownames(Reference_10X)==rownames(Query_TRF))

Correlation <- cbind(Reference_10X,Query_TRF)

comp <- cor(Correlation,method = "pearson")

comp <- comp[1:9,10:20]

rownames(comp) <- c("C0","C1","C2","C3","C4","C5","C6","C7","C8")
colnames(comp) <- c("Cho","DC","EC","Hep","HSC","IS","Kup","Lym","MitoHigh","Mono","Myof")

comp <- t(comp)

pheatmap::pheatmap(comp,cluster_rows = T, cluster_cols = T, angle_col = 0)

write.csv(comp,"Fig.5B_Correlation.csv")
#-------------------------------------------------------------------------------
#########
#Fig.5C
#########

VlnPlot(data,features = c("nCount_RNA","nFeature_RNA","Percent.mt","Percent.rp"), pt.size = 0, ncol = 1)

#ModuleScore
Genelist <- read.csv("../Figure 4/Fig.4F_Plasma protein.csv", header = T)
Genelist <- as.data.frame(Genelist[,1], ncol = 1)
Genelist <- as.list(Genelist)

DefaultAssay(data) <- "SCT"
data <- AddModuleScore(data, features = Genelist, assay = "SCT")
colnames(data@meta.data)[length(colnames(data@meta.data))] <- 'Secretion_Score'

VlnPlot(data,features = c("Percent.mt","Secretion_Score"), pt.size = 0, group.by = "seurat_clusters", ncol = 1)

write.csv(data@meta.data,"Fig.5C_Metadata.csv")

aggregate(data$Percent.mt, by = list(type = data$seurat_clusters), FUN = mean)
aggregate(data$Percent.mt, by = list(type = data$seurat_clusters), FUN = sd)

aggregate(data$Secretion_Score, by = list(type = data$seurat_clusters), FUN = mean)
aggregate(data$Secretion_Score, by = list(type = data$seurat_clusters), FUN = sd)
#-------------------------------------------------------------------------------
##########
#Fig.5D/E
##########

rm(list =ls()) 
gc()

Hep_H <- readRDS("Fig.5D_snRNA_Human_GSE223581.rds")
Hep_M <- readRDS("Fig.5E_snRNA_Mouse_GSE223558.rds")

VlnPlot(Hep_H,features = c("nCount_RNA","nFeature_RNA","percent.mito"), group.by = "seurat_clusters", pt.size = 0, ncol = 1)
VlnPlot(Hep_M,features = c("nCount_RNA","nFeature_RNA","percent.mito"), group.by = "seurat_clusters", pt.size = 0, ncol = 1)
#-------------------------------------------------------------------------------
#Fig.5D(left)
FeaturePlot(Hep_H,features = "percent.mito")
Hep_H$seurat_clusters <- paste0("C",Hep_H$seurat_clusters)
DimPlot(Hep_H,group.by = "seurat_clusters",label = T)
#Fig.5E(left)
FeaturePlot(Hep_M,features = "percent.mito")
Hep_M$seurat_clusters <- paste0("C",Hep_M$seurat_clusters)
DimPlot(Hep_M,group.by = "seurat_clusters",label = T)
#-------------------------------------------------------------------------------
#Fig.5D(right)
#Human
Idents(Hep_H) <- "seurat_clusters"
Marker <- FindAllMarkers(Hep_H, only.pos = T, logfc.threshold = 0.5, min.pct = 0.5)
Marker <- subset(Marker,p_val_adj < 0.05)

write.csv(Marker,"Fig.5D_Markers.csv")
#-------------------------------------------------------------------------------
#library(homologene)
#H2M <- human2mouse(rownames(Hep_H))

library(biomaRt)

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

human <- useMart('ensembl',dataset = "hsapiens_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")
mouse <- useMart('ensembl',dataset = "mmusculus_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")

H2M.g <- getLDS(attributes = c("hgnc_symbol"),filters = "hgnc_symbol",
                values = rownames(Hep_H), mart = human,
                attributesL = c("mgi_symbol"), martL = mouse,uniqueRows = T)

#remove duplicates
H2M.g <- H2M.g[!duplicated(H2M.g$HGNC.symbol),]
H2M.g <- H2M.g[!duplicated(H2M.g$MGI.symbol),]
#-------------------------------------------------------------------------------
#Secretion_Score
Genelist <- read.csv("../Figure 4/Fig.4F_Plasma protein.csv", header = T)

#how many human genes in this list
table(Genelist$Liver.Secretion%in%H2M.g$MGI.symbol)

#human liver addmodulescore
Genelist <- H2M.g[H2M.g$MGI.symbol%in%Genelist$Liver.Secretion,]$HGNC.symbol
Genelist <- as.data.frame(Genelist)
Genelist <- as.list(Genelist)

Hep_H <- SCTransform(Hep_H)
Hep_H <- AddModuleScore(Hep_H, features = Genelist, assay = "SCT")
colnames(Hep_H@meta.data)[length(colnames(Hep_H@meta.data))] <- 'Human_Liver_Secretion'
VlnPlot(Hep_H,features = "Human_Liver_Secretion",group.by = "seurat_clusters",pt.size = 0,sort = T)
#-------------------------------------------------------------------------------
##Correlation_Human
sc <- qread("../Figure 3/Zonation.qs")
sc$Heps <- as.character(sc$Subtype)
sc$Heps[sc$Heps!="IS"&sc$Heps!="MitoHigh"] <- "Hep"

AE_Mouse <- AverageExpression(sc, group.by = "Heps", features = sc@assays[["SCT"]]@var.features, assays = "SCT", layer = 'scale.data')
AE_Mouse <- as.data.frame(AE_Mouse)

AE_Human <- AverageExpression(Hep_H, group.by = "seurat_clusters", features = Hep_H@assays[["SCT"]]@var.features, assays = "SCT", layer = 'scale.data')
AE_Human <- as.data.frame(AE_Human)

AE_Human <- AE_Human[rownames(AE_Human)%in%H2M.g$HGNC.symbol,]
rownames(AE_Human) <- H2M.g[H2M.g$HGNC.symbol%in%rownames(AE_Human),]$MGI.symbol

overlap <- intersect(rownames(AE_Human),rownames(AE_Mouse))
AE_Human <- AE_Human[overlap,]
AE_Mouse <- AE_Mouse[overlap,]

table(rownames(AE_Human)==rownames(AE_Mouse))
Correlation <- cbind(AE_Human,AE_Mouse)

#sg <- names(tail(sort(apply(Correlation,1,sd)),3000))

comp <- cor(Correlation,method = "pearson")

comp <- comp[1:9,10:12]

rownames(comp) <- c("C0","C1","C2","C3","C4","C5","C6","C7","C8")
colnames(comp) <- c("Hep","IS","MitoHigh")
pheatmap::pheatmap(comp,cluster_rows = F, cluster_cols = F, angle_col = 0)

write.csv(comp,"Fig.5D_Correlation.csv")
#-------------------------------------------------------------------------------
#Fig.5E(right)
#Mouse
Idents(Hep_M) <- "seurat_clusters"
Marker <- FindAllMarkers(Hep_M, only.pos = T, logfc.threshold = 0.5, min.pct = 0.5)
Marker <- subset(Marker,p_val_adj < 0.05)

write.csv(Marker,"Fig.5E_Markers.csv")
#-------------------------------------------------------------------------------
Genelist <- read.csv("../Figure 4/Fig.4F_Plasma protein.csv", header = T)
Genelist <- as.list(Genelist)

Hep_M <- SCTransform(Hep_M)

Hep_M <- AddModuleScore(Hep_M, features = Genelist, assay = "SCT")
colnames(Hep_M@meta.data)[length(colnames(Hep_M@meta.data))] <- 'Mouse_Liver_Secretion'
VlnPlot(Hep_M,features = "Mouse_Liver_Secretion",group.by = "seurat_clusters",pt.size = 0,sort = T)
VlnPlot(Hep_M,features = c("nCount_RNA","nFeature_RNA","percent.mito"), pt.size = 0, group.by = "seurat_clusters", ncol = 1)
#-------------------------------------------------------------------------------
##Correlation_Mouse
sc <- qread("../Figure 3/Zonation.qs")
sc$Heps <- as.character(sc$Subtype)
sc$Heps[sc$Heps!="IS"&sc$Heps!="MitoHigh"] <- "Hep"
#-------------------------------------------------------------------------------
Reference_Mouse <- AverageExpression(Hep_M, group.by = "seurat_clusters", features = Hep_M@assays[["SCT"]]@var.features, assays = "SCT", layer = 'scale.data')
Reference_Mouse <- as.data.frame(Reference_Mouse)

Query_Mouse <- AverageExpression(sc, group.by = "Heps", features = sc@assays[["SCT"]]@var.features, assays = "SCT", layer = 'scale.data')
Query_Mouse <- as.data.frame(Query_Mouse)

overlap <- intersect(rownames(Reference_Mouse),rownames(Query_Mouse))
Reference_Mouse <- Reference_Mouse[overlap,]
Query_Mouse <- Query_Mouse[overlap,]

table(rownames(Reference_Mouse)==rownames(Query_Mouse))

Correlation <- cbind(Reference_Mouse,Query_Mouse)
#-------------------------------------------------------------------------------
comp <- cor(Correlation,method = "pearson")

comp <- comp[1:9,10:12]

rownames(comp) <- c("C0","C1","C2","C3","C4","C5","C6","C7","C8")
colnames(comp) <- c("Hep","IS","MitoHigh")

pheatmap::pheatmap(comp,cluster_rows = F, cluster_cols = F, angle_col = 0)

write.csv(comp,"Fig.5E_Correlation.csv")
#-------------------------------------------------------------------------------
########
#Fig.5F
########
#clean environmental variables in case of unnecessary troubles
rm(list =ls()) 
gc()

sp <- readRDS("Fig.5F_DY1_D0_stereo-seq.rds")
DimPlot(sp, reduction = "spatial", group.by = "annotation")
#-------------------------------------------------------------------------------
########
#Fig.5G
########

sp <- CreateSeuratObject(sp@assays$Spatial@counts,meta.data = sp@meta.data)
sp <- SCTransform(sp)

sp <- qread("Fig.5F_DY1_D0_SCTransform.qs")

Markers <- FindAllMarkers(sp, only.pos = T, logfc.threshold = 0.25)
Markers <- subset(Markers, p_val_adj<0.05)

write.csv(Markers,"Fig.5F_Markers.csv")
#-------------------------------------------------------------------------------
#ModuleScore
Genelist <- read.csv("../Figure 4/Fig.4F_Plasma protein.csv", header = T)
Genelist <- as.data.frame(Genelist[,1], ncol = 1)
Genelist <- as.list(Genelist)

Idents(sp) <- "annotation"

sp <- AddModuleScore(sp, features = Genelist, assay = "SCT")
colnames(sp@meta.data)[length(colnames(sp@meta.data))] <- 'Secretion_Score'

VlnPlot(sp,features = c("Secretion_Score"), pt.size = 0, 
        group.by = "annotation", sort = F)+NoLegend()

write.csv(sp@meta.data,"Fig.5G_MetaData.csv")
#-------------------------------------------------------------------------------
kruskal.test(sp$Secretion_Score ~ sp$annotation)
DescTools::DunnTest(sp$Secretion_Score ~ sp$annotation)
#-------------------------------------------------------------------------------
aggregate(sp$Secretion_Score, by = list(type = sp$annotation), FUN = mean)
aggregate(sp$Secretion_Score, by = list(type = sp$annotation), FUN = sd)
#-------------------------------------------------------------------------------
library(ggpubr)

ggviolin(sp@meta.data, x="annotation", y="Secretion_Score", width = 0.6,
         color = "black",
         fill="annotation",
         palette = "npg",
         add = 'mean_sd',
         xlab = F,
         bxp.errorbar=T,
         bxp.errorbar.width=0.5,
         size=1,
         outlier.shape=NA,
         legend = "right",
         title = "ST_Liver_Secretion")
#-------------------------------------------------------------------------------
#Another choice for data visulization
x <- as.data.frame(aggregate(sp$Secretion_Score, by = list(type = sp$annotation), FUN = mean))
x <- x[order(x$x,decreasing = F),]

x$Zonation[1:3] <- "PC zone"
x$Zonation[4:6] <- "MZ zone"
x$Zonation[7:9] <- "PP zone"
x$Zonation <- factor(x$Zonation,levels = c("PP zone","MZ zone","PC zone"))

x$color[x$Zonation=="PC zone"] <- "red"
x$color[x$Zonation=="MZ zone"] <- "blue"
x$color[x$Zonation=="PP zone"] <- "darkgreen"

dotchart(x$x,
         labels = x$type,
         cex = .7,
         groups = x$Zonation,
         gcolor = "black",
         color = x$color,
         pch = 19,
         main = "Average secretion score of each zonation layer",
         xlab = "")
#-------------------------------------------------------------------------------
########
#Fig.5H
########
Reference_Mouse <- AverageExpression(sp, group.by = "annotation", features = sp@assays[["SCT"]]@var.features, assays = "SCT", layer = 'scale.data')
Reference_Mouse <- as.data.frame(Reference_Mouse)

sc <- qread("../Figure 3/Zonation.qs")
Query_Mouse <- AverageExpression(sc, group.by = "Subtype", features = sc@assays[["SCT"]]@var.features, assays = "SCT", layer = 'scale.data')
Query_Mouse <- as.data.frame(Query_Mouse)

overlap <- intersect(rownames(Reference_Mouse),rownames(Query_Mouse))
Reference_Mouse <- Reference_Mouse[overlap,]
Query_Mouse <- Query_Mouse[overlap,]

table(rownames(Reference_Mouse)==rownames(Query_Mouse))

Correlation <- cbind(Reference_Mouse,Query_Mouse)

comp <- cor(Correlation,method = "pearson")

comp <- comp[1:9,10:18]

rownames(comp) <- c("L1","L2","L3","L4","L5","L6","L7","L8","L9")
colnames(comp) <- c("IS","MitoHigh","Hep_Pcdh11x+","Hep_Egr1+","Hep_Car3+","Hep_Cidec+","Hep_PP","Hep_MZ",
                    "Hep_PC")

pheatmap::pheatmap(t(comp),cluster_rows = T, cluster_cols = F, angle_col = 0)

write.csv(comp,"Fig.5H_Correlation.csv")
#-------------------------------------------------------------------------------