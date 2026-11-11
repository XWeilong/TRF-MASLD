#clean environmental variables in case of unnecessary troubles
rm(list =ls()) 
gc()

#set working path
setwd("C:/Users/xiaoweilong/Desktop/Now!/Code/Figure 6/")
getwd()

library("Seurat")
library("qs")
library("dior")
library("dplyr")
library("scRNAtoolVis")


Integration <- qread("../Figure 3/Zonation.qs")

#------------------------------------FindMarkers--------------------------------
#Subtype-specific DEGs
Integration$Subtype_Group <- paste(Integration$Subtype, Integration$Group, sep = "_")
Idents(Integration) <- "Subtype_Group"
table(Integration$Subtype_Group)
#-------------------------------------------------------------------------------
Hep_PC <- FindMarkers(Integration, ident.1 = "Hep_PC_TRF", ident.2 = "Hep_PC_HFD", min.pct = 0.25, logfc.threshold = 0.5)
Hep_PC <- subset(Hep_PC, p_val_adj < 0.05)
Hep_PC$Direction[Hep_PC$avg_log2FC>0] <- "up"
Hep_PC$Direction[Hep_PC$avg_log2FC<0] <- "down"
#-------------------------------------------------------------------------------
Hep_MZ <- FindMarkers(Integration, ident.1 = "Hep_MZ_TRF", ident.2 = "Hep_MZ_HFD", min.pct = 0.25, logfc.threshold = 0.5)
Hep_MZ <- subset(Hep_MZ, p_val_adj < 0.05)
Hep_MZ$Direction[Hep_MZ$avg_log2FC>0] <- "up"
Hep_MZ$Direction[Hep_MZ$avg_log2FC<0] <- "down"
#-------------------------------------------------------------------------------
Hep_PP <- FindMarkers(Integration, ident.1 = "Hep_PP_TRF", ident.2 = "Hep_PP_HFD", min.pct = 0.25, logfc.threshold = 0.5)
Hep_PP <- subset(Hep_PP, p_val_adj < 0.05)
Hep_PP$Direction[Hep_PP$avg_log2FC>0] <- "up"
Hep_PP$Direction[Hep_PP$avg_log2FC<0] <- "down"
#-------------------------------------------------------------------------------
MitoHigh <- FindMarkers(Integration, ident.1 = "MitoHigh_TRF", ident.2 = "MitoHigh_HFD", min.pct = 0.25, logfc.threshold = 0.5)
MitoHigh <- subset(MitoHigh, p_val_adj < 0.05)
MitoHigh$Direction[MitoHigh$avg_log2FC>0] <- "up"
MitoHigh$Direction[MitoHigh$avg_log2FC<0] <- "down"
#-------------------------------------------------------------------------------
Hep_Cidec <- FindMarkers(Integration, ident.1 = "Hep_Cidec+_TRF", ident.2 = "Hep_Cidec+_HFD", min.pct = 0.25, logfc.threshold = 0.5)
Hep_Cidec <- subset(Hep_Cidec, p_val_adj < 0.05)
Hep_Cidec$Direction[Hep_Cidec$avg_log2FC>0] <- "up"
Hep_Cidec$Direction[Hep_Cidec$avg_log2FC<0] <- "down"
#-------------------------------------------------------------------------------
Hep_Car3 <- FindMarkers(Integration, ident.1 = "Hep_Car3+_TRF", ident.2 = "Hep_Car3+_HFD", min.pct = 0.25, logfc.threshold = 0.5)
Hep_Car3 <- subset(Hep_Car3, p_val_adj < 0.05)
Hep_Car3$Direction[Hep_Car3$avg_log2FC>0] <- "up"
Hep_Car3$Direction[Hep_Car3$avg_log2FC<0] <- "down"
#-------------------------------------------------------------------------------
Hep_Pcdh11x <- FindMarkers(Integration, ident.1 = "Hep_Pcdh11x+_TRF", ident.2 = "Hep_Pcdh11x+_HFD", min.pct = 0.25, logfc.threshold = 0.5)
Hep_Pcdh11x <- subset(Hep_Pcdh11x, p_val_adj < 0.05)
Hep_Pcdh11x$Direction[Hep_Pcdh11x$avg_log2FC>0] <- "up"
Hep_Pcdh11x$Direction[Hep_Pcdh11x$avg_log2FC<0] <- "down"
#-------------------------------------------------------------------------------
Hep_Egr1 <- FindMarkers(Integration, ident.1 = "Hep_Egr1+_TRF", ident.2 = "Hep_Egr1+_HFD", min.pct = 0.25, logfc.threshold = 0.5)
Hep_Egr1 <- subset(Hep_Egr1, p_val_adj < 0.05)
Hep_Egr1$Direction[Hep_Egr1$avg_log2FC>0] <- "up"
Hep_Egr1$Direction[Hep_Egr1$avg_log2FC<0] <- "down"
#-------------------------------------------------------------------------------
IS <- FindMarkers(Integration, ident.1 = "IS_TRF", ident.2 = "IS_HFD", min.pct = 0.25, logfc.threshold = 0.5)
IS <- subset(IS, p_val_adj < 0.05)
IS$Direction[IS$avg_log2FC>0] <- "up"
IS$Direction[IS$avg_log2FC<0] <- "down"
#-------------------------------------------------------------------------------
#unique DEGs
a <- union(rownames(Hep_PC),rownames(Hep_MZ))
a <- union(a,rownames(Hep_PP))
a <- union(a,rownames(MitoHigh))
a <- union(a,rownames(Hep_Cidec))
a <- union(a,rownames(Hep_Car3))
a <- union(a,rownames(Hep_Egr1))
a <- union(a,rownames(Hep_Pcdh11x))
a <- union(a,rownames(IS))
length(a)
#-------------------------------pySCENIC workflow-------------------------------
#RNA-counts
data <- Integration
DefaultAssay(data) <- "RNA"
data@assays$SCT <- NULL

data <- data[a,] #Subtype-specific DEGs

data[["RNA"]] <- as(object = data[["RNA"]], Class = "Assay")

# Convert SeuratObj to h5 by dior
write_h5(data,
         file = "TRF_Liver_pySCENIC_250825.h5", 
         object.type = 'seurat', 
         assay.name = 'RNA',
         save.graphs = F,
         save.scale = F)
#-------------------------------------------------------------------------------
##pySCENIC in python!
#-------------------------------------------------------------------------------
########
#Fig.6A
########

df <- data.table::fread("pySCENIC/auc_mtx.csv", stringsAsFactors = TRUE) %>% as.data.frame()
rownames(df) <- df[,1]
df <- df[,-1]
df[1:3, 1:3]

TFs <- gsub("[(+)]*", "", colnames(df))

meta_scenic <- Integration@meta.data

mtx <- matrix(0, nrow = 9 * 3, ncol = ncol(df))
rownames(mtx) <- paste0(c("Hep_PC", "Hep_MZ", "Hep_PP","MitoHigh",
                          "Hep_Cidec+","Hep_Car3+","Hep_Pcdh11x+","Hep_Egr1+","IS"), "_",rep(c("ND", "HFD", "TRF"), each = 9))
colnames(mtx) <- colnames(df)
mtx <- mtx %>% t()
#-------------------------------------------------------------------------------
# Mean value of each CellType per Group
for (i in c("Hep_PC", "Hep_MZ", "Hep_PP","MitoHigh",
            "Hep_Cidec+","Hep_Car3+","Hep_Pcdh11x+","Hep_Egr1+","IS")) {
  CellType.df <- filter(meta_scenic, Subtype == i)
  
  for (g in rownames(mtx)) {
    CellType.df$value <- df[rownames(CellType.df), g]
    
    r <- CellType.df %>%
      group_by(Group) %>%
      summarise(across(value, mean, na.rm = TRUE)) %>%
      ungroup()
    
    mtx[g, paste0(i, "_", r$Group)] <- r$value
  }
}

#options(repr.plot.width =10, repr.plot.height = 20)
mtx <- mtx[,c(1,10,19,2,11,20,3,12,21,4,13,22,5,14,23,6,15,24,7,16,25,8,17,26,9,18,27)]
#-------------------------------------------------------------------------------
annotation_col <- data.frame(row.names = colnames(mtx), CellType = rep(c("Hep_PC", "Hep_MZ", "Hep_PP","MitoHigh",
                                                                         "Hep_Cidec+","Hep_Car3+","Hep_Pcdh11x+","Hep_Egr1+","IS"), each = 3))

annotation_col$CellType <- factor(annotation_col$CellType, levels = c("Hep_PC", "Hep_MZ", "Hep_PP","MitoHigh",
                                                                      "Hep_Cidec+","Hep_Car3+","Hep_Pcdh11x+","Hep_Egr1+","IS")) 

ann_colors <- list(CellType = c("#E3AD68","#edd05e","#F6F5B4","#ACD45E","#2F528F","#E54924","#9ECABE","#1a918f", "#e07233"))
names(ann_colors$CellType) <- c("Hep_PC", "Hep_MZ", "Hep_PP","MitoHigh",
                                "Hep_Cidec+","Hep_Car3+","Hep_Pcdh11x+","Hep_Egr1+","IS")
#-------------------------------------------------------------------------------
pheatmap::pheatmap(mtx, show_rownames = T, angle_col = 45,
                   cluster_cols = F, scale = "row", 
                   gaps_col = c(1:9) * 3,
                   labels_col = rep(c("ND", "HFD","TRF"), 9),
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(200), border_color = NA,
                   annotation_col = annotation_col,
                   annotation_colors = ann_colors,
                   cellheight = 12, cellwidth = 14, fontsize = 6, main = "SCENIC Regulons")

write.csv(mtx,"Fig.6A_Regulon.csv")
#-------------------------------------------------------------------------------
########
#Fig.6B
########
library(Mfuzz)

expression_matrix <- t(df)

group_info <- Integration@meta.data$Group

grouped_matrix <- sapply(unique(group_info), function(Group) {
  rowMeans(expression_matrix[, group_info == Group, drop = FALSE])
})

colnames(grouped_matrix) <- c("ND","HFD","TRF")

eset <- ExpressionSet(assayData = as.matrix(grouped_matrix))
eset <- standardise(eset)
#-------------------------------------------------------------------------------
errors <- sapply(2:10, function(c) {
  cl <- mfuzz(eset, c = c, m = 1.25)
  sum(cl$withinerror)
})

plot(2:10, errors, type = "b", xlab = "Number of clusters (c)", ylab = "Clustering error")

cl <- mfuzz(eset, c = 6, m = 1.25)

mfuzz.plot(eset, cl, mfrow = c(3, 2), new.window = FALSE)
#-------------------------------------------------------------------------------
cluster_results <- cl$cluster

membership <- cl$membership
gene_clusters <- apply(membership, 1, which.max)
cluster_genes <- split(names(gene_clusters), gene_clusters)
print(cluster_genes)
#-------------------------------------------------------------------------------
write.csv(as.data.frame(eset@assayData[["exprs"]]),"Fig.6B_Exprs.csv")
write.csv(as.data.frame(cl[["centers"]]),"Fig.6B_Centers.csv")
write.csv(as.data.frame(cl[["membership"]]),"Fig.6B_Membership.csv")
#-------------------------------------------------------------------------------
########
#Fig.6C
########
gene=c("Clock","Mlxipl","Onecut1","Bach2","Srebf2",
       "Srebf1","Bhlhe40","Foxa3","Tef","Esrra","Dbp")

options(repr.plot.width=18, repr.plot.height=10)
p=DotPlot(Integration,features=gene, dot.scale=11, scale.min = 0, scale.max = 100, group.by = "Group")+scale_y_discrete("")+scale_x_discrete("")#+coord_flip()
p1=p+theme(panel.grid=element_blank())+theme_bw(base_size = 13)+#RotatedAxis()+
  scale_color_viridis_c()+labs(x=NULL,y=NULL)+theme(aspect.ratio=0.3)+ggtitle("")
plot(p1)
#-------------------------------------------------------------------------------
########
#Fig.6D
########
#Metascape.
#-------------------------------------------------------------------------------
########
#Fig.6E
########

#GO biological process enrichment
Targets <- read.csv("Fig.6E_Bhlhe40_targets.csv")

library(clusterProfiler)
library(org.Mm.eg.db)

GO <- enrichGO(gene=Targets$gene, keyType = "SYMBOL", ont = "BP", 'org.Mm.eg.db')
#GO <- simplify(GO)

dotplot(GO, showCategory = 10, font.size = 10)

cnetplot(GO, circular = TRUE, colorEdge = TRUE, cex_label_category = 0.6,cex_label_gene = 0.5,
         showCategory = c("circadian rhythm", "response to nutrient levels", "lipid homeostasis","carbohydrate catabolic process",
                                                                 "fatty acid metabolic process","triglyceride metabolic process","cholesterol metabolic process"))

write.csv(GO@result,"Fig.6E_GO terms.csv")
#-------------------------------------------------------------------------------
########
#Fig.6F
########
#Prism Graphpad.
#-------------------------------------------------------------------------------