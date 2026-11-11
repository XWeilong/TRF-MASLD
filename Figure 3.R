#clean environmental variables in case of unnecessary troubles
rm(list =ls()) 

#set working path
setwd("C:/Users/xiaoweilong/Desktop/Now!/Code/Figure 3/")
getwd()

#load needed R package
if(!require(Seurat))install.packages("Seurat")
if(!require(qs))install.packages("qs")
if(!require(Augur))devtools::install_github("neurorestore/Augur")
if(!require(scRNAtoolVis))devtools::install_github('junjunlab/scRNAtoolVis')
#-------------------------------------------------------------------------------
########
#Fig.3A
########
Integration <- qread("../Figure 2/Integration.qs") #dataset for all cells
#-------------------------------------------------------------------------------
#TRF
Integration <- subset(Integration, Group!="ND")
#-------------------------------------------------------------------------------
#HFD
Integration <- subset(Integration, Group!="TRF")
#-------------------------------------------------------------------------------
DefaultAssay(Integration) <- "RNA"
#Integration <- split(Integration,f = Integration$orig.ident)
Integration <- NormalizeData(Integration)
Integration <- ScaleData(Integration)
#-------------------------------------------------------------------------------
DimPlot(Integration, label = T, group.by = "CellType")+NoLegend()

Integration[["RNA"]] <- JoinLayers(Integration[["RNA"]])

table(Integration$Group)

augur = calculate_auc(Integration, cell_type_col = "CellType", label_col = "Group", show_progress = T)
#-------------------------------------------------------------------------------
AUC_HFD <- as.data.frame(augur[["AUC"]])
write.csv(AUC_HFD,"Fig.3A_Augur_HFD.csv")

AUC_TRF <- as.data.frame(augur[["AUC"]])
write.csv(AUC_TRF,"Fig.3A_Augur_TRF.csv")
#-------------------------------------------------------------------------------
data <- augur$AUC
data$cell_type <- factor(data$cell_type,levels = rev(data$cell_type))
data <- data[order(data$auc),]

names(data)

data$Category <- paste0(round(data$auc,2)," – ",data$cell_type)

data$Category

ggplot(data = data, mapping = aes(x = reorder(cell_type,auc), y = auc, fill = cell_type)) +
  
  geom_bar(stat = "identity", position = "dodge") +
  
  ylim(c(0, 1.3))+
  
  xlab("") + ylab("") +
  
  coord_polar(theta = "y")+
  
  geom_text(data = data, hjust = 0, size = 3,
            
            aes(x = reorder(Category,auc, decreasing =T), y = 0, label = Category)) +
  
  theme_minimal() +
  
  theme(legend.position = "none",
        
        panel.grid.major = element_blank(),
        
        panel.grid.minor = element_blank(),
        
        axis.line = element_blank(),
        
        axis.text.y = element_blank(),
        
        axis.text.x = element_blank(),
        
        axis.ticks = element_blank())
#-------------------------------------------------------------------------------
########
#Fig.3B
########
#Illustrations created with BioRender.com.
#-------------------------------------------------------------------------------
########
#Fig.3C
########
Integration <- qread("Zonation.qs")
#-------------------------------------------------------------------------------
DefaultAssay(Integration) <- "RNA"
Integration@assays$SCT <- NULL
Integration <- JoinLayers(Integration,assay = "RNA")
Integration[["RNA"]] <- split(Integration[["RNA"]],f=Integration$orig.ident)

Integration <- SCTransform(Integration)
Integration <- RunPCA(Integration)

Integration <- IntegrateLayers(object = Integration, method = HarmonyIntegration,
                               orig.reduction = "pca", new.reduction = 'harmony',
                               assay = "SCT", verbose = TRUE)

Integration <- RunUMAP(Integration, reduction = "harmony", dims = 1:30)
#-------------------------------------------------------------------------------
Idents(Integration) <- "Subtype"
Integration$CellType <- factor(Integration$Subtype, levels = c("IS","Hep_Cidec+","Hep_Car3+","Hep_Pcdh11x+","Hep_Egr1+","MitoHigh","Hep_PP","Hep_MZ","Hep_PC"))
DimPlot(Integration, group.by = "Subtype", label = T, label.size = 3, repel = T)

#Subtype markers
Integration <- PrepSCTFindMarkers(Integration)
Markers <- FindAllMarkers(Integration, logfc.threshold = 0.5, return.thresh = 0.05, 
                         test.use = "wilcox", min.pct = 0.5, only.pos = T)

Markers <- subset(Markers, Markers$p_val_adj < 0.05)
write.csv(Markers,"Fig.3C_Markers.csv")
#--------------------------------------------------------------------
########
#Fig.3D
########
Idents(Integration) <- "Subtype"

Integration$Subtype <- factor(Integration$Subtype, levels = c("IS","MitoHigh",
        "Hep_Pcdh11x+","Hep_Egr1+","Hep_Car3+","Hep_Cidec+","Hep_PP","Hep_MZ","Hep_PC"))

cellRatioPlot(object = Integration,
              sample.name = "Group",
              celltype.name = "Subtype",
              col.width = 0.6, flow.alpha = 0.15,
              flow.curve = 0.5)

CellRatio <- as.data.frame(prop.table(table(Integration$Group,Integration$Subtype),margin = 1))
write.csv(CellRatio,"Fig.3D_CellRatio.csv")
#--------------------------------------------------------------------
########
#Fig.3E
########
Genelist <- readxl::read_excel("Fig.3E_Landmarkers.xlsx")

#RNA-counts based
DefaultAssay(Integration) <- "RNA"

CV <- as.data.frame(t(Integration@assays$RNA$counts[Genelist$Central[1:20],]))
PV <- as.data.frame(t(Integration@assays$RNA$counts[Genelist$Portal[1:15],]))

table(rownames(CV)==rownames(PV))
Correlation <- cbind(CV,PV)

comp <- cor(Correlation,method = "spearman") #pearson or spearman

pheatmap::pheatmap(comp,cluster_rows = F, cluster_cols = F, show_rownames = T, 
                   show_colnames = T, angle_col = 90, fontsize_col = 8,fontsize_row = 8 )

write.csv(comp, "Fig.3E_Correlation.csv")
#--------------------------------------------------------------------
########
#Fig.3F
########
seu <- SplitObject(Integration, split.by = "orig.ident")
df <- list()

for (i in 1:9){
  seu[[i]]@meta.data <- seu[[i]]@meta.data[,c(1:9)]
  DefaultAssay(seu[[i]]) <- "RNA"
  seu[[i]] <- SCTransform(seu[[i]])
}

Genelist <- readxl::read_excel("Fig.3E_Landmarkers.xlsx")
Genelist <- as.data.frame(Genelist[,c(1,2)], ncol = 2)
Genelist <- as.list(Genelist)

for (s in names(seu)) {
  for (i in names(Genelist)){
    seu[[s]] <- AddModuleScore(seu[[s]], assay = "SCT", features = Genelist[i], 
                               nbin = 12, name = paste0(i,"_Score"))
    colnames(seu[[s]]@meta.data)[ncol(seu[[s]]@meta.data)] <- paste0(i,"_Score")
  }
#-------------------------------------------------------------------------------

  # Perform scaling to (0, 1)
  seu[[s]]$Central_Score <- scales::rescale(seu[[s]]$Central_Score, to = c(0,1))
  seu[[s]]$Portal_Score <- scales::rescale(seu[[s]]$Portal_Score, to = c(0,1))
  seu[[s]]$Zonation_Score <- ((seu[[s]]$Portal_Score) / (seu[[s]]$Central_Score + seu[[s]]$Portal_Score + 0.01))
  
  df[[s]] <- data.table::data.table(sample = s, value = seu[[s]]$Zonation_Score, 
                              CellType=seu[[s]]$Subtype, Group = seu[[s]]$Group)
  x <- df[[s]]$value
  
  idx <- x < 1e-2
  if (sum(idx > 0)) {
    x[idx] <- 1e-2 + stats::runif(n = sum(idx), min = 0, max = 1e-4)
  }
  idx <- x > 1 - 1e-2
  if (sum(idx > 0)) {
    x[idx] <- 1 - 1e-2 - stats::runif(n = sum(idx), min = 0, max = 1e-4)
  }
  df[[s]]$value <- x
}

# Join all samples
df_plot <- data.table::rbindlist(df) |> as.data.frame()
#-------------------------------------------------------------------------------

condition_col_pal <- c("#38639E","#BEA881","#E7E2DA")

print(ggplot(df_plot, aes(x = value, fill = Group, colour = Group)) +
        geom_histogram(aes(y=..density..), bins = 50, alpha = 0.4, position = "identity") + #..count../sum(..count..))
        scale_color_manual(values = condition_col_pal) +
        scale_fill_manual(values = condition_col_pal) +
        xlab("Zonation specificity score") +
        ylab(NULL) + theme_classic() + theme(legend.position = c(0.9, 0.85)))
#-------------------------------------------------------------------------------
#Kernel Density Plot
library(sm)
par(lwd=2)

zonation.f <- factor(df_plot$Group,levels = c("ND","HFD","TRF"),
                     labels = c("ND","HFD","TRF"))

sm.density.compare(df_plot$value,df_plot$Group,xlab="Zonation score",
                   col=c("#38639E","#BEA881","#E7E2DA"))

title(main = "Hepatocyte Distribution")
legend(locator(1),levels(zonation.f),fill = c("#38639E","#BEA881","#E7E2DA"))
#-------------------------------------------------------------------------------
qsave(seu,"Fig.3F_Zonation.qs")
#-------------------------------------------------------------------------------
library(Seurat)
library(ggplot2)

PCAPlot(Integration)

# 提取 PC2
pca_emb <- Embeddings(Integration, reduction = "pca")

# 构造 data.frame（关键：所有变量都放在一个数据框里）
df <- data.frame(
  Zonation_Score = Integration$Zonation_Score,
  PC2            = pca_emb[, 2]
)

# 绘图
ggplot(df, aes(x = Zonation_Score, y = PC2)) +
  geom_point(size = 0.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
  theme_classic() +
  labs(#title = "Zonation Score vs PC2",
       x = "Zonation Score",
       y = "PC2")


# 1. 提取 PCA 载荷
pca_loadings <- Loadings(Integration, reduction = "pca")

# 2. 确认 PC_2 列名（根据你之前的情况，应该是 "PC_2"）
pc2_col <- grep("PC.?2$", colnames(pca_loadings), value = TRUE)
cat("使用的列:", pc2_col, "\n")

# 3. 构造 data.frame，按 PC_2 loading 降序排列
pc2_df <- data.frame(
  gene    = rownames(pca_loadings),
  loading = as.numeric(pca_loadings[, pc2_col])
)
pc2_df <- pc2_df[order(pc2_df$loading, decreasing = TRUE), ]
pc2_df$rank <- seq_len(nrow(pc2_df))

# 4. 标记目标基因
target_gene <- c("Sds","Hal","Cyp2f2","Gls2","Hsd17b13","Sfxn1","Aspg","Uroc1",
                 "Esr1","Amdhd1","Gldc","Noct","Asl","Ftcd","Slc25a51","Slc1a2",
                 "Glul","Cyp2c29","Tbx3","Oat","Gulo","Slc16a10","Cyp2c67","Akr1c6",
                 "Rnf43","Cyp7a1","Cyp2c54","Aldh3a2","Cyp7b1","Slco1b2","Slc25a21",
                 "Cyp27a1","Cyp2e1","Gpcpd1","Sult2a8")

pc2_df$is_target <- pc2_df$gene %in% target_gene

# 5. 查看目标基因的排名
target_ranks <- pc2_df[pc2_df$is_target, c("gene", "loading", "rank")]
print(target_ranks)


library(ggplot2)
library(ggrepel)

# 从目标基因中挑出 loading 最大和最小的各 3 个
target_df <- pc2_df[pc2_df$is_target, ]
target_top3 <- target_df[order(target_df$loading, decreasing = TRUE), ][1:3, ]
target_bot3 <- target_df[order(target_df$loading, decreasing = FALSE), ][1:3, ]
label_df <- rbind(target_top3, target_bot3)

ggplot(pc2_df, aes(x = rank, y = loading)) +
  geom_line(color = "grey50") +
  # 全部目标基因点（可选，保留红点）
  geom_point(data = target_df,
             aes(x = rank, y = loading),
             color = "red", size = 1.5) +
  # 只给前3和后3加标签
  geom_text_repel(
    data = label_df,
    aes(label = gene),
    color = "red", size = 3, max.overlaps = 20
  ) +
  theme_classic() +
  labs(x = "Gene rank",
       y = "PC_2 Loadings")
#-------------------------------------------------------------------------------
########
#Fig.3G
########
#Zonation score
Genelist <- readxl::read_excel("Fig.3E_Landmarkers.xlsx")

CV <- as.data.frame(Genelist[,1], ncol = 1)
CV <- as.list(CV)

PV <- as.data.frame(Genelist[,2], ncol = 1)
PV <- as.list(PV)

DefaultAssay(Integration) <- "SCT"

#Central score
Integration <- AddModuleScore(Integration, assay = "SCT", features = CV, name = "Central Score")
colnames(Integration@meta.data)[ncol(Integration@meta.data)] <- "Central_Score"

#Portal score
Integration <- AddModuleScore(Integration, assay = "SCT", features = PV, name = "Portal Score")
colnames(Integration@meta.data)[ncol(Integration@meta.data)] <- "Portal_Score"

VlnPlot(Integration,features = c("Central_Score","Portal_Score"), sort = T, group.by =  "Subtype", pt.size = 0, ncol=1) +NoLegend()

#DataScale
Integration$CS_Scale <- scales::rescale(Integration$Central_Score, to = c(0,1))
Integration$PS_Scale <- scales::rescale(Integration$Portal_Score, to = c(0,1))

#Zonation Score
Integration$Zonation_Score <- ((Integration$PS_Scale) / (Integration$CS_Scale + Integration$PS_Scale + 0.01))

VlnPlot(Integration,features = "Zonation_Score", group.by =  "Subtype", sort = T, pt.size = 0) +NoLegend()

write.csv(Integration@meta.data,"Fig.3G_Zonation.csv")
#-------------------------------------------------------------------------------
qsave(Integration,"Zonation.qs")
#-------------------------------------------------------------------------------
#########
#Fig.3H&I
#########
Genelist <- read.csv("Fig.3H-I_Geneset.csv", header = T)
Genelist <- as.data.frame(Genelist[,1], ncol = 1)
Genelist <- as.list(Genelist)
DefaultAssay(Integration) <- "SCT"
Integration <- AddModuleScore(Integration, assay = "SCT", features = Genelist, name = "Xenobiotic metabolic process")
colnames(Integration@meta.data)[ncol(Integration@meta.data)] <- "Xenobiotic metabolic process"
#-------------------------------------------------------------------------------
Genelist <- read.csv("Fig.3H-I_Geneset.csv", header = T)
Genelist <- as.data.frame(Genelist[,2], ncol = 1)
Genelist <- as.list(Genelist)
DefaultAssay(Integration) <- "SCT"
Integration <- AddModuleScore(Integration, assay = "SCT", features = Genelist, name = "Amino acid catabolic process")
colnames(Integration@meta.data)[ncol(Integration@meta.data)] <- "Amino acid catabolic process"
#-------------------------------------------------------------------------------
Hep <- Integration[,Integration$Subtype=="Hep_PC"|Integration$Subtype=="Hep_MZ"|Integration$Subtype=="Hep_PP"]

Hep$Subtype <- factor(Hep$Subtype, levels = c("Hep_PC","Hep_MZ","Hep_PP"))

p1 <- VlnPlot(Hep,features = "Xenobiotic metabolic process", group.by = "Subtype", sort = F, pt.size = 0) +NoLegend()
p2 <- VlnPlot(Hep,features = "Amino acid catabolic process", group.by = "Subtype", sort = F, pt.size = 0) +NoLegend()

p1|p2
#-------------------------------------------------------------------------------
#P value adjustment method: holm
kruskal.test(Hep$`Xenobiotic metabolic process` ~ Hep$Group)
DescTools::DunnTest(Hep$`Xenobiotic metabolic process` ~ Hep$Group)
#-------------------------------------------------------------------------------
kruskal.test(Hep$`Amino acid catabolic process` ~ Hep$Group)
DescTools::DunnTest(Hep$`Amino acid catabolic process` ~ Hep$Group)
#-------------------------------------------------------------------------------
#Statistic
aggregate(Hep$`Xenobiotic metabolic process`, by = list(type = Hep$Subtype), FUN = mean)
aggregate(Hep$`Xenobiotic metabolic process`, by = list(type = Hep$Subtype), FUN = sd)

aggregate(Hep$`Amino acid catabolic process`, by = list(type = Hep$Subtype), FUN = mean)
aggregate(Hep$`Amino acid catabolic process`, by = list(type = Hep$Subtype), FUN = sd)
#-------------------------------------------------------------------------------