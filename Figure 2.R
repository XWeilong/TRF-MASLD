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
Integration <- qread("Integration.qs") #dataset for all cells
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

Integration <- RunUMAP(Integration, reduction = "harmony", dims = 1:35)
#-------------------------------------------------------------------------------
########
#Fig.2A
########
Idents(Integration) <- "CellType"
DimPlot(Integration, group.by = "CellType", label = T, label.size = 3, repel = T)

#Number of cells
table(Integration$orig.ident)

#Quality control
qc_summary <- Integration@meta.data %>%
  group_by(orig.ident) %>%
  summarise(
    Avg_UMIs = mean(nCount_RNA),
    Avg_Genes = mean(nFeature_RNA),
    Avg_MT = mean(percent.mt)
  )
#-------------------------------------------------------------------------------
########
#Fig.2B
########
table(Integration$orig.ident)
Idents(Integration) <- "CellType"
Integration <- PrepSCTFindMarkers(Integration)
Markers <- FindAllMarkers(Integration, logfc.threshold = 1, return.thresh = 0.05, 
                         test.use = "wilcox", min.pct = 0.5, only.pos = T)

Markers <- subset(Markers, Markers$p_val_adj < 0.05)
write.csv(Markers,"Fig.2B_CellType_Markers.csv")
#-------------------------------------------------------------------------------
#Celltype-specific markers
Hep=c("Hnf4a","Cyp2e1","Hamp","Slc7a2","Apob")
MitoHigh=c("Mup20","Scd1","Apof","Igfbp2")
IS=c("Rsad2","Ifit1","Ifit3","Cmpk2")
HSC=c("Reln","Hgf")
Myof=c("Dcn","Thsd4")
EC=c("Ptprb","Stab2")
Cho=c("Pkhd1","Spp1")
Kup=c("Cd5l","Clec4f")  
Lym=c("Bcl2","Skap1")
Mono=c("Ctss","Ccr2")
DC=c("Irf8","Naaa")

genes_to_check =list(
  Hep=Hep,
  MitoHigh=MitoHigh,
  IS=IS,
  HSC=HSC,
  Myof=Myof,
  EC=EC,
  Cho=Cho,
  Kup=Kup,
  Lym=Lym,
  Mono=Mono,
  DC=DC)

#order
Idents(Integration) <- "CellType"
My_levels<-c("DC","Mono","Lym","Kup","Cho","EC","Myof","HSC","IS","MitoHigh","Hep")
Integration@active.ident<-factor(Integration@active.ident,levels= rev(My_levels))

DefaultAssay(Integration) <- "SCT"
p1=DotPlot(Integration,#scale=F,
           features=genes_to_check,cols = c("#B2BEB5","#D72F24"))+
  theme(
    axis.text.x  = element_text(color="black",size=12,angle = 45,vjust = 0.5, hjust=0.5),
    panel.border = element_rect(color="black"),
    panel.spacing = unit(2, "mm"),
    axis.line = element_blank(),
  )+labs(x="", y=""); p1
#-------------------------------------------------------------------------------
########
#Fig.2C
########
cellRatioPlot(object = Integration,
              sample.name = "Group", #orig.ident
              celltype.name = "CellType",
              col.width = 0.6, flow.alpha = 0.15,
              flow.curve = 0.5)

mydata <- as.matrix(table(Integration$CellType,Integration$orig.ident))
write.csv(mydata, "Fig.2C_CellNum.csv")
mydata <- as.matrix(prop.table(table(Integration$CellType, Integration$orig.ident), margin = 2))
write.csv(mydata, "Fig.2C_CellRatio.csv")

table(Integration$orig.ident)
table(Integration$CellType)

#Fig.2C Made by GraphPad Prism 10.
#-------------------------------------------------------------------------------
########
#Fig.2D
########
DefaultAssay(Integration) <- "RNA"
Integration[["RNA"]] <- split(Integration[["RNA"]],f=Integration$orig.ident)

Integration <- NormalizeData(Integration)
Integration <- FindVariableFeatures(Integration)

data <- AggregateExpression(Integration, group.by = "orig.ident", slot = "counts", features = VariableFeatures(Integration), assays = "RNA")

data <- as.data.frame(data)

colnames(data) <- c("HFD1","HFD2","HFD3","ND1","ND2","ND3","TRF1","TRF2","TRF3")

comp <- cor(data, method = "spearman")
pheatmap::pheatmap(comp, cluster_rows = T, cluster_cols = T, angle_col = 0, fontsize_col = 8 )

write.csv(comp, "Fig.2D_Correlation.csv")
#-------------------------------------------------------------------------------
########
#Fig.2E
########
#NAFLD-related gene modules
GeneSet <- read.csv("Fig.2E_GeneSet.csv")
Genelist <- as.data.frame(GeneSet[,])
Genelist <- as.list(Genelist)

ref.1=list(Genelist$Cholesterol.metabolic.process_GO_0008203)
ref.2=list(Genelist$Triglyceride.metabolic.process_GO_0006641)
ref.3=list(Genelist$Fatty.acid.metabolic.process_GO_0006631)
ref.4=list(Genelist$Cellular.response.to.oxidative.stress_GO_0034599)
ref.5=list(Genelist$Cellular.response.to.ER.stress_GO_0034976)
ref.6=list(Genelist$Inflammatory.response_GO_0006954)
ref.7=list(Genelist$Extracellular.matrix.assembly_GO_0085029)
#-------------------------------------------------------------------------------
DefaultAssay(Integration) <- "SCT"

for (i in ls(pattern = "ref.")) {
  Integration <- AddModuleScore(object = Integration, assay = "SCT",
                                features = get(i), name = gsub("ref.", "", i))
}; rm(i)

Idents(Integration) <- "Group"

#diffcellmean module score
diffcellmeanModule=function(pathways,x){
  meta = x@meta.data
  celltypes = unique(meta$CellType)
  res.df  = as.data.frame(matrix(nrow=length(celltypes),ncol=length(pathways)))
  rownames(res.df) <- celltypes
  colnames(res.df)=pathways
  for (type in celltypes) {
    cell.data = subset(x, subset = CellType == type)
    matrix = cell.data@meta.data[,pathways]
    result = apply(matrix, 2, mean)
    res.df[type,pathways] <- result
  }
  return(res.df)
}


ND=subset(Integration,idents='ND')
pathways_ND <- c("11","21","31","41","51","61","71")
info_ND=diffcellmeanModule(pathways_ND,ND)
info2_ND=info_ND[c("Hep","MitoHigh","IS","HSC","Myof","EC","Cho","Kup","Lym","Mono","DC"),]

HFD=subset(Integration,idents='HFD')
pathways_HFD <- c("11","21","31","41","51","61","71")
info_HFD=diffcellmeanModule(pathways_HFD,HFD)
info2_HFD=info_HFD[c("Hep","MitoHigh","IS","HSC","Myof","EC","Cho","Kup","Lym","Mono","DC"),]

TRF=subset(Integration,idents='TRF')
pathways_TRF <- c("11","21","31","41","51","61","71")
info_TRF=diffcellmeanModule(pathways_TRF,TRF)
info2_TRF=info_TRF[c("Hep","MitoHigh","IS","HSC","Myof","EC","Cho","Kup","Lym","Mono","DC"),]

info_total=cbind(t(info2_ND),t(info2_HFD),t(info2_TRF))
info_total

colnames(info_total) <- c("Hep_ND","MitoHigh_ND","IS_ND","HSC_ND","Myof_ND","EC_ND","Cho_ND","Kup_ND","Lym_ND","Mono_ND","DC_ND",
                          "Hep_HFD","MitoHigh_HFD","IS_HFD","HSC_HFD","Myof_HFD","EC_HFD","Cho_HFD","Kup_HFD","Lym_HFD","Mono_HFD","DC_HFD",
                          "Hep_TRF","MitoHigh_TRF","IS_TRF","HSC_TRF","Myof_TRF","EC_TRF","Cho_TRF","Kup_TRF","Lym_TRF","Mono_TRF","DC_TRF")

info_total2 <- info_total[,c(1,12,23,2,13,24,3,14,25,4,15,26,5,16,27,6,17,28,7,18,29,8,19,30,9,20,31,10,21,32,11,22,33)]

row.names(info_total2)=c("Cholesterol metabolic process_GO_0008203",
                         "Triglyceride metabolic process_GO_0006641",
                         "Fatty acid metabolic process_GO_0006631",
                         "Cellular response to oxidative stress_GO_0034599",
                         "Cellular response to ER stress_GO_0034976",
                         "Inflammatory response_GO_0006954",
                         "Extracellular matrix assembly_GO_0085029")

options(repr.plot.width = 10, repr.plot.height = 5)
annotation_col=data.frame(Group = factor(rep(c("ND","HFD","TRF"), 11),
                                         levels=c("ND","HFD","TRF")), row.names = colnames(info_total2))
ann_colors = list(Group = c(ND="#440154FF",HFD= "#414487FF",TRF="#2A788EFF"))  
map=pheatmap::pheatmap(info_total2, show_rownames = T, scale = "row", gaps_col = c(1:11)*3,
                       cluster_rows =F,border_color=NA,annotation_col = annotation_col,annotation_colors=ann_colors,
                       cluster_cols = F,cellheight = 14,cellwidth = 14,fontsize = 7, color = colorRampPalette(c("navy",
                                                                                                                 "white", "firebrick3"))(200),main="NAFLD-related pathways module gene sets")
write.csv(info_total2,"Fig.2E_ModuleScore.csv")
#-------------------------------------------------------------------------------
colnames(Integration@meta.data)[10:16]=c("Cholesterol metabolic process_GO_0008203",
                                        "Triglyceride metabolic process_GO_0006641",
                                        "Fatty acid metabolic process_GO_0006631",
                                        "Cellular response to oxidative stress_GO_0034599",
                                        "Cellular response to ER stress_GO_0034976",
                                        "Inflammatory response_GO_0006954",
                                        "Extracellular matrix assembly_GO_0085029")
#-------------------------------------------------------------------------------
########
#Fig.2F
########
library("CellChat")

Integration <- qread("Integration.qs")
Integration <- PrepSCTFindMarkers(Integration)

# Subset for disease
s_ND  <-  subset(Integration, Group == "ND")
s_HFD <- subset(Integration, Group == "HFD")
s_TRF <- subset(Integration, Group == "TRF")

# Run cellchat
cellchat_ND <- createCellChat(object = s_ND,   assay = "SCT",
                             group.by = "CellType")
cellchat_HFD <- createCellChat(object = s_HFD, assay = "SCT",
                             group.by = "CellType")
cellchat_TRF <- createCellChat(object = s_TRF, assay = "SCT",
                             group.by = "CellType")

# Define ligand-receptor pair database
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)

# Set the used database in the object
cellchat_ND@DB  <- CellChatDB
cellchat_HFD@DB <- CellChatDB
cellchat_TRF@DB <- CellChatDB

# Preprocessing
# subset the expression data of signaling genes for saving computation cost
cellchat_ND  <- subsetData(cellchat_ND)
cellchat_HFD <- subsetData(cellchat_HFD)
cellchat_TRF <- subsetData(cellchat_TRF)

# Do parallel
#future::plan("multisession", workers = 4)

cellchat_ND <- identifyOverExpressedGenes(cellchat_ND)
cellchat_ND <- identifyOverExpressedInteractions(cellchat_ND)
cellchat_HFD <- identifyOverExpressedGenes(cellchat_HFD)
cellchat_HFD <- identifyOverExpressedInteractions(cellchat_HFD)
cellchat_TRF <- identifyOverExpressedGenes(cellchat_TRF)
cellchat_TRF <- identifyOverExpressedInteractions(cellchat_TRF)
#-------------------------------------------------------------------------------
# Compute the communication probability and infer cellular communication network
cellchat_ND <- computeCommunProb(cellchat_ND, type = "triMean")
cellchat_ND <- filterCommunication(cellchat_ND, min.cells = 10)
# Infer the cell-cell communication at a signaling pathway level
cellchat_ND <- computeCommunProbPathway(cellchat_ND)
# Calculate the aggregated cell-cell communication network
cellchat_ND <- aggregateNet(cellchat_ND)
cellchat_ND <- netAnalysis_computeCentrality(cellchat_ND, slot.name = "netP", thresh = 0.05)
#-------------------------------------------------------------------------------
# Compute the communication probability and infer cellular communication network
cellchat_HFD <- computeCommunProb(cellchat_HFD, type = "triMean")
cellchat_HFD <- filterCommunication(cellchat_HFD, min.cells = 10)
# Infer the cell-cell communication at a signaling pathway level
cellchat_HFD <- computeCommunProbPathway(cellchat_HFD)
# Calculate the aggregated cell-cell communication network
cellchat_HFD <- aggregateNet(cellchat_HFD)
cellchat_HFD <- netAnalysis_computeCentrality(cellchat_HFD, slot.name = "netP", thresh = 0.05)
#-------------------------------------------------------------------------------
# Compute the communication probability and infer cellular communication network
cellchat_TRF <- computeCommunProb(cellchat_TRF, type = "triMean")
cellchat_TRF <- filterCommunication(cellchat_TRF, min.cells = 10)
# Infer the cell-cell communication at a signaling pathway level
cellchat_TRF <- computeCommunProbPathway(cellchat_TRF)
# Calculate the aggregated cell-cell communication network
cellchat_TRF <- aggregateNet(cellchat_TRF)
cellchat_TRF <- netAnalysis_computeCentrality(cellchat_TRF, slot.name = "netP", thresh = 0.05)
#-------------------------------------------------------------------------------
#  Create cellchat object
object.list <- list(ND = cellchat_ND, HFD = cellchat_HFD, TRF = cellchat_TRF)
cellchat <- mergeCellChat(object.list, add.names = c("ND", "HFD", "TRF"))

cellchat <- qread("Fig.2F_mergeCellChat.qs")
#-------------------------------------------------------------------------------
## Identify and visualize the conserved and context-specific signaling pathways
#  Compare the overall information flow of each signaling pathway
gg1 <- rankNet(cellchat, comparison = c(1, 2), mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, comparison = c(1, 2), mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

gg3 <- rankNet(cellchat, comparison = c(2, 3), mode = "comparison", stacked = T, do.stat = TRUE)
gg4 <- rankNet(cellchat, comparison = c(2, 3), mode = "comparison", stacked = F, do.stat = TRUE)
gg3 + gg4

gg1+gg3

HFD <- as.data.frame(gg1[["data"]])
write.csv(HFD,"Fig.2F_CellChat_HFD.csv")
TRF <- as.data.frame(gg3[["data"]])
write.csv(TRF,"Fig.2F_CellChat_TRF.csv")
#-------------------------------------------------------------------------------
########
#Fig.2G
########
library(superheat)

Integration$CellType_Group <- paste(Integration$CellType, Integration$Group, sep = "_")
Idents(Integration) <- "CellType_Group"
table(Integration$CellType_Group)

#DEG identification
CellType <- FindMarkers(Integration, ident.1 = "CellType_TRF", ident.2 = "CellType_HFD", min.pct = 0.25, logfc.threshold = 0.5)
#CellType <- FindMarkers(Integration, ident.1 = "CellType_HFD", ident.2 = "CellType_ND", min.pct = 0.25, logfc.threshold = 0.5)
CellType <- subset(CellType, p_val_adj < 0.05)
CellType$Direction[CellType$avg_log2FC>0] <- "up"
CellType$Direction[CellType$avg_log2FC<0] <- "down"
#-------------------------------------------------------------------------------
##up-DEGs
DEGs <- read.csv("Fig.2G_DEGs.csv", header = T)

DEGs <- subset(DEGs,Direction=="up")
t <- as.data.frame(table(DEGs$GeneID))

t$Hep <- t$Var1 %in% DEGs$GeneID[DEGs$CellType=="Hep"]
t$MitoHigh <- t$Var1 %in% DEGs$GeneID[DEGs$CellType=="MitoHigh"]
t$IS <- t$Var1 %in% DEGs$GeneID[DEGs$CellType=="IS"]
t$EC <- t$Var1 %in% DEGs$GeneID[DEGs$CellType=="EC"]
t$HSC <- t$Var1 %in% DEGs$GeneID[DEGs$CellType=="HSC"]
t$Myof <- t$Var1 %in% DEGs$GeneID[DEGs$CellType=="Myof"]
t$Cho <- t$Var1 %in% DEGs$GeneID[DEGs$CellType=="Cho"]
t$Kup <- t$Var1 %in% DEGs$GeneID[DEGs$CellType=="Kup"]
t$Lym <- t$Var1 %in% DEGs$GeneID[DEGs$CellType=="Lym"]
t$Mono <- t$Var1 %in% DEGs$GeneID[DEGs$CellType=="Mono"]
t$DC <- t$Var1 %in% DEGs$GeneID[DEGs$CellType=="DC"]
write.csv(t,"Fig.2G_up-DEGs.csv")
#-------------------------------------------------------------------------------
##dw-DEGs
DEGs <- read.csv("Fig.2G_DEGs.csv", header = T)

DEGs <- subset(DEGs,Direction=="down")
t <- as.data.frame(table(DEGs$GeneID))

t$Hep <- t$Var1 %in% DEGs$GeneID[DEGs$CellType=="Hep"]
t$MitoHigh <- t$Var1 %in% DEGs$GeneID[DEGs$CellType=="MitoHigh"]
t$IS <- t$Var1 %in% DEGs$GeneID[DEGs$CellType=="IS"]
t$EC <- t$Var1 %in% DEGs$GeneID[DEGs$CellType=="EC"]
t$HSC <- t$Var1 %in% DEGs$GeneID[DEGs$CellType=="HSC"]
t$Myof <- t$Var1 %in% DEGs$GeneID[DEGs$CellType=="Myof"]
t$Cho <- t$Var1 %in% DEGs$GeneID[DEGs$CellType=="Cho"]
t$Kup <- t$Var1 %in% DEGs$GeneID[DEGs$CellType=="Kup"]
t$Lym <- t$Var1 %in% DEGs$GeneID[DEGs$CellType=="Lym"]
t$Mono <- t$Var1 %in% DEGs$GeneID[DEGs$CellType=="Mono"]
t$DC <- t$Var1 %in% DEGs$GeneID[DEGs$CellType=="DC"]
write.csv(t,"Fig.2G_dw-DEGs.csv")
#-------------------------------------------------------------------------------
DEGs <- read.csv("Fig.2G_up-DEGs.csv", header = T)

DEGs <- DEGs[1:34,]
rownames(DEGs) <- DEGs$Gene
DEGs <- DEGs[,1:11]
DEGs <- as.matrix(DEGs)
rownames(DEGs) <- factor(rownames(DEGs),levels = rownames(DEGs))

p1 <- as.ggplot(~superheat(DEGs, left.label.size = 0.2, left.label.text.size = 2,
          bottom.label.size = 0.1, bottom.label.text.size = 2,
          heat.col.scheme = "red", legend = F,bottom.label.text.angle = 90,
          grid.hline.col = "white",grid.hline.size = 1,
          grid.vline.col = "white",grid.vline.size = 2))
#-------------------------------------------------------------------------------
DEGs <- read.csv("Fig.2G_dw-DEGs.csv", header = T)

DEGs <- DEGs[1:40,]
rownames(DEGs) <- DEGs$Gene
DEGs <- DEGs[,1:11]
DEGs <- as.matrix(DEGs)
rownames(DEGs) <- factor(rownames(DEGs),levels = rownames(DEGs))

p2 <- as.ggplot(~superheat(DEGs, left.label.size = 0.2, left.label.text.size = 2,
          bottom.label.size = 0.1, bottom.label.text.size = 2,
          heat.col.scheme = "blue", legend = F,bottom.label.text.angle = 90,
          grid.hline.col = "white",grid.hline.size = 1,
          grid.vline.col = "white",grid.vline.size = 2))
#-------------------------------------------------------------------------------
p1|p2
