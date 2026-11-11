########
#Fig.S2
########
#clean environmental variables in case of unnecessary troubles
rm(list =ls()) 

#set working path
setwd("C:/Users/xiaoweilong/Desktop/Now!/Code/Figure 2/")
getwd()

#load needed R package
if(!require(Seurat))install.packages("Seurat")
if(!require(qs))install.packages("qs")
if(!require(qs))install.packages("CellChat")
#-------------------------------------------------------------------------------
#########
#Fig.S2A
#########
cellchat <- qread("Fig.2F_mergeCellChat.qs")

#Compare the total number of interactions and interaction strength
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3), measure = "count")
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3), measure = "weight")
gg1 + gg2
#-------------------------------------------------------------------------------
#########
#Fig.S2B
#########
#Compare the major sources and targets in 2D space
object.list <- qread("Fig.S2B_ObjectList.qs")

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link))
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}

# plots
patchwork::wrap_plots(plots = gg)

write.csv(gg[[1]][["data"]],"Fig.S2B_ND.csv")
write.csv(gg[[2]][["data"]],"Fig.S2B_HFD.csv")
write.csv(gg[[3]][["data"]],"Fig.S2B_TRF.csv")
#-------------------------------------------------------------------------------
#########
#Fig.S2C
#########
#Differential number of interactions or interaction strength in a greater details using a heatmap
p1 <- netVisual_heatmap(cellchat, comparison = c(1,2), color.heatmap = c("#2166ac","#b2182b"), measure = "count")
p3 <- netVisual_heatmap(cellchat, comparison = c(1,2), color.heatmap = c("#2166ac","#b2182b"), measure = "weight")
p2 <- netVisual_heatmap(cellchat, comparison = c(2,3), color.heatmap = c("#2166ac","#b2182b"), measure = "count")
p4 <- netVisual_heatmap(cellchat, comparison = c(2,3), color.heatmap = c("#2166ac","#b2182b"), measure = "weight")

p1+p2+p3+p4

write.csv(p1@matrix,"Fig.S2C_1.csv")
write.csv(p2@matrix,"Fig.S2C_2.csv")
write.csv(p3@matrix,"Fig.S2C_3.csv")
write.csv(p4@matrix,"Fig.S2C_4.csv")
#-------------------------------------------------------------------------------
#########
#Fig.S2D
#########
library(ggplot2)
library(tidyverse)
library(ggrepel)

df <- read.csv("Fig.2G_DEGs.csv", header = T) 
head(df)

df$label <- ifelse(df$avg_log2FC>0,"Upregulated genes", "Downregulated genes")
#-------------------------------------------------------------------------------
top10sig1 <- filter(df, CellType == "Hep")  %>% distinct(GeneID, .keep_all = T) %>% top_n(5, avg_log2FC)
top10sig2 <- filter(df, CellType == "MitoHigh")   %>% distinct(GeneID, .keep_all = T) %>% top_n(5, avg_log2FC)
top10sig3 <- filter(df, CellType == "IS")  %>% distinct(GeneID, .keep_all = T) %>% top_n(4, avg_log2FC)
top10sig4 <- filter(df, CellType == "HSC") %>% distinct(GeneID, .keep_all = T) %>% top_n(5, avg_log2FC)
top10sig5 <- filter(df, CellType == "Myof")  %>% distinct(GeneID, .keep_all = T) %>% top_n(5, avg_log2FC)
top10sig6 <- filter(df, CellType == "EC")  %>% distinct(GeneID, .keep_all = T) %>% top_n(5, avg_log2FC)
top10sig7 <- filter(df, CellType == "Cho")  %>% distinct(GeneID, .keep_all = T) %>% top_n(5, avg_log2FC)
top10sig8 <- filter(df, CellType == "Kup")  %>% distinct(GeneID, .keep_all = T) %>% top_n(5, avg_log2FC)
top10sig9 <- filter(df, CellType == "Lym")  %>% distinct(GeneID, .keep_all = T) %>% top_n(5, avg_log2FC)
top10sig10 <- filter(df, CellType == "Mono")  %>% distinct(GeneID, .keep_all = T) %>% top_n(5, avg_log2FC)
top10sig11 <- filter(df, CellType == "DC")  %>% distinct(GeneID, .keep_all = T) %>% top_n(2, avg_log2FC)
#-------------------------------------------------------------------------------
top10sig <- rbind(top10sig1,top10sig2,top10sig3,top10sig4,top10sig5,top10sig6,top10sig7,top10sig8,top10sig9,top10sig10,top10sig11)
#-------------------------------------------------------------------------------
top10sig1 <- filter(df, CellType == "Hep")  %>% distinct(GeneID, .keep_all = T) %>% top_n(-3, avg_log2FC)
top10sig2 <- filter(df, CellType == "MitoHigh")   %>% distinct(GeneID, .keep_all = T) %>% top_n(-3, avg_log2FC)
top10sig3 <- filter(df, CellType == "IS")  %>% distinct(GeneID, .keep_all = T) %>% top_n(-1, avg_log2FC)
top10sig4 <- filter(df, CellType == "HSC") %>% distinct(GeneID, .keep_all = T) %>% top_n(-3, avg_log2FC)
top10sig5 <- filter(df, CellType == "Myof")  %>% distinct(GeneID, .keep_all = T) %>% top_n(-3, avg_log2FC)
top10sig6 <- filter(df, CellType == "EC")  %>% distinct(GeneID, .keep_all = T) %>% top_n(-3, avg_log2FC)
top10sig7 <- filter(df, CellType == "Cho")  %>% distinct(GeneID, .keep_all = T) %>% top_n(-3, avg_log2FC)
top10sig8 <- filter(df, CellType == "Kup")  %>% distinct(GeneID, .keep_all = T) %>% top_n(-3, avg_log2FC)
top10sig9 <- filter(df, CellType == "Lym")  %>% distinct(GeneID, .keep_all = T) %>% top_n(-3, avg_log2FC)
top10sig10 <- filter(df, CellType == "Mono")  %>% distinct(GeneID, .keep_all = T) %>% top_n(-3, avg_log2FC)
top10sig11 <- filter(df, CellType == "DC")  %>% distinct(GeneID, .keep_all = T) %>% top_n(0, avg_log2FC)
#-------------------------------------------------------------------------------
top10sig <- rbind(top10sig,top10sig1,top10sig2,top10sig3,top10sig4,top10sig5,top10sig6,top10sig7,top10sig8,top10sig9,top10sig10,top10sig11)
#-------------------------------------------------------------------------------
df$size <- case_when(!(df$GeneID) %in% top10sig$GeneID ~ 1,
                     df$GeneID %in% top10sig$GeneID ~ 2)

dt <- filter(df, size==1)
head(dt)

dt$CellType <- factor(dt$CellType, levels = c("Hep","MitoHigh","IS","HSC","Myof","EC","Cho","Kup","Lym","Mono","DC"))
dt$label <- factor(dt$label, levels = c("Upregulated genes","Downregulated genes"))
#-------------------------------------------------------------------------------
p1 <- ggplot()+
  scale_color_manual(values = c("#E64B35","#00468B"))+
  geom_jitter(data = dt,
              aes(x = CellType, y = avg_log2FC, color = label),
              size = 0.85,width = 0.4)+
  
  geom_jitter(data = top10sig,
              aes(x = CellType, y = avg_log2FC, color = label),
              size = 1.2, width = 0.4)

p1
#-------------------------------------------------------------------------------
dfbar1 <- data.frame(x=c("Hep","MitoHigh","IS","HSC","Myof","EC","Cho","Kup","Lym","Mono","DC"),
                     y=c(3.4,3.6,3,5.1,4.8,4.8,4.8,4,4,4,2.6))
dfbar1$x <- factor(dfbar1$x, levels = c("Hep","MitoHigh","IS","HSC","Myof","EC","Cho","Kup","Lym","Mono","DC"))
dfbar2 <- data.frame(x=c("Hep","MitoHigh","IS","HSC","Myof","EC","Cho","Kup","Lym","Mono","DC"),
                     y=c(-3.5,-2.9,-3.2,-3.2,-5,-3.4,-3.2,-2.3,-3,-2.4,0))
dfbar2$x <- factor(dfbar2$x, levels = c("Hep","MitoHigh","IS","HSC","Myof","EC","Cho","Kup","Lym","Mono","DC"))

p2 <- ggplot()+
  geom_col(data = dfbar1,
           mapping = aes(x = x, y = y),
           fill = "#dcdcdc", alpha = 0.6)+
  
  geom_col(data = dfbar2,
           mapping = aes(x = x, y = y),
           fill = "#dcdcdc", alpha = 0.6)+
  
  scale_color_manual(values = c("#E64B35","#00468B"))+
  geom_jitter(data = dt,
              aes(x = CellType, y = avg_log2FC, color = label),
              size = 0.85,width = 0.4)+
  
  geom_jitter(data = top10sig,
              aes(x = CellType, y = avg_log2FC, color = label),
              size = 1.5, width = 0.4)


p2
#-------------------------------------------------------------------------------
dfcol <- data.frame(x=c("Hep","MitoHigh","IS","HSC","Myof","EC","Cho","Kup","Lym","Mono","DC"),y=0,
                    label=c("Hep","MitoHigh","IS","HSC","Myof","EC","Cho","Kup","Lym","Mono","DC"))
mycol <- c("purple","#E3AD68","#F6F5B4","#ACD45E","#edd05e","orange","#2F528F","#E54924","#9ECABE","lightgreen","skyblue")

p3 <- p2+geom_tile(data = dfcol, aes(x=x,y=y),
                   height=0.5, color="black", fill=mycol, alpha=0.6,show.legend = F)
p3
#-------------------------------------------------------------------------------
p4 <- p3+geom_text_repel(data = top10sig,
                         aes(x=CellType,y=avg_log2FC,label=GeneID),
                         force = 1,
                         arrow = arrow(length=unit(0.008,"npc"),
                                       type = "open", ends = "last"))
p4
#-------------------------------------------------------------------------------
p5 <- p4+labs(x="CellType",y="avg_log2FC")+
  geom_text(data=dfcol,aes(x=x,y=y,label=label),
            size=3,
            color="white")
p5
#-------------------------------------------------------------------------------
p6 <- p5+
  theme_minimal()+
  theme(axis.title = element_text(size = 11,color="black",face="bold"),
        axis.line.y = element_line(color="black",size=0.6),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.direction = "vertical",
        legend.justification = c(1,0),
        legend.text = element_text(size = 10)
  )
p6
#-------------------------------------------------------------------------------
#########
#Fig.S2E
#########
##Made in Excel.
#-------------------------------------------------------------------------------
#########
#Fig.S2F
#########
##Made in Metascape.
#-------------------------------------------------------------------------------
#########
#Fig.S2G
#########
##Made by Metascape.
#-------------------------------------------------------------------------------