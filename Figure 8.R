#clean environmental variables in case of unnecessary troubles
rm(list =ls()) 
gc()

#set working path
setwd("C:/Users/xiaoweilong/Desktop/Now!/Code/Figure 8/")
getwd()
#-------------------------------------------------------------------------------
##########
#Fig.8A-B
##########
##MetaboAnalyst.
#-------------------------------------------------------------------------------
########
#Fig.8C
########
##GraphPad Prism.
#-------------------------------------------------------------------------------
########
#Fig.8D
########

table_S <- read.csv("Fig.8D_Correlations.csv",row.names = 1)
Result <- cor(table_S,method = "pearson",use="pairwise.complete.obs")

library(Hmisc)
rcorr_result <- rcorr(as.matrix(table_S), type = "pearson")
p <- as.data.frame(rcorr_result[["P"]])

myCor=cor.test(table_S$BCAAs,table_S$TG_Liver,
               alternative = c("two.sided"),method="pearson")

myCor$estimate
myCor$p.value
#-------------------------------------------------------------------------------
library(ggplot2)

ggplot(table_S,aes(y=TG_Liver,x=BCAAs))+
  geom_smooth(method="lm",  
              se=T,         
              color="red")+
  geom_point()+           
  annotate(
    "text", x = min(table_S$BCAAs,na.rm = T), y = max(table_S$TG_Liver,na.rm = T)+5,
    label = paste0("R^2=",round(myCor$estimate,2),
                   ",P=",round(myCor$p.value,3)),
    size = 5,hjust=0,vjust=1)+
  theme_classic()
#-------------------------------------------------------------------------------
########
#Fig.8E
########

library(Seurat)
library(qs)
library(ggpubr)

Integration <- qread("../Figure 3/Zonation.qs")

Integration$Group <- factor(Integration$Group, levels = c("ND","HFD","TRF"))
#-------------------------------------------------------------------------------
#BCAAs
Genelist <- read.csv("Fig.8E_BCAA_genes.csv")
Genelist <- as.data.frame(Genelist$BCAA.catabolism, ncol = 1)
Genelist <- as.list(Genelist)

DefaultAssay(Integration) <- "SCT"
Integration <- AddModuleScore(Integration, assay = "SCT", features = Genelist, name = "BCAAs catabolic process")
colnames(Integration@meta.data)[ncol(Integration@meta.data)] <- "BCAAs catabolic process"

VlnPlot(Integration,features = "BCAAs catabolic process", group.by = "Group", sort = T, pt.size = 0) +NoLegend()
#-------------------------------------------------------------------------------
ggpubr::ggviolin(Integration@meta.data,x="Group", y="BCAAs catabolic process", width = 0.6, #x="group"
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
                 title = "")+NoLegend()
#-------------------------------------------------------------------------------
kruskal.test(Integration$`BCAAs catabolic process` ~ Integration$Group)
DescTools::DunnTest(Integration$`BCAAs catabolic process` ~ Integration$Group)

aggregate(Integration$`BCAAs catabolic process`, by = list(type = Integration$Group), FUN = mean)
aggregate(Integration$`BCAAs catabolic process`, by = list(type = Integration$Group), FUN = sd)
#-------------------------------------------------------------------------------
########
#Fig.8F
########

library(microeco)
library(magrittr)

set.seed(123)

otu <- read.csv("Fig.8F_OTU.csv",row.names = 1)
group <- read.csv("Fig.8F_Group.csv",header = T,row.names = 1)
class <- read.csv("Fig.8F_Class.csv",sep = ";",header = F)
colnames(class) <- c("Domain","Phylum","Class","Order","Family","Genus","Species")
rownames(class) <- rownames(otu)
#-------------------------------------------------------------------------------
#Only TRF and HFD group
group <- group[10:31,]
otu <- otu[,group$Sample]

# make the taxonomic information unified, very important
class %<>% tidy_taxonomy

#create a microtable object
df <- microtable$new(otu_table = otu,
                     tax_table = class,
                     sample_table = group)

df #microtable-class object

#make the OTU and sample information consistent
df$tidy_dataset()
print(df)

# use default parameters
df$cal_abund()

write.csv(df[["taxa_abund"]][["Species"]],"Fig.8F_Species_Abundance.csv")
#-------------------------------------------------------------------------------
t1 <- trans_diff$new(dataset = df, method = "lefse", 
                     taxa_level = "Species", group = "Group", alpha = 0.05)

# From v0.8.0, threshold is used for the LDA score selection.
t1$plot_diff_bar(threshold = 3, group_order = c("TRF", "HFD"))

write.csv(t1[["res_diff"]],"Fig.8F_LEfSe.csv")
#-------------------------------------------------------------------------------
########
#Fig.8G
########

#Microbiota-Metabolites correlation
table_S <- read.csv("Fig.8G_Correlations.csv",row.names = 1)
table_S <- table_S[,c(4,5)]
table_S <- as.data.frame(table_S)

myCor=cor.test(table_S$BCAAs,table_S$P.vulgatus,
               alternative = c("two.sided"),method="pearson")

myCor$estimate
myCor$p.value

ggplot(table_S,aes(y=BCAAs,x=P.vulgatus))+
  geom_smooth(method="lm",
              se=T,
              color="red")+
  geom_point()+
  annotate(
    "text", x = min(table_S$P.vulgatus,na.rm = T), y = max(table_S$BCAAs,na.rm = T)+100,
    label = paste0("R^2=",round(myCor$estimate,2),
                   ",P=",round(myCor$p.value,3)),
    size = 5,hjust=0,vjust=1)+
  theme_classic()
#-------------------------------------------------------------------------------
car::scatterplotMatrix(table_S,spread=F,lty.smooth=2,main="Scatter Plot Matrix")
#-------------------------------------------------------------------------------