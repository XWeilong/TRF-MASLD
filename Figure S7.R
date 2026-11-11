#clean environmental variables in case of unnecessary troubles
rm(list =ls()) 
gc()

#set working path
setwd("C:/Users/xiaoweilong/Desktop/Now!/Code/Figure 8/")
getwd()

library(Seurat)
library(qs)
if(!require(scRNAtoolVis))devtools::install_github('junjunlab/scRNAtoolVis')
#-------------------------------------------------------------------------------
###########
#Fig.S7A-D
###########
##Made in MetaboAnalyst.
#-------------------------------------------------------------------------------
#########
#Fig.S7E
#########
Integration <- qread("../Figure 3/Zonation.qs")


Genelist <- read.csv("Fig.8E_BCAA_genes.csv", header = T)
gene=unique(Genelist$BCAA.catabolism)

options(repr.plot.width=18, repr.plot.height=10)
p=DotPlot(Integration,features=rev(gene), dot.scale=8, scale.min = 0, scale.max = 100, group.by = "Group")+scale_y_discrete("")+scale_x_discrete("")+coord_flip()
p1=p+theme(panel.grid=element_blank())+theme_bw(base_size = 13)+#RotatedAxis()+
  scale_color_viridis_c()+labs(x=NULL,y=NULL)+theme(aspect.ratio=2)+ggtitle("")
plot(p1)
#-------------------------------------------------------------------------------
###########
#Fig.S7F-H
###########
#clean environmental variables in case of unnecessary troubles
rm(list =ls()) 

library(microeco)
library(magrittr)

set.seed(123)

otu <- read.csv("Fig.8F_OTU.csv",row.names = 1)
group <- read.csv("Fig.8F_Group.csv",header = T,row.names = 1)
class <- read.csv("Fig.8F_Class.csv",sep = ";",header = F)
colnames(class) <- c("Domain","Phylum","Class","Order","Family","Genus","Species")
rownames(class) <- rownames(otu)
#-------------------------------------------------------------------------------
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
#-------------------------------------------------------------------------------
#########
#Fig.S7F
#########
#beta diversity analysis
df$cal_betadiv()
class(df$beta_diversity)

t1 <- trans_beta$new(dataset = df, group = "Group", measure = "bray")
t1$cal_ordination(method = "PCoA")

class(t1$res_ordination)

# manova for each paired groups
t1$cal_manova(manova_all = F)
t1$res_manova

my_colors <- c(
  "ND" = "#D94A35",
  "HFD" = "#4DB1C9",
  "TRF" = "#109982")

t1[["sample_table"]][["Group"]] <- factor(
  t1[["sample_table"]][["Group"]],
  levels = c("ND", "HFD", "TRF"))

# plot the PCoA result with confidence ellipse
t1$plot_ordination(plot_color = "Group", plot_shape = "Group", 
                   plot_type = c("point", "ellipse"),color_values = my_colors)
#-------------------------------------------------------------------------------
#########
#Fig.S7G
#########
#alpha diversity analysis
df$cal_alphadiv(PD = FALSE)
t1 <- trans_alpha$new(dataset = df, group = "Group")

t1$cal_diff(method = "wilcox")

t1[["data_alpha"]]$Group <- factor(
  t1[["data_alpha"]]$Group,
  levels = c("ND", "HFD", "TRF"))

my_colors <- c(
  "ND" = "#D94A35",
  "HFD" = "#4DB1C9",
  "TRF" = "#109982")


t1$plot_alpha(measure = "Shannon", shape = "Group", color_values = my_colors,
              add_sig = F, add = "jitter", xtext_angle = 0)

write.csv(t1[["data_alpha"]],"Fig.S7G_Shannon_index.csv")
#-------------------------------------------------------------------------------
#########
#Fig.S7H
#########
#-------------------------------------------------------------------------------
t1 <- trans_abund$new(dataset = df, taxrank = "Phylum")

t1[["data_abund"]]$Group <- factor(t1[["data_abund"]]$Group, levels = c("ND","HFD","TRF"))

my_colors <- c(
  "ND" = "#D94A35",
  "HFD" = "#4DB1C9",
  "TRF" = "#109982")

t1$plot_box(group = "Group", xtext_angle = 30, xtext_size = 7, ytitle_size = 7, color_values = my_colors)#plot_flip = T
#t1$plot_bar(others_color = "grey70", facet = "Group", xtext_keep = FALSE, legend_text_italic = FALSE)

t1 <- trans_abund$new(dataset = df, taxrank = "Phylum")
data <- as.data.frame(t1[["data_abund"]])
write.csv(data,"Fig.S7H_Relative_Abundance.csv")
#-------------------------------------------------------------------------------