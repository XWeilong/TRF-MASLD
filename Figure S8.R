########
#Fig.S8
########
#clean environmental variables in case of unnecessary troubles
rm(list =ls()) 
gc()

metabolites <- read.csv("Fig.S8_Differential Metabolites.csv",row.names = 1)
microbiota  <- read.csv("Fig.S8_Differential Species.csv",row.names = 1)
ID  <- read.csv("Fig.S8_ID.csv",row.names = 1, header = F)
#-------------------------------------------------------------------------------
colnames(metabolites) <- ID[2,1:27]
colnames(microbiota) <- ID[1,1:14]

cor_matrix <- matrix(NA, nrow = ncol(metabolites), ncol = ncol(microbiota))
p_matrix <- matrix(NA, nrow = ncol(metabolites), ncol = ncol(microbiota))
#-------------------------------------------------------------------------------
for (i in 1:ncol(metabolites)) {
  for (j in 1:ncol(microbiota)) {
    cor_test <- cor.test(metabolites[, i], microbiota[, j], method = "spearman")
    cor_matrix[i, j] <- cor_test$estimate
    p_matrix[i, j] <- cor_test$p.value
  }
}

rownames(cor_matrix) <- colnames(metabolites)
colnames(cor_matrix) <- colnames(microbiota)
rownames(p_matrix) <- colnames(metabolites)
colnames(p_matrix) <- colnames(microbiota)

write.csv(cor_matrix,"Fig.S8_cor_matrix.csv")
write.csv(p_matrix,"Fig.S8_p_matrix.csv")
#-------------------------------------------------------------------------------
significance_matrix <- ifelse(p_matrix < 0.05, "*", "")

display_matrix <- matrix(sprintf("%.2f", cor_matrix), nrow = nrow(cor_matrix), ncol = ncol(cor_matrix))
display_matrix[p_matrix >= 0.05] <- ""

pheatmap::pheatmap(cor_matrix, 
                   cluster_rows = TRUE, cluster_cols = TRUE, 
                   angle_col = 45, fontsize_col = 7, fontsize_row = 7,
                   display_numbers = display_matrix,
                   number_format = "%.2f",
                   cellnote = significance_matrix,
                   fontsize_number = 8,
                   color = colorRampPalette(c("blue", "white", "red"))(100),
                   legend_breaks = seq(-1, 1, 0.2),
                   legend_labels = seq(-1, 1, 0.2))
#-------------------------------------------------------------------------------
table_S <- cbind(metabolites,microbiota)
ggplot(table_S,aes(table_S$`N-(p-Coumaroyl)serotonin`,table_S$`Lawsonibacter sp009911595`))+
  geom_smooth(method="lm",
              se=T,
              color="red")+
  geom_point()+
  theme_classic()

myCor=cor.test(table_S$`N-(p-Coumaroyl)serotonin`,table_S$`Lawsonibacter sp009911595`,
               alternative = c("two.sided"),method="pearson")

myCor$estimate
myCor$p.value
#-------------------------------------------------------------------------------