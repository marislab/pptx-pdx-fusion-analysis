################################################
# Author: Komal S Rathi
# Function: Annotate Fusions with Expression
# Date: 02/18/2019
# Step 5
################################################

library(reshape2)
library(tidyverse)
library(xlsx)
setwd('~/Projects/Maris-lab/PPTC_fusion_analysis/')

# annotate further on FPKM
clin <- read.delim('data/2019-07-25-pdx-clinical-final-for-paper.txt')
load('data/pptc_rnaseq_hg38_matrix_v2.RData')
rna.mat$not_expressed <- apply(rna.mat[,2:ncol(rna.mat)], 1, FUN = function(x) all(x < 1))
df <- read.delim('results/Filtered_Annotated_Fusions.txt', stringsAsFactors = F)
df <- cbind(colsplit(df$Fused_Genes, pattern = '--', names = c("Gene1","Gene2")), df)
genes <- unique(c(df$Gene1, df$Gene2))
to.check <- setdiff(genes, rna.mat$gene_short_name) # 39
rna.mat <- rna.mat[which(rna.mat$gene_short_name %in% genes),]
rna.mat <- melt(rna.mat)
rna.mat <- merge(rna.mat, clin[,c("Model","Histology.Broad")], by.x = 'variable', by.y = "Model")
rna.mat$Histology.Broad <- as.character(rna.mat$Histology.Broad)

# now add filter
df$Gene1_model <- NA
df$Gene1_hist_mean <- NA
df$Gene1_mean <- NA
df$Gene1_expr <- NA
df$Gene2_model <- NA
df$Gene2_hist_mean <- NA
df$Gene2_mean <- NA
df$Gene2_expr <- NA
df$Gene1_not_expressed <- NA
df$Gene2_not_expressed <- NA
for(i in 1:nrow(df)){
  print(i)
  genea <- df[i,'Gene1']
  geneb <- df[i,'Gene2']
  model <- df[i,'Model']
  hist <- df[i,'Histology.Broad.Label']
  hist <- gsub(' [(].*','', hist)
  genea.expr <- unique(rna.mat[which(rna.mat$gene_short_name == genea),'not_expressed'])
  geneb.expr <- unique(rna.mat[which(rna.mat$gene_short_name == geneb),'not_expressed'])
  genea.val <- rna.mat[which(rna.mat$variable %in% model & rna.mat$gene_short_name == genea),'value']
  geneb.val <- rna.mat[which(rna.mat$variable %in% model & rna.mat$gene_short_name == geneb),'value']
  genea.hist.mean <- mean(rna.mat[which(rna.mat$gene_short_name == genea & rna.mat$Histology.Broad == hist),'value'])
  geneb.hist.mean <- mean(rna.mat[which(rna.mat$gene_short_name == geneb & rna.mat$Histology.Broad == hist),'value'])
  genea.mean <- mean(rna.mat[which(rna.mat$gene_short_name == genea),'value'])
  geneb.mean <- mean(rna.mat[which(rna.mat$gene_short_name == geneb),'value'])
  df[i,'Gene1_not_expressed'] <- ifelse(length(genea.expr) == 0, NA, genea.expr)
  df[i,'Gene1_model'] <- ifelse(is.na(genea.val) || length(genea.val) == 0, NA, genea.val)
  df[i,'Gene1_hist_mean'] <- ifelse(is.na(genea.hist.mean) || length(genea.hist.mean) == 0, NA, genea.hist.mean)
  df[i,'Gene1_mean'] <- ifelse(is.na(genea.mean) || length(genea.mean) == 0, NA, genea.mean)
  df[i,'Gene1_expr'] <- ifelse(df[i,'Gene1_model'] < df[i,'Gene1_mean'],'Decreased',ifelse(df[i,'Gene1_model'] == df[i,'Gene1_mean'], 'Same','Increased'))
  df[i,'Gene2_not_expressed'] <- ifelse(length(geneb.expr) == 0, NA, geneb.expr)
  df[i,'Gene2_model'] <- ifelse(is.na(geneb.val) || length(geneb.val) == 0, NA, geneb.val)
  df[i,'Gene2_hist_mean'] <- ifelse(is.na(geneb.hist.mean) || length(geneb.hist.mean) == 0, NA, geneb.hist.mean)
  df[i,'Gene2_mean'] <- ifelse(is.na(geneb.mean) || length(geneb.mean) == 0, NA, geneb.mean)
  df[i,'Gene2_expr'] <- ifelse(df[i,'Gene2_model'] < df[i,'Gene2_mean'],'Decreased',ifelse(df[i,'Gene2_model'] == df[i,'Gene2_mean'], 'Same','Increased'))
}
df[is.na(df)] <- NA
df$Gene1_not_expressed[is.na(df$Gene1_not_expressed)] <- "Not Reported"
df$Gene2_not_expressed[is.na(df$Gene2_not_expressed)] <- "Not Reported"

# filters
write.table(df, 'results/Filtered_Annotated_Fusions_Expr', quote = F, sep = "\t", row.names = F)
write.xlsx(x = df, file = 'results/FusionTable.xlsx', sheetName = "Filtered_Annotated_Fusions_Expr", row.names = F, append = TRUE)
