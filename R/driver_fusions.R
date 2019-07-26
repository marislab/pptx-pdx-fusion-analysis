################################################
# Author: Komal S Rathi
# Function: Driver fusion results for cBio and paper main
# Date: 02/18/2019
# Step 2 
################################################

library(tidyr)
library(dplyr)
library(reshape2)
detach('package:plyr')

setwd('~/Projects/Maris-lab/PPTC_fusion_analysis/')

# clinical data
clin <- read.delim('data/2019-07-25-pdx-clinical-final-for-paper.txt', stringsAsFactors = F)
clin$RNA.human.bam.filename <- gsub('-R.human.bam$|-R-human.bam$|-R_star_hg19_final.bam$|_star_hg19_final.bam$','', clin$RNA.human.bam.filename)
clin <- clin[which(clin$RNA.Part.of.PPTC == "yes"),]
clin <- clin[,c("Model","Histology.Detailed","Histology.Broad","RNA.human.bam.filename")]
clin.ct <- plyr::count(clin$Histology.Broad)
clin.ct$freq <- paste0(clin.ct$x,' (n=',clin.ct$freq,')')
colnames(clin.ct)[2] <- "Histology.Broad.Label"
clin <- merge(clin, clin.ct, by.x = "Histology.Broad", by.y = "x")

# fusion data
# maria - key gene list
mc <- read.delim('data/PPTC_Fusion_Human-MC_v2.txt', stringsAsFactors = F)
mc <- mc[which(mc$key.gene.list == "yes"),]
mc <- unique(mc[,c('pairs','sample_RNA',"down_fusion_part_frame.shift_or_not")])
mc$sample_RNA <- gsub('-R$|','',mc$sample_RNA)
mc <- merge(mc, clin, by.x = 'sample_RNA', by.y = 'RNA.human.bam.filename')
mc$down_fusion_part_frame.shift_or_not[mc$down_fusion_part_frame.shift_or_not %in% c("inframe-shift")] <- "in-frame"
mc$down_fusion_part_frame.shift_or_not <- ifelse(mc$down_fusion_part_frame.shift_or_not %in% c("in-frame"), "in-frame", "frameshift") 
colnames(mc)[3] <- 'Frame'
colnames(mc)[2] <- 'Fusion'
mc2 <- mc
mc$Method <- 'deFuse'
mc2$Method <- 'SOAPFuse'
mc <- unique(rbind(mc, mc2))
rm(mc2)

# chelsea
cm <- read.delim('data/ALL_cytogenetics_translocations_CM_v2.txt', stringsAsFactors = F)
cm <- cm[which(cm$Search_Type == "Fusion"),]
cm$Key.translocations <- paste0(cm$GeneA,'--',cm$GeneB)
cm <- unique(cm[,c('Xenograft.ID','Key.translocations','Search_Type')])
setdiff(cm$Xenograft.ID, clin$Model) # 11
cm <- merge(cm, clin, by.x = 'Xenograft.ID', by.y = 'Model')
colnames(cm)[1] <- 'Model'
colnames(cm)[2] <- 'Fusion'
cm$Frame <- "in-frame"
cm$Fusion <- gsub('--','_',cm$Fusion)
cm$Method <- 'Cytogenetics'

# total (list #1)
common.cols <- intersect(colnames(cm), colnames(mc))
cm <- cm[,common.cols]
mc <- mc[,common.cols]
total <- unique(rbind(cm, mc))
rm(cm, mc)

# extract these gene fusions from other files
fc <- data.table::fread('data/fusioncatcher_PPTC_RNASeq_hg19_v2.txt', stringsAsFactors = F)
colnames(fc)[1:2] <- c("Gene1","Gene2")
colnames(fc)[9:10] <- c('Gene1_pos','Gene2_pos')
fc$Fusion_Type <- ifelse(fc$Predicted_effect == "in-frame",'In-Frame','Other')
fc$Caller <- 'FusionCatcher'
fc$Fused_Genes <- paste0(fc$Gene1,'--',fc$Gene2)
fc.total <- unique(fc[,c('Fused_Genes','Sample','Caller','Fusion_Type')])
fc.total$Fused_Genes <- gsub('@','',fc.total$Fused_Genes)

# star fusion
sf <- data.table::fread('data/starfusion_PPTC_RNASeq_hg19_v2.txt', stringsAsFactors = F)
sf$LeftBreakpoint <- gsub('^chr','',sf$LeftBreakpoint)
sf$RightBreakpoint <- gsub('^chr','',sf$RightBreakpoint)
colnames(sf)[c(6,8)] <- c('Gene1_pos','Gene2_pos')
sf$Fusion_Type <- ifelse(sf$PROT_FUSION_TYPE == "INFRAME",'In-Frame','Other')
sf$Caller <- 'STARFusion'
sf.total <- unique(sf[,c('Fused_Genes','Sample','Caller','Fusion_Type')])
sf.total$Fused_Genes <- gsub('@','',sf.total$Fused_Genes)

# soapfuse
soapfuse <- data.table::fread('data/SOAPfuse_PPTC_RNASeq_hg19_v2.txt', stringsAsFactors = F)
soapfuse$`down_fusion_part_frame-shift_or_not`[is.na(soapfuse$`down_fusion_part_frame-shift_or_not`)] <- ""
soapfuse$Fusion_Type <- ifelse(soapfuse$`down_fusion_part_frame-shift_or_not` == "inframe-shift",'In-Frame','Other')
soapfuse$Gene1_pos <- paste0(gsub('^chr','',soapfuse$up_chr),':',soapfuse$up_Genome_pos,':',soapfuse$up_strand)
soapfuse$Gene2_pos <- paste0(gsub('^chr','',soapfuse$dw_chr),':',soapfuse$dw_Genome_pos,':',soapfuse$dw_strand)
soapfuse$Caller <- 'SOAPFuse'
soapfuse$Fused_Genes <- paste0(soapfuse$up_gene,'--', soapfuse$dw_gene)
soapfuse.total <- unique(soapfuse[,c('Fused_Genes','Sample','Caller','Fusion_Type')])

# defuse
defuse <- data.table::fread('data/Defuse_PPTC_RNASeq_hg19_v2.txt', stringsAsFactors = F)
defuse <- defuse[-which(defuse$expression1 == 0 & defuse$expression2 == 0),]
defuse$Gene1_pos <- paste(defuse$gene_chromosome1,defuse$genomic_break_pos1, defuse$gene_strand1, sep = ":")
defuse$Gene2_pos <- paste(defuse$gene_chromosome2,defuse$genomic_break_pos2, defuse$gene_strand2, sep = ":")
defuse$Caller <- "deFuse"
defuse$Fused_Genes <- paste(defuse$gene_name1, defuse$gene_name2, sep = "--")
defuse$Fusion_Type <- ifelse(defuse$gene_location1 == "coding" & defuse$gene_location2 == "coding", "In-Frame", "Other")
defuse.total <- unique(defuse[,c('Fused_Genes','Sample','Caller','Fusion_Type')])

# get read-throughs
sf.rt <- unique(sf[grep('readthrough|neighbors|GTEx_Recurrent',sf$annots),'Fused_Genes'])
fc.rt <- unique(fc[grep('readthrough|adjacent|healthy',fc$Fusion_description),'Fused_Genes'])
soapfuse.rt <- unique(soapfuse[grep('INTRACHR-SS-OGO-0GAP', soapfuse$Fusion_Cat),'Fused_Genes'])
defuse.rt <- unique(defuse[which(defuse$adjacent == "Y" | defuse$read_through == "Y"),'Fused_Genes'])
rts <- unique(c(sf.rt$Fused_Genes, fc.rt$Fused_Genes, soapfuse.rt$Fused_Genes, defuse.rt$Fused_Genes))
rts.rev <- unique(unlist(lapply(strsplit(rts, '--'), FUN = function(x) paste0(x[2],'--',x[1]))))
rts <- unique(c(rts, rts.rev))
rts <- gsub('--','_',rts)

# per cancer how many unique fusions were called and how many are inframe
all.callers <- rbind(fc.total, sf.total, soapfuse.total, defuse.total)
all.callers$Sample <- gsub('-R-human$|-R$|-human$','',all.callers$Sample)
rm(fc.total, fc, fc.rt, sf.total, sf, sf.rt,
   soapfuse.total, soapfuse, defuse.total, defuse)
to.remove <- setdiff(all.callers$Sample, clin$RNA.human.bam.filename)
all.callers <- all.callers[-which(all.callers$Sample %in% to.remove),]
length(unique(all.callers$Sample)) # 244
all.callers$Fused_Genes <- gsub('--','_',all.callers$Fused_Genes)
if(length(setdiff(all.callers$Sample, clin$RNA.human.bam.filename)) == 0){
  all.callers <- merge(all.callers, clin, by.x = 'Sample', by.y = 'RNA.human.bam.filename')
  all.callers$Sample <- NULL
  all.callers$Fusion_Type <- ifelse(all.callers$Fusion_Type == "In-Frame", "in-frame", "frameshift")
}

# get literature genes
# these are both 5' and 3' genes but includes frameshift
lit.genes <- read.delim('data/2019-02-14-driver-fusions_v2.txt', stringsAsFactors = F)
lit.genes <- lit.genes[!is.na(lit.genes$FusionPartner),]
for(i in 1:nrow(lit.genes)){
  genes.to.search <- lit.genes[i,2]
  fusions.to.search <- lit.genes[i,3]
  genes.to.search <- unlist(strsplit(genes.to.search, ','))
  genes.to.search <- c(paste0('^',genes.to.search,'_'), paste0('_',genes.to.search,'$'))
  if(fusions.to.search == ""){
    print("no fusions to check")
  } else {
    fusions.to.search <- paste0('^',fusions.to.search,'$')
    genes.to.search <- c(genes.to.search, fusions.to.search)
  }
  genes.to.search <- paste0(genes.to.search, collapse = '|')
  hist.to.search <- lit.genes[i,1]
  getfusions <- all.callers[grep(genes.to.search, all.callers$Fused_Genes),]
  getfusions <- getfusions[which(getfusions$Histology.Detailed %in% hist.to.search),]
  getfusions <- unique(getfusions)
  if(nrow(getfusions) == 0){
    print(hist.to.search)
    print(genes.to.search)
  }
  if(i == 1){
    to.add <- getfusions
  } else {
    to.add <- rbind(to.add, getfusions)
  }
}
to.add <- unique(to.add)
colnames(to.add) <- c("Fusion", "Method", "Frame", "Histology.Broad", "Model", "Histology.Detailed", "Histology.Broad.Label")

# subset all.callers by literature genes (list #3)
all.callers <- all.callers[which(all.callers$Fused_Genes %in% total$Fusion),]
colnames(all.callers) <- colnames(to.add)

# merge the three lists
final <- rbind(all.callers, total, to.add)
final <- unique(final)

# remove read-throughs
final <- final[-which(final$Fusion %in% rts),]

####### separate low expressing fusions and fusions where gene expression not reported
load('data/pptc_rnaseq_hg38_matrix_v2.RData')
rna.mat$not_expressed <- apply(rna.mat[,2:ncol(rna.mat)], 1, FUN = function(x) all(x < 1))
df <- final
df <- cbind(colsplit(df$Fusion, pattern = '_', names = c("Gene1","Gene2")), df)
genes <- unique(c(df$Gene1, df$Gene2))
to.check <- setdiff(genes, rna.mat$gene_short_name) # 9
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
  hist <- df[i,'Histology.Broad']
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
separate.fusions <- df[(df$Gene1_not_expressed %in% c(TRUE,"Not Reported") & df$Gene2_not_expressed %in% c(TRUE,"Not Reported")),]
if(nrow(separate.fusions) > 0){
  print("Fusions to be separated")
  write.table(separate.fusions, file = 'results/Driver_Fusions_noExprReported.txt', quote = F, sep = "\t", row.names = F)
  df <- df[-which(df$Fusion %in% separate.fusions$Fusion),]
}
final <- final[which(final$Fusion %in% df$Fusion),]
####### separate low expressing fusions and fusions where gene expression not reported

# format for cBio
final <- final %>% 
  group_by(Fusion, Model, Frame, Histology.Broad, Histology.Broad.Label) %>%
  summarise(Method = toString(Method), n = n()) %>% unique() %>%
  as.data.frame()
ct <- plyr::count(final, c('Fusion','Model'))
ct <- ct[which(ct$freq > 1),]
final <- final[-which(final$Fusion %in% ct$Fusion & 
                        final$Model %in% ct$Model & 
                        final$Frame == "frameshift"),]
final$Hugo_Symbol <- gsub('_.*','', final$Fusion)
final$DNA_support <- 'unknown'
final$RNA_support <- 'yes'
colnames(final)[2] <- 'Tumor_Sample_Barcode'
final$Center <- 'BCM'
final$n <- NULL

# required columns
# Hugo_Symbol .
# Entrez_Gene_Id
# Center . 
# Tumor_Sample_Barcode . 
# Fusion .
# DNA support .
# RNA support .
# Method . 
# Frame . 

# now add entrez id
genes <- read.delim('~/Projects/Archive/cBio/genes.tsv')
genes <- unique(genes[,2:1])
final <- merge(final, genes, by.x = 'Hugo_Symbol',by.y = 'HUGO_GENE_SYMBOL', all.x = TRUE)
colnames(final)[11] <- 'Entrez_Gene_Id'
unique(final[is.na(final$Entrez_Gene_Id),'Hugo_Symbol'])
final$Entrez_Gene_Id[final$Hugo_Symbol == "CFAP300"] <- 85016
final <- unique(final[,c("Hugo_Symbol","Entrez_Gene_Id","Center","Tumor_Sample_Barcode","Fusion","DNA_support","RNA_support","Method","Frame","Histology.Broad","Histology.Broad.Label")])
for.driver.fusion <- final
final$Fusion <- gsub("_","-",final$Fusion)

# for cbio 
final[is.na(final)] <- ''
final <- final[-which(final$Fusion %in% c("SRP9-EPHX1")),]
intersect(final$Fusion, gsub('_','-',rts))
dim(unique(final[,c('Tumor_Sample_Barcode','Fusion')])) # 166
final$Histology.Broad.Label <- NULL
write.table(final, file = '~/Projects/Maris-lab/PPTC/data/pedcbio/pptc/data_fusions.txt', quote = F, sep = "\t", row.names = F)

# driver-fusions
driver.fusions <- for.driver.fusion[,c("Fusion","Tumor_Sample_Barcode","Histology.Broad.Label","Method","Frame")]
driver.fusions <- driver.fusions[-which(driver.fusions$Fusion %in% c("SRP9_EPHX1")),]
colnames(driver.fusions) <- c("Fused_Genes","Model","Histology.Broad.Label","Method","Fusion_Type")
driver.fusions$Fused_Genes <- sub('_','--',driver.fusions$Fused_Genes)
driver.fusions$Cytogenetics <- 'No'
driver.fusions[grep('Cytogenetics', driver.fusions$Method),'Cytogenetics'] <- "Yes"
write.table(driver.fusions, file = 'results/Driver_Fusions.txt', quote = F, sep = "\t", row.names = F)

# driver fusions collapsed
# fusions (n = 99)
clin <- read.delim('data/2019-07-25-pdx-clinical-final-for-paper.txt', stringsAsFactors = F)
clin$EXPRESSION <- clin$RNA.Part.of.PPTC
clin <- clin[which(clin$EXPRESSION == "yes"),]
ct <- plyr::count(clin$Histology.Detailed)
ct$label <- paste0(ct$x, ' (n=',ct$freq,')')
clin <- merge(clin, ct, by.x = 'Histology.Detailed', by.y = 'x')
clin <- clin[,c("Model","label")]

dat <- read.delim('~/Projects/Maris-lab/PPTC/data/pedcbio/pptc/data_fusions.txt')
dat <- dat[,c("Fusion","Tumor_Sample_Barcode")]
dat <- merge(dat, clin, by.x = 'Tumor_Sample_Barcode', by.y = 'Model')

dat <- dat %>% group_by(Fusion, label) %>% 
  summarise(toString(Tumor_Sample_Barcode)) %>%
  as.data.frame()
colnames(dat) <- c('Fused.Genes','Histology.Detailed','Models')
write.table(dat, file = 'results/DriverFusions_Collapsed.txt', quote = F, sep = "\t", row.names = F)
