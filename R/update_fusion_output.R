################################################
# Author: Komal S Rathi
# Function: Update gene symbols in Fusion output
# Date: 01/08/2019
# Step 1
################################################

setwd('~/Projects/Maris-lab/PPTC_fusion_analysis/')
library(dplyr)
library(tidyr)

################ read raw data from four callers ################ 
# fusion catcher
fc <- data.table::fread('data/fusioncatcher_PPTC_RNASeq_hg19.txt', stringsAsFactors = F)
fc$`Gene_1_symbol(5end_fusion_partner)` <- gsub('@','',fc$`Gene_1_symbol(5end_fusion_partner)`)
fc$`Gene_2_symbol(3end_fusion_partner)` <- gsub('@','',fc$`Gene_2_symbol(3end_fusion_partner)`)
for(i in 1:nrow(fc)){
  if(length(grep('^C',fc[i,1])) > 0){
    print(fc[i,1])
    fc[i,1] <- gsub('ORF','orf',fc[i,1])
  }
  if(length(grep('^C',fc[i,2])) > 0){
    print(fc[i,2])
    fc[i,2] <- gsub('ORF','orf',fc[i,2])
  }
}
fc[which(fc$`Gene_2_symbol(3end_fusion_partner)` == "TGIF2-C20ORF24"),2] <- "TGIF2-C20orf24"

# star fusion
sf <- data.table::fread('data/starfusion_PPTC_RNASeq_hg19.txt', stringsAsFactors = F)
sf$LeftGene <- gsub('\\^.*','',sf$LeftGene)
sf$RightGene <- gsub('\\^.*','',sf$RightGene)
sf$LeftGene <- gsub('@','',sf$LeftGene)
sf$RightGene <- gsub('@','',sf$RightGene)

# soapfuse
soapfuse <- data.table::fread('data/SOAPfuse_PPTC_RNASeq_hg19.txt', stringsAsFactors = F)

# defuse
defuse <- data.table::fread('data/Defuse_PPTC_RNASeq_hg19.txt', stringsAsFactors = F)

# mc
mc <- read.delim('data/PPTC_Fusion_Human-MC.txt', stringsAsFactors = F)
mc$GeneA <- gsub('_.*','',mc$pairs)
mc$GeneB <- gsub('.*_','',mc$pairs)

# cm
cm <- read.delim('data/ALL_cytogenetics_translocations_CM.txt', stringsAsFactors = F)
cm$Key.translocations <- gsub('ORF','orf',cm$Key.translocations)
cm$GeneA <- gsub('-.*','',cm$Key.translocations)
cm$GeneB <- gsub('.*-','',cm$Key.translocations)

# literature genes
lit.genes <- read.delim('data/2019-02-14-driver-fusions.txt', stringsAsFactors = F)
lit.genes <- lit.genes %>%
  group_by(Histology) %>%
  mutate(FusionPartner = strsplit(FusionPartner, ",")) %>%
  unnest(FusionPartner)

# fix gene names from old to new
dat <- read.delim('data/2019-02-14-Hugo-Symbols-approved.txt', stringsAsFactors = F)
dat <- unique(dat[,c('Approved.symbol','Previous.symbols')])
genes.to.remove <- intersect(dat$Approved.symbol, dat$Previous.symbols)
dat <- dat[-which(dat$Previous.symbols %in% genes.to.remove),]
dat <- dat[which(dat$Previous.symbols != ""),]
dat$Previous.symbols <- gsub(',.*','',dat$Previous.symbols)
dat <- dat %>%
  mutate(Previous.symbols = strsplit(Previous.symbols, ",")) %>%
  unnest(Previous.symbols)
dat <- dat %>% 
  mutate(Previous.symbols = trimws(Previous.symbols, which = "both"),
         Approved.symbol = trimws(Approved.symbol, which = "both")) %>%
  mutate(Previous.symbols = gsub('@','',Previous.symbols),
         Approved.symbol = gsub('@','',Approved.symbol))
dat <- dat[which(dat$Approved.symbol != dat$Previous.symbols),]
tmp <- plyr::count(dat$Previous.symbols)
tmp <- tmp[which(tmp$freq == 1),]
dat <- dat[which(dat$Previous.symbols %in% tmp$x),]
rm(tmp, genes.to.remove)

# replace in driver fusions
for(i in 1:nrow(lit.genes)){
  print(i)
  if(lit.genes[i,3] %in% dat$Previous.symbols == TRUE){
    new.gene <- dat[dat$Previous.symbols %in% lit.genes[i,3],'Approved.symbol']
    lit.genes[i,3] <- new.gene
  }
}
lit.genes <- lit.genes %>% 
  group_by(Histology, Fusion) %>%
  summarise(FusionPartner = paste(FusionPartner, collapse = ",")) %>%
  as.data.frame()
lit.genes <- lit.genes[,c(1,3,2)]
write.table(lit.genes, file = 'data/2019-02-14-driver-fusions_v2.txt', quote = F, sep = "\t", row.names = F)

# replace in fusion catcher
for(i in 1:nrow(fc)){
  print(i)
  if(fc[i,1] %in% dat$Previous.symbols == TRUE){
    new.gene <- dat[dat$Previous.symbols %in% fc[i,1],'Approved.symbol'] 
    fc[i,1] <- new.gene
  }
  if(fc[i,2] %in% dat$Previous.symbols == TRUE){
    new.gene <- dat[dat$Previous.symbols %in% fc[i,2],'Approved.symbol'] 
    fc[i,2] <- new.gene
  }
}
write.table(fc, file = 'data/fusioncatcher_PPTC_RNASeq_hg19_v2.txt', quote = F, sep = "\t", row.names = F)


# replace in star fusion
for(i in 1:nrow(sf)){
  print(i)
  if(sf[i,5] %in% dat$Previous.symbols == TRUE){
    new.gene <- dat[dat$Previous.symbols %in% sf[i,5],'Approved.symbol'] 
    sf[i,5] <- new.gene
  }
  if(sf[i,7] %in% dat$Previous.symbols == TRUE){
    new.gene <- dat[dat$Previous.symbols %in% sf[i,7],'Approved.symbol'] 
    sf[i,7] <- new.gene
  }
}
sf$Fused_Genes <- paste0(sf$LeftGene, '--', sf$RightGene)
write.table(sf, file = 'data/starfusion_PPTC_RNASeq_hg19_v2.txt', quote = F, sep = "\t", row.names = F)


# replace in defuse
for(i in 1:nrow(defuse)){
  print(i)
  if(defuse[i,31] %in% dat$Previous.symbols == TRUE){
    new.gene <- dat[dat$Previous.symbols %in% defuse[i,31],'Approved.symbol'] 
    defuse[i,31] <- new.gene
  }
  if(defuse[i,32] %in% dat$Previous.symbols == TRUE){
    new.gene <- dat[dat$Previous.symbols %in% defuse[i,32],'Approved.symbol'] 
    defuse[i,32] <- new.gene
  }
}
write.table(defuse, file = 'data/Defuse_PPTC_RNASeq_hg19_v2.txt', quote = F, sep = "\t", row.names = F)

# replace in soapfuse
for(i in 1:nrow(soapfuse)){
  print(i)
  if(soapfuse[i,1] %in% dat$Previous.symbols == TRUE){
    new.gene <- dat[dat$Previous.symbols %in% soapfuse[i,1],'Approved.symbol'] 
    soapfuse[i,1] <- new.gene
  }
  if(soapfuse[i,6] %in% dat$Previous.symbols == TRUE){
    new.gene <- dat[dat$Previous.symbols %in% soapfuse[i,6],'Approved.symbol'] 
    soapfuse[i,6] <- new.gene
  }
}
write.table(soapfuse, file = 'data/SOAPfuse_PPTC_RNASeq_hg19_v2.txt', quote = F, sep = "\t", row.names = F)

# replace in mc
for(i in 1:nrow(mc)){
  print(i)
  if(mc[i,17] %in% dat$Previous.symbols == TRUE){
    new.gene <- dat[dat$Previous.symbols %in% mc[i,17],'Approved.symbol'] 
    mc[i,17] <- new.gene
  }
  if(mc[i,18] %in% dat$Previous.symbols == TRUE){
    new.gene <- dat[dat$Previous.symbols %in% mc[i,18],'Approved.symbol'] 
    mc[i,18] <- new.gene
  }
}
mc$pairs <- paste0(mc$GeneA,'_',mc$GeneB)
mc$GeneA <- NULL
mc$GeneB <- NULL
write.table(mc, file = 'data/PPTC_Fusion_Human-MC_v2.txt', quote = F, sep = "\t", row.names = F)

# replace in cm
for(i in 1:nrow(cm)){
  if(cm[i,'Search_Type'] == "Fusion"){
    print(i)
    if(cm[i,9] %in% dat$Previous.symbols == TRUE){
      new.gene <- dat[dat$Previous.symbols %in% cm[i,9],'Approved.symbol'] 
      cm[i,9] <- new.gene
    }
    if(cm[i,10] %in% dat$Previous.symbols == TRUE){
      new.gene <- dat[dat$Previous.symbols %in% cm[i,10],'Approved.symbol'] 
      cm[i,10] <- new.gene
    }
    cm[i,3] <- paste0(cm[i,9],"-",cm[i,10])
  }
}
write.table(cm, file = 'data/ALL_cytogenetics_translocations_CM_v2.txt', quote = F, sep = "\t", row.names = F)

# replace symbols in rna matrix
load('../PPTC/data/pedcbio/2019-02-14-PPTC_FPKM_matrix_withModelID-244.rda')
rna.mat$gene_short_name <- as.character(rna.mat$gene_short_name)
for(i in 1:nrow(rna.mat)){
  print(i)
  if(rna.mat[i,1] %in% dat$Previous.symbols == TRUE){
    new.gene <- dat[dat$Previous.symbols %in% rna.mat[i,1],'Approved.symbol'] 
    rna.mat[i,1] <- new.gene
  }
}
save(rna.mat, file = 'data/2019-02-14-PPTC_FPKM_matrix_withModelID-244_v2.rda')

# replace symbols in rna matrix (hg38)
load('data/pptc_rnaseq_hg38_matrix.RData')
rna.mat$gene_short_name <- as.character(rna.mat$gene_short_name)
for(i in 1:nrow(rna.mat)){
  print(i)
  if(rna.mat[i,1] %in% dat$Previous.symbols == TRUE){
    new.gene <- dat[dat$Previous.symbols %in% rna.mat[i,1],'Approved.symbol'] 
    rna.mat[i,1] <- new.gene
  }
}
save(rna.mat, file = 'data/pptc_rnaseq_hg38_matrix_v2.RData')
