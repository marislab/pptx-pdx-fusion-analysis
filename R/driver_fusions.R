################################################
# Author: Komal S Rathi
# Function: Driver fusion results for cBio and paper main
# Date: 12/18/2018
# literature fusions 
################################################

library(tidyr)
library(dplyr)

setwd('~/Projects/Maris-lab/PPTC_fusion_analysis/')

# clinical data
clin <- read.delim('data/2018-12-28-pdx-clinical-final-for-paper.txt', stringsAsFactors = F)
clin$RNA.human.bam.filename <- gsub('-R-human.bam$|-R_star_hg19_final.bam$|_star_hg19_final.bam$','', clin$RNA.human.bam.filename)
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
cm <- unique(cm[,c('Xenograft.ID','Key.translocations','Search_Type')])
cm <- cm[which(cm$Search_Type == "Fusion"),]
setdiff(cm$Xenograft.ID, clin$Model)
cm <- merge(cm, clin, by.x = 'Xenograft.ID', by.y = 'Model')
colnames(cm)[1] <- 'Model'
colnames(cm)[2] <- 'Fusion'
cm$Frame <- "in-frame"
cm$Fusion <- gsub('-','_',cm$Fusion)
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
all.callers$Fused_Genes <- gsub('--','_',all.callers$Fused_Genes)
if(length(setdiff(all.callers$Sample, clin$RNA.human.bam.filename)) == 0){
  all.callers <- merge(all.callers, clin, by.x = 'Sample', by.y = 'RNA.human.bam.filename')
  all.callers$Sample <- NULL
  all.callers$Fusion_Type <- ifelse(all.callers$Fusion_Type == "In-Frame", "in-frame", "frameshift")
}

# get literature genes
# these are both 5' and 3' genes but includes frameshift
lit.genes <- read.delim('data/2019-01-09-driver-fusions.txt', stringsAsFactors = F)
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
final$Entrez_Gene_Id[final$Hugo_Symbol == "C11orf95"] <- 65998
final$Entrez_Gene_Id[final$Hugo_Symbol == "CFAP300"] <- 85016
final <- unique(final[,c("Hugo_Symbol","Entrez_Gene_Id","Center","Tumor_Sample_Barcode","Fusion","DNA_support","RNA_support","Method","Frame","Histology.Broad","Histology.Broad.Label")])
final$Fusion <- gsub("_","-",final$Fusion)

# for cbio 
final[is.na(final)] <- ''
final <- final[-which(final$Fusion %in% c("SRP9-EPHX1")),]
intersect(final$Fusion, gsub('_','-',rts))
dim(unique(final[,c('Tumor_Sample_Barcode','Fusion')])) # 159

# driver-fusions
driver.fusions <- final[,c("Fusion","Tumor_Sample_Barcode","Histology.Broad.Label","Method","Frame")]
colnames(driver.fusions) <- c("Fused_Genes","Model","Histology.Broad.Label","Method","Fusion_Type")
driver.fusions$Fused_Genes <- sub('-','--',driver.fusions$Fused_Genes)
driver.fusions$Cytogenetics <- 'No'
driver.fusions[grep('Cytogenetics', driver.fusions$Method),'Cytogenetics'] <- "Yes"
write.table(driver.fusions, file = 'results/DriverFusions.txt', quote = F, sep = "\t", row.names = F)

# for cbio
final$Histology.Broad.Label <- NULL
write.table(final, file = '~/Projects/Maris-lab/PPTC/data/pedcbio/pptc/data_fusions.txt', quote = F, sep = "\t", row.names = F)

##### for JLH ##### 
# fusions (n = 95)
clin <- read.delim('data/2018-12-28-pdx-clinical-final-for-paper.txt', stringsAsFactors = F)
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
