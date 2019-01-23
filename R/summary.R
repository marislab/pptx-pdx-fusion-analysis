################################################
# Author: Komal S Rathi
# Function: Summarize fusion results for SV figure
# Date: 12/18/2018
################################################

library(ggplot2)
library(reshape2)
library(dplyr)
library(RColorBrewer)
library(tidyr)

setwd('~/Projects/Maris-lab/PPTC_fusion_analysis/')
source('~/Projects/Maris-lab/PPTC/R/pubTheme.R')
source('R/themes.R')

################ read raw data from four callers ################ 
# fusion catcher
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

# defuse (only defuse has expression values so we can use that filter)
defuse <- data.table::fread('data/Defuse_PPTC_RNASeq_hg19_v2.txt', stringsAsFactors = F)
defuse <- defuse[-which(defuse$expression1 == 0 & defuse$expression2 == 0),]
defuse$Gene1_pos <- paste(defuse$gene_chromosome1,defuse$genomic_break_pos1, defuse$gene_strand1, sep = ":")
defuse$Gene2_pos <- paste(defuse$gene_chromosome2,defuse$genomic_break_pos2, defuse$gene_strand2, sep = ":")
defuse$Caller <- "deFuse"
defuse$Fused_Genes <- paste(defuse$gene_name1, defuse$gene_name2, sep = "--")
defuse$Fusion_Type <- ifelse(defuse$gene_location1 == "coding" & defuse$gene_location2 == "coding", "In-Frame", "Other")
defuse.total <- unique(defuse[,c('Fused_Genes','Sample','Caller','Fusion_Type')])

# merge all data
all.callers <- rbind(fc.total, sf.total, soapfuse.total, defuse.total)
all.callers$Sample <- gsub('-R-human$|-R$|','',all.callers$Sample)

# this is temporary
# all.callers$GeneA <- gsub('--.*','',all.callers$Fused_Genes)
# all.callers$GeneB <- gsub('.*--','',all.callers$Fused_Genes)

# fix gene names from old to new
# dat <- read.delim('data/2019-01-03-Hugo-Symbols-edited.txt', stringsAsFactors = F)
# dat <- unique(dat[,c('Approved.symbol','Previous.symbols')])
# dat <- dat[which(dat$Previous.symbols != ""),]
# dat <- dat %>%
#   mutate(Previous.symbols = strsplit(Previous.symbols, ",")) %>%
#   unnest(Previous.symbols)
# dat <- dat %>% 
#   mutate(Previous.symbols = trimws(Previous.symbols, which= "both"),
#          Approved.symbol = trimws(Approved.symbol, which = "both")) %>%
#   mutate(Previous.symbols = gsub('@','',Previous.symbols),
#          Approved.symbol = gsub('@','',Approved.symbol))
# dat <- dat[which(dat$Approved.symbol != dat$Previous.symbols),]
# 
# tmp <- plyr::count(dat$Previous.symbols)
# tmp <- tmp[which(tmp$freq == 1),]
# dat <- dat[which(dat$Previous.symbols %in% tmp$x),]
# rm(tmp)
# all.callers$A.map <- ifelse(all.callers$GeneA %in% dat$Previous.symbols, 'Yes', 'No')
# all.callers$B.map <- ifelse(all.callers$GeneB %in% dat$Previous.symbols, 'Yes', 'No')

# for(i in 1:nrow(all.callers)){
#   print(i)
#   if(all.callers[i,'A.map'] == "Yes"){
#     new.gene <- dat[dat$Previous.symbols %in% all.callers[i,'GeneA'],'Approved.symbol'] 
#     all.callers[i,'GeneA'] <- new.gene
#   }
#   if(all.callers[i,'B.map'] == "Yes"){
#     new.gene <- dat[dat$Previous.symbols %in% all.callers[i,'GeneB'],'Approved.symbol'] 
#     all.callers[i,'GeneB'] <- new.gene
#   }
# }
# all.callers$Fused_Genes <- paste0(all.callers$GeneA,'--',all.callers$GeneB)
# all.callers <- unique(all.callers[,1:4])

# remove @ genes
# all.callers <- all.callers[grep('@', all.callers$Fused_Genes, invert = T),]

# remove the samples that we don't need
clin <- read.delim('data/2018-12-28-pdx-clinical-final-for-paper.txt', stringsAsFactors = F)
clin <- clin[which(clin$RNA.Part.of.PPTC == "yes"),]
clin <- clin[,c('RNA.human.bam.filename','Model','Histology.Broad')]
clin$RNA.human.bam.filename <- gsub('-R-human.bam$|-R_star_hg19_final.bam$|_star_hg19_final.bam$|-human.bam$','', clin$RNA.human.bam.filename)

# histology labels
hist.dt.ct <- unique(clin[,c('Histology.Broad','Model')])
hist.dt.ct <- plyr::count(hist.dt.ct, 'Histology.Broad')
hist.dt.ct$freq <- paste0(hist.dt.ct$Histology.Broad,' (n=', hist.dt.ct$freq, ')')
colnames(hist.dt.ct)[2] <- 'Histology.Broad.Label'

# no fusions reported for NCH-RBD-1 so n = 247 here
to.remove <- setdiff(all.callers$Sample, clin$RNA.human.bam.filename)
if(length(to.remove) > 0){
  all.callers <- all.callers[-which(all.callers$Sample %in% to.remove),]
}
length(unique(all.callers$Sample))

# add Model & Histology at the end
all.callers <- merge(all.callers, clin, by.x = 'Sample', by.y = 'RNA.human.bam.filename')
all.callers <- unique(all.callers[,c('Model','Fused_Genes','Caller','Fusion_Type','Histology.Broad')])
to.add <- all.callers[,c('Model','Fused_Genes','Caller','Fusion_Type')] # for final annotation

# Gene fusion should be in-frame
# Called by at least 2 callers
all.callers.summary <- all.callers %>% 
  filter(Fusion_Type != "Other") %>%
  group_by(Fused_Genes, Model, Histology.Broad) %>% 
  unique() %>%
  mutate(Caller = toString(Caller), caller.count = n()) %>%
  filter(caller.count >= 2) %>% 
  select(-caller.count, -Caller, -Fusion_Type) %>%
  unique() %>%
  as.data.frame()

# or found in at least 2 samples of the same histology (n = 545)
sample.count <- all.callers %>% 
  filter(Fusion_Type != "Other") %>%
  group_by(Fused_Genes, Histology.Broad) %>% 
  unique() %>%
  mutate(sample.count = n()) %>%
  filter(sample.count > 1) %>%
  select(-Caller, -sample.count, -Fusion_Type) %>%
  unique() %>%
  as.data.frame() 
length(unique(sample.count$Fused_Genes))

# or 3' or 5' gene recurrently fused within a histology (>= 5 genes)
rec <- cbind(all.callers, colsplit(all.callers$Fused_Genes, pattern = '--', names = c("GeneA","GeneB")))
rec2 <- rec %>% group_by(Histology.Broad) %>% 
  select(Histology.Broad,GeneA,GeneB) %>% 
  unique() %>% group_by(Histology.Broad, GeneA) %>% 
  summarise(GeneA.ct = n()) %>%
  filter(GeneA.ct >= 5) %>% as.data.frame()
rec3 <- rec %>% group_by(Histology.Broad) %>% 
  select(Histology.Broad,GeneA,GeneB) %>% 
  unique() %>% group_by(Histology.Broad, GeneB) %>% 
  summarise(GeneB.ct = n()) %>%
  filter(GeneB.ct >= 5) %>% as.data.frame()
rec2 <- merge(rec2, rec, by = c('GeneA','Histology.Broad'))
rec3 <- merge(rec3, rec, by = c('GeneB','Histology.Broad'))
rec2 <- unique(rec2[,c("Model","Fused_Genes","Histology.Broad")])
rec3 <- unique(rec3[,c("Model","Fused_Genes","Histology.Broad")])
res <- unique(rbind(rec2, rec3))

# merge these (n = 6287)
total <- unique(rbind(all.callers.summary, sample.count, res))
nrow(total)

# remove fusions that are in > 1 histology (n = 1183)
hist.count <- total %>% 
  select(Fused_Genes, Histology.Broad) %>%
  unique() %>%
  group_by(Fused_Genes) %>%
  summarise(hist.count = n()) %>%
  filter(hist.count == 1)
total <- total[which(total$Fused_Genes %in% hist.count$Fused_Genes),]
length(unique(total$Fused_Genes))

# remove read-throughs or neighbouring genes (n = 944)
sf.rt <- unique(sf[grep('readthrough|neighbors|GTEx_Recurrent',sf$annots),'Fused_Genes'])
fc.rt <- unique(fc[grep('readthrough|adjacent|healthy',fc$Fusion_description),'Fused_Genes'])
soapfuse.rt <- unique(soapfuse[grep('INTRACHR-SS-OGO-0GAP', soapfuse$Fusion_Cat),'Fused_Genes'])
defuse.rt <- unique(defuse[which(defuse$adjacent == "Y" | defuse$read_through == "Y"),'Fused_Genes'])
rts <- unique(c(sf.rt$Fused_Genes, fc.rt$Fused_Genes, soapfuse.rt$Fused_Genes, defuse.rt$Fused_Genes))
rts.rev <- unique(unlist(lapply(strsplit(rts, '--'), FUN = function(x) paste0(x[2],'--',x[1]))))
rts <- unique(c(rts, rts.rev))
total <- total[-which(total$Fused_Genes %in% rts),]
length(unique(total$Fused_Genes))

# add back MLL fusions by Chelsea, 
# Pubmed fusions by Maria and 
# literature fusions by Jo Lynne (generated using driver_fusions.R)
drivers <- read.delim('~/Projects/Maris-lab/PPTC/data/pedcbio/pptc/data_fusions.txt', stringsAsFactors = F)
drivers <- unique(drivers[,c('Fusion','Tumor_Sample_Barcode')])
drivers <- merge(drivers, clin, by.x = 'Tumor_Sample_Barcode', by.y = 'Model')
drivers$RNA.human.bam.filename <- NULL
drivers$Fusion <- gsub('-','--',drivers$Fusion)
drivers$Fusion <- gsub('RP11--', 'RP11-', drivers$Fusion)
colnames(drivers) <- colnames(total)
# check for read-throughs (none should be present)
length(intersect(drivers$Fused_Genes, rts))
# check how many we missed or captured
length(unique(setdiff(drivers$Fused_Genes, total$Fused_Genes))) # 37/94 unique fusions were missed
total <- unique(rbind(total, drivers))
length(unique(total$Fused_Genes)) # (n = 981)

# split in gene a and gene b
total <- cbind(total, colsplit(total$Fused_Genes, '--', names = c('geneA', 'geneB')))

# fusions in normal cells (none - this is good)
normals.fc <- unique(as.character(fc[grep('gtex', fc$Fusion_description),'Fused_Genes']$Fused_Genes))
normals.sf <- unique(as.character(sf[grep('GTEx', sf$annots),'Fused_Genes']$Fused_Genes))
normals <- unique(c(normals.fc, normals.sf))
if(length(intersect(total$Fused_Genes, normals)) > 0){
  total <- total[-which(total$Fused_Genes %in% normals),]
} else {
  total <- total
}
length(unique(total$Fused_Genes)) # 981

# add Cytogenetics info
cyto <- read.delim('results/DriverFusions.txt', stringsAsFactors = F)
cyto <- unique(cyto[which(cyto$Cytogenetics == "Yes"),1:2])
cyto$Cytogenetics <- "Yes"

# excel tab 1
extab1 <- merge(total, to.add, by = c('Model','Fused_Genes'))
extab1 <- merge(extab1, hist.dt.ct, by = 'Histology.Broad')
extab1 <- extab1[,c("Model","Fused_Genes","Histology.Broad.Label","Caller","Fusion_Type")]
extab1 <- extab1[order(extab1$Fused_Genes, extab1$Model),]
extab1 <- merge(extab1, cyto, all.x = T)
extab1$Cytogenetics[is.na(extab1$Cytogenetics)] <- "No"
extab1 <- extab1[,c("Fused_Genes","Model","Histology.Broad.Label","Caller","Fusion_Type","Cytogenetics")]

# annotate using COSMIC 
for.table <- total
cosmic <- read.csv('~/Projects/Maris-lab/NBL-Fusion-Analysis/data/Cosmic_gene_census.csv', stringsAsFactors = F)
genes.in.cosmic <- unlist(strsplit(cosmic$Synonyms, split = ","))
genes.in.cosmic <- unique(c(cosmic$Gene.Symbol, genes.in.cosmic))
rm(cosmic)
for.table$Cosmic <- ifelse(for.table$geneA %in% genes.in.cosmic | for.table$geneB %in% genes.in.cosmic, 'Yes', 'No')

# annotate using kinases
kinase <- read.delim('~/Projects/Maris-lab/NBL-Fusion-Analysis/data/kinases.txt')
for.table$Kinase <- ifelse(for.table$geneA %in% kinase$gene_symbol | for.table$geneB %in% kinase$gene_symbol, 'Yes', 'No')

# annotate using transcription factors
tf <- read.delim('~/Projects/Maris-lab/NBL-Fusion-Analysis/data/TRANSFAC_TF.txt', header = F)
for.table$TF <- ifelse(for.table$geneA %in% tf$V1 | for.table$geneB %in% tf$V1, 'Yes', 'No')

# n = 306
for.table <- unique(for.table[which(for.table$Cosmic == "Yes" | for.table$Kinase == "Yes" | for.table$TF == "Yes"),])
setdiff(drivers$Fused_Genes, for.table$Fused_Genes) # 5 fusions are missing
for.table$geneA <- NULL
for.table$geneB <- NULL

# excel tab 2
extab2 <- merge(for.table, hist.dt.ct, by = 'Histology.Broad')
extab2 <- extab2[,c("Model","Fused_Genes","Histology.Broad.Label","Cosmic","TF","Kinase")]
extab2 <- extab2 %>% 
  group_by(Fused_Genes, Histology.Broad.Label, Cosmic, TF, Kinase) %>% 
  summarise(Model = toString(Model), Models.Found.In = n()) %>% 
  unique() %>% as.data.frame()
extab2 <- unique(extab2[,c("Fused_Genes","Cosmic","TF","Kinase")])
extab1 <- merge(extab1, extab2, by = 'Fused_Genes', all.x = T)
extab1[is.na(extab1)] <- "No"
# 2884
colnames(extab1)[4] <- 'Method'
write.table(extab1, file = 'results/Filtered_Annotated_Fusions.txt', quote = F, sep = "\t", row.names = F)

#########
# plot for fusions only
#########
total$geneA <- NULL
total$geneB <- NULL

clin <- read.delim('data/2018-12-28-pdx-clinical-final-for-paper.txt', stringsAsFactors = F)
clin <- clin[which(clin$RNA.Part.of.PPTC == "yes"),]
clin <- clin[,c('RNA.human.bam.filename','Model','Histology.Detailed')]
clin$RNA.human.bam.filename <- gsub('-R-human.bam$|-R_star_hg19_final.bam$|_star_hg19_final.bam$','', clin$RNA.human.bam.filename)

# histology detailed labels
hist.ct <- unique(clin[,c('Histology.Detailed','Model')])
hist.ct <- plyr::count(hist.ct, 'Histology.Detailed')
hist.ct$freq <- paste0(hist.ct$Histology.Detailed,' (n=', hist.ct$freq, ')')
colnames(hist.ct)[2] <- 'Histology.Detailed.Label'

final <- merge(total, clin[,c('Model','Histology.Detailed')], by = 'Model')
final <- merge(final, hist.ct, by = 'Histology.Detailed')
final <- final %>% group_by(Model, Histology.Detailed) %>% summarise(value = n()) 
final <- final %>% group_by(Histology.Detailed) %>% mutate(median = median(value)) %>% as.data.frame()
to.include <- setdiff(clin$Histology.Detailed, final$Histology.Detailed)
final <- rbind(final, data.frame(Model = c(rep(NA,length(to.include))), 
                                 Histology.Detailed = to.include, 
                                 value = c(rep(0,length(to.include))), 
                                 median = c(rep(0,length(to.include)))))
final$Histology.Detailed <- reorder(final$Histology.Detailed, final$median)
write.table(final, file = 'results/FusionPlot_rawdata.txt', quote = F, sep = "\t", row.names = F)

cols <- read.delim('~/Projects/Maris-lab/PPTC/data/figures/2018-08-23-all-hist-colors', header = F, stringsAsFactors = F)
cols <- cols[which(cols$V1 %in% final$Histology.Detailed),]
cols <- cols[match(levels(final$Histology.Detailed), cols$V1),] # reorder colors to match histology
group.colors <- setNames(cols$V2, nm = cols$V1)

p <- ggplot(final, aes(x = Histology.Detailed, y = value, color = Histology.Detailed, alpha = 0.5)) + 
  geom_boxplot(outlier.shape = 21, fill = 'white') + 
  geom_jitter(position=position_jitter(width=.1, height=0), shape = 21) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  theme_bw() +
  guides(alpha = FALSE, fill = FALSE) + 
  theme_Publication() + xlab("Histology") + ylab('Number of Fusions') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  guides(color = F) +
  scale_color_manual(values = group.colors) +
  scale_y_continuous(breaks = seq(0, 26, by = 5))
ggsave(filename = 'results/FusionPlot.pdf', plot = p, device = 'pdf', height = 6, width = 10)

medians <- unique(final[,c("Histology.Detailed","median")])
write.table(medians, file = 'results/FusionPlot_medians.txt', quote = F, sep = "\t", row.names = F)
