################################################
# Author: Komal S Rathi
# Function: Create Venn from Tab1 and Tab2 of Excel Sheet
# Date: 02/18/2018
# Step 4
################################################

library(dplyr)
library(grid)
library(gridExtra)
library(VennDiagram)
library(ggplot2)
library(extrafont)
loadfonts()
# font_import()
setwd('~/Projects/Maris-lab/PPTC_fusion_analysis/')

tab1 <- read.delim('results/Filtered_Annotated_Fusions.txt', stringsAsFactors = F)
tab2 <- read.delim('results/Driver_Fusions.txt', stringsAsFactors = F)
setdiff(tab2$Fused_Genes, tab1$Fused_Genes) # n = 4

# venn diagram 
tab1$Fusion_Type <- NULL
tab1 <- unique(tab1)
tab1 <- tab1[order(tab1$Fused_Genes, tab1$Method),]
tab1 <- tab1[,c('Fused_Genes','Model','Histology.Broad.Label','Method')]
tab2 <- tab2[,c('Fused_Genes','Model','Histology.Broad.Label','Method')]
tab2 <- tab2 %>% 
  group_by(Fused_Genes, Model, Histology.Broad.Label) %>%
  mutate(Method = strsplit(Method, ", ")) %>%
  unnest(Method) %>%
  filter(Method != 'Cytogenetics') %>%
  as.data.frame()

# total of tab1 and tab2
total <- unique(rbind(tab1, tab2))

for.venn <- total %>% 
  group_by(Fused_Genes, Model, Histology.Broad.Label) %>% 
  unique() %>%
  summarise(Method = toString(Method), Method.count = n()) %>%
  as.data.frame()
for.venn.count <- plyr::count(for.venn, 'Method')
a <- paste0('deFuse\n(n=',sum(for.venn.count[grep('deFuse', for.venn.count$Method),2]),')')
b <- paste0('SOAPFuse\n(n=',sum(for.venn.count[grep('SOAPFuse', for.venn.count$Method),2]),')')
c <- paste0('FusionCatcher\n(n=',sum(for.venn.count[grep('FusionCatcher', for.venn.count$Method),2]),')')
d <- paste0('STARFusion\n(n=',sum(for.venn.count[grep('STARFusion', for.venn.count$Method),2]),')')
g <- draw.quad.venn(area1 = sum(for.venn.count[grep('deFuse',for.venn.count$Method),'freq']),
                    area2 = sum(for.venn.count[grep('SOAPFuse',for.venn.count$Method),'freq']), 
                    area3 = sum(for.venn.count[grep('FusionCatcher',for.venn.count$Method),'freq']), 
                    area4 = sum(for.venn.count[grep('STARFusion',for.venn.count$Method),'freq']),
                    n12 = sum(for.venn.count[grep('deFuse.*SOAPFuse',for.venn.count$Method),'freq']), 
                    n13 = sum(for.venn.count[grep('deFuse.*FusionCatcher',for.venn.count$Method),'freq']),
                    n14 = sum(for.venn.count[grep('deFuse.*STARFusion',for.venn.count$Method),'freq']),
                    n23 = sum(for.venn.count[grep('FusionCatcher.*SOAPFuse',for.venn.count$Method),'freq']), 
                    n24 = sum(for.venn.count[grep('SOAPFuse.*STARFusion',for.venn.count$Method),'freq']),
                    n34 = sum(for.venn.count[grep('FusionCatcher.*STARFusion',for.venn.count$Method),'freq']),
                    n123 = sum(for.venn.count[grep('deFuse.*FusionCatcher.*SOAPFuse',for.venn.count$Method),'freq']),
                    n124 = sum(for.venn.count[grep('deFuse.*SOAPFuse.*STARFusion',for.venn.count$Method),'freq']),
                    n134 = sum(for.venn.count[grep('deFuse.*FusionCatcher.*STARFusion',for.venn.count$Method),'freq']),
                    n234 = sum(for.venn.count[grep('FusionCatcher.*SOAPFuse.*STARFusion',for.venn.count$Method),'freq']),
                    n1234 = sum(for.venn.count[grep('deFuse.*FusionCatcher.*SOAPFuse.*STARFusion',for.venn.count$Method),'freq']),
                    fill = c("skyblue", "pink1", "mediumorchid", "blue"), 
                    category = c(a,b,c,d),
                    cex = 1, cat.cex = 1, cat.pos = c(-10, 10, -10, -10), alpha = 0.5, 
                    cat.fontfamily = "Arial")


g <- gridExtra::grid.arrange(gTree(children = g), top=textGrob("Gene Fusions", gp=gpar(fontsize=15,font=8)))
ggsave(plot = g, filename = 'results/FusionPlot_InframeVenn.pdf', device = 'pdf', height = 6, width = 8)

