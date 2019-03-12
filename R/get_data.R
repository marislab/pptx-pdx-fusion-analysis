################################################
# Author: Komal S Rathi
# Function: Get rawdata to run fusion scripts
# Date: 03/12/2019
# Step 0
################################################

# RNA-seq expression matrix in hg38
system("wget --output-document='data/pptc_rnaseq_hg38_matrix.RData' \
	https://ndownloader.figshare.com/files/14452982")

# RNA-seq expression matrix in hg19 (for pedcbio)
system("wget --output-document='data/pedcbio/2019-02-14-PPTC_FPKM_matrix_withModelID-244.rda' \
	https://ndownloader.figshare.com/files/14452985")

# fusion input
# defuse
system("wget --output-document='data/Defuse_PPTC_RNASeq_hg19.txt' \
	https://ndownloader.figshare.com/files/14460299")
# soapfuse
system("wget --output-document='data/SOAPfuse_PPTC_RNASeq_hg19.txt' \
	https://ndownloader.figshare.com/files/14460305")
# starfusion
system("wget --output-document='data/starfusion_PPTC_RNASeq_hg19.txt' \
	https://ndownloader.figshare.com/files/14460308")
# fusioncatcher
system("wget --output-document='data/fusioncatcher_PPTC_RNASeq_hg19.txt' \
	https://ndownloader.figshare.com/files/14460302")
