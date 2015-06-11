# ----------------------------------------------------------------------
# Description : R script for microarray analysis to identify differentiated 
#               expressed genes in the WAT of patDP/+ mice compared with WT
# Author      : Xiaoxi Liu
# Date        : 2015-06-11
# ----------------------------------------------------------------------


# ----------------------------------------------------------------------
# Load the limma package, install from bioconductor if needed
# source("http://bioconductor.org/biocLite.R")
# biocLite("limma")
rm(list=ls())
library("limma")
# ----------------------------------------------------------------------


# ----------------------------------------------------------------------
# Download the GSE58191_RAW.tar from GEO repository at 
# http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE58191
# untar the file
# Create a text file named "targets.txt" inside the GSE58191_RAW folder
# add the following lines to "targets.txt"
# SampleNumber  FileName	Condition
# 1	GSM1402845_WT1.txt.gz	WT
# 2	GSM1402846_WT2.txt.gz	WT
# 3	GSM1402847_WT3.txt.gz	WT
# 4	GSM1402848_pat1.txt.gz	PAT
# 5	GSM1402849_pat2.txt.gz	PAT
# 6	GSM1402850_pat3.txt.gz	PAT
# ----------------------------------------------------------------------


# ----------------------------------------------------------------------
setwd("~/GSE58191_RAW")
targets <- readTargets("targets.txt")
# read single color image
raw_signal <- read.maimages(targets, path="./", source="agilent",green.only=TRUE)

# background correction
processed_signal <- backgroundCorrect(raw_signal, method="normexp", offset=16)

# quantile normalization
processed_signal <- normalizeBetweenArrays(processed_signal, method="quantile")

# this is the averaged signal intensity to probe level 
processed_signal_aveaged <- avereps(processed_signal, ID=processed_signal$genes$ProbeName)

# Export Expression for GEO data deposit 
processed_signal_aveaged$out <- cbind(processed_signal_aveaged$genes, processed_signal_aveaged$E)
write.table(processed_signal_aveaged$out,file="Probe-level-expression-data-with-cBind.txt",sep="\t",quote=FALSE)


# Visualization each sample's expression distribution
# pairs(processed_signal_aveaged$E[,1:3])
# pairs(processed_signal_aveaged$E[,4:6])
# boxplot(processed_signal_aveaged$E)   


# Set up design matrix 
f <- factor(targets$Condition, levels = unique(targets$Condition))
design <- model.matrix(~ 0 + f)
colnames(design) <- levels(f)

# fit the linear model 
fit <- lmFit(processed_signal_aveaged, design)

contrast.matrix <- makeContrasts("PAT-WT", levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)

# eBayes modest T statistics 
fit2 <- eBayes(fit2)

# Save the final results to file "Pat_vs_WT_DEG.txt"
output <- topTable(fit2, adjust="BH", coef="PAT-WT", genelist=processed_signal_aveaged$genes, number=Inf)
write.table(output, file="Pat_vs_WT_DEG.txt", sep="\t", quote=TRUE)

top <- subset(output,P.Value <= 0.05 & abs(logFC) >= log2(1.3))
write.table(top, file="Pat_vs_WT_top_p_0.05_FC_1.3.txt", sep="\t", quote=TRUE)
nrow(top) #340

# Only retain the coding genes
coding <- grep("^NM.", top$SystematicName)
top_coding <- top[coding,]
nrow(top_coding)  #262

# Save the top DEGs which are coding genes 
write.table(top_coding, file="top_coding_Pat_vs_WT.txt", sep="\t", quote=TRUE)


# Prepare for Table 1
top_coding_ordered_ascending <- top_coding[order(top_coding$logFC),]
top_coding_ordered_ascending <- subset(top_coding_ordered_ascending,logFC<0)


top_coding_ordered_decending <- top_coding[order(-top_coding$logFC),]
top_coding_ordered_decending <- subset(top_coding_ordered_decending, logFC>0)


length(unique(top_coding_ordered_ascending$GeneName)) #85 genes
length(unique(top_coding_ordered_decending$GeneName)) #145 genes # Total: 230

gene_ascending <- paste(unique(top_coding_ordered_ascending$GeneName),sep=",",collapse=", ")
gene_decending <- paste(unique(top_coding_ordered_decending$GeneName),sep=",",collapse=", ")

gene_ascending
gene_decending
