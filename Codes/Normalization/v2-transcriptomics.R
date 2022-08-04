####### Transcriptomics & Epigenomics 

##  Install the library required 

# Install Bioconductor 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

install.packages("devtools")
library(devtools)

install.packages("readr")
library(readr)

install.packages("data.table")
library(data.table)

install.packages("tidyr")
library(tidyr)

install.packages("tidyverse")
library(tidyverse)

install.packages("stringr")
library(stringr)

# Installing the data set 

data_Raw <- read_tsv("SRP029880.raw_counts.tsv")

data_Colon <- read_tsv("SRP029880.colData.tsv")
sapply(data_colon,class)

counts <- as.matrix(data_Raw, header = T,sep = '\t')

col_Data <- as.matrix(data_Colon,header=T,sep = '\t')

gene_ID <- counts[,1]

counts_A <- subset(counts, select= c(-width,-1))
sapply(counts_A,class)

### Procedure to remove the reads with gene length of zero 

gene_length <- as.data.frame(subset(counts,select = c(width)))

dim(gene_length)
class(gene_length)
sapply(gene_length,class)

as.data.frame(gene_length) %>% 
 select(width)

as.data.frame(gene_length)%>%
 separate(width,c("Length","id"))

#gene_length %>% seperate(width)

gene_length_N <- separate(as.data.frame(gene_length),width,into = c("Length","id"))
sapply(gene_length_N,class)

gene_length_only <- subset(gene_length_N,select = c(Length))

sapply(gene_length_only,class)

gene_Length_numeric <- as.numeric(unlist(gene_length_only))
sapply(gene_Length_numeric,class)

sapply(gene_Length_numeric,class)

counts_A_1 <- cbind(counts_A,gene_Length_numeric)
sapply(counts_A_1,class)

counts_A_0 <-counts_A_1[apply(counts_A_1, 1, function(row) all(row !=0 )), ]

dim(counts_A_0)
sapply(counts_A_0,class)



counts_A_2 <- subset(counts_A_0,select = c(-gene_Length_numeric))
sapply(counts_A_2,class)

counts_N <- matrix(as.numeric(counts_A_2),
                   ncol=ncol(counts_A_2))

colnames(counts_N) <- c("CASE_2","CASE_3","CASE_4","CASE_5","CTRL_1","CTRL_2","CTRL_3",
                        "CTRL_4","CTRL_5")

class(counts_N)
dim(counts_N)

# Sequence Depth Bias Normalization 

# calculate the " Counts Per Million Reads (CPM)"

CPM_Fun <-function(x) {
  x* 10^6/sum(x)
}

cpm <- apply(counts_N,2,CPM_Fun)

colSums(cpm)
################################################################################
# Calculate RPKM 


gene_long <- subset(counts_A_0,select = c("gene_Length_numeric"))
gene_size <- as.numeric(unlist(gene_long))
sapply(gene_size,class)

RPKM_FUN <- function(x) {
  x*10^9/(sum(x)*gene_size)
}

RPKM <- apply(counts_N,2,RPKM_FUN)

head(RPKM)

colSums(RPKM)
##########################################################################

# TPM 

RPK_FUN <- function(x) {
  x/gene_size/1000
}



rpk <- apply(counts_N,2,RPK_FUN)
class(rpk)

TPM_FUN <- function(x) {
  x/(sum(x)*10^6)
}

TPM <- apply(rpk,2,TPM_FUN)

head(TPM)
colSums(TPM)
###########################################################################

#### Data Exploratory Analysis 

#data <- cbind(gene_ID,TPM)

install.packages("pheatmap")
library(pheatmap)

var_gene  <- apply(TPM,1,var)
head(var_gene)

top_genes <-order(var_gene,decreasing = TRUE)[1:100]
head(top_genes)

selected_genes <- TPM[top_genes,]
pheatmap(selected_genes,
         scale ="row",
         show_rownames = TRUE)


#######################################################################################

# Dimensionality Reduction 

install.packages("stats")
library(stats)

install.packages("ggplot2")
library(ggplot2)

T_genes < t(selected_genes)

pca <- prcomp(T_genes)

summary(pca)


core_result <-cor(TPM)
install.packages("corplot")

colData <- read.table(data_Colon, header = T, sep = '\t', 
                      stringsAsFactors = TRUE)
pheatmap(core_result,
        annotation_col=colData,
        cutree_cols = 2)
############################################################################

# Differential Expression Analysis 

BiocManager::install("DESeq2")

if(!require('DESeq2')) {
  install.packages('DESeq2')
  library('DESeq2')
}

data_Colon <- read_tsv("SRP029880.colData.tsv")

counts <- as.matrix(data_Raw, header = T,sep = '\t')


counts_Data <- subset(counts, select= c(-width))

col_Data <- as.matrix(data_Colon[c(-1),],header=T,sep = '\t')



design_Formula <- "~ group"

dds <- DESeqDataSetFromMatrix(counts_N, 
                              colData = col_Data, 
                              design = as.formula(design_Formula))



print(dds)

dds <- dds[ rowSums(DESeq2::counts(dds)) > 1, ]

dds <- DESeq(dds)

#compute the contrast for the 'group' variable where 'CTRL' 
#samples are used as the control group. 

DEresults = results(dds, contrast = c("group", 'CASE', 'CTRL'))

#sort results by increasing p-value
DEresults <- DEresults[order(DEresults$pvalue),]

# extract differential expression results
DEresults <- results(dds, contrast = c('group', 'CASE', 'CTRL'))
