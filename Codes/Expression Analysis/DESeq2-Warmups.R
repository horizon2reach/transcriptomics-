
# Read in the raw read counts

rawCounts <- read.delim("http://genomedata.org/gen-viz-workshop/intro_to_deseq2/tutorial/E-GEOD-50760-raw-counts.tsv")

#head(rawCounts,3)

#rawCounts <- read_tsv("E-GEOD-50760-raw-counts.tsv")
  
# Read in the sample mappings
sampleData <- read.delim("http://genomedata.org/gen-viz-workshop/intro_to_deseq2/tutorial/E-GEOD-50760-experiment-design.tsv")
#head(sampleData,3)

#sampleData <- read_tsv("E-GEOD-50760-experiment-design (1).tsv")
# Also save a copy for later
sampleData_v2 <- sampleData

# Install the package 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(DESeq2)

# Convert count data to a matrix of appropriate form that DESeq2 can read

geneID <- rawCounts$`Gene ID`

sampleIndex <- grepl("SRR\\d+", colnames(rawCounts))

rawCounts <- as.matrix(rawCounts[,sampleIndex])

#rawcounts <- as.matrix(rawCounts)

rownames(rawCounts) <- geneID

head(rawCounts,3)


# Convert sample variable mappings to an appropriate form that DESeq2 can read
head(sampleData,3)

rownames(sampleData) <- sampleData$Run

keep <- c("Sample.Characteristic.biopsy.site.", "Sample.Characteristic.individual.")

sampleData <- sampleData[,keep]

colnames(sampleData) <- c("tissueType", "individualID")

sampleData$individualID <- factor(sampleData$individualID)

head(sampleData)


# Put the columns of the count data in the same order as rows names of the sample mapping, then make sure it worked
rawCounts <- rawCounts[,unique(rownames(sampleData))]
all(colnames(rawCounts) == rownames(sampleData))


# rename the tissue types
rename_tissues <- function(x){
  x <- switch(as.character(x), "normal"="normal-looking surrounding colonic epithelium", "primary tumor"="primary colorectal cancer",  "colorectal cancer metastatic in the liver"="metastatic colorectal cancer to the liver")
  return(x)
}
sampleData$tissueType <- unlist(lapply(sampleData$tissueType, rename_tissues))

# Order the tissue types so that it is sensible and make sure the control sample is first: normal sample -> primary tumor -> metastatic tumor

sampleData$tissueType <- factor(sampleData$tissueType, levels=c("normal-looking surrounding colonic epithelium", "primary colorectal cancer", "metastatic colorectal cancer to the liver"))

# Create the DEseq2DataSet object

deseq2Data <- DESeqDataSetFromMatrix(countData=rawCounts, colData=sampleData, design= ~ individualID + tissueType)


#######################################################################################

# Data Cleaning 

dim(deseq2Data)

dim(deseq2Data[rowSums(counts(deseq2Data)) > 10, ])

deseq2Data <- deseq2Data[rowSums(counts(deseq2Data)) > 10, ]


#######################################################################################

# Parallel Computations 

library(BiocParallel)
register(MulticoreParam(4))


######################################################################################

# Differential Expression Analysis 

deseq2Data <- DESeq(deseq2Data, parallel = TRUE)

# Result Extraction 

deseq2Results <- results(deseq2Data, contrast=c("tissueType", "primary colorectal cancer", "normal-looking surrounding colonic epithelium"))

summary(deseq2Results)

plotMA(deseq2Results)



deseq2ResDF <- as.data.frame(deseq2Results)

###### HEATMAP 

# Transform count data using the variance stablilizing transform
deseq2VST <- vst(deseq2Data)

# Convert the DESeq transformed object to a data frame
deseq2VST <- assay(deseq2VST)

deseq2VST <- as.data.frame(deseq2VST)

deseq2VST$Gene <- rownames(deseq2VST)

head(deseq2VST)

# Keep only the significantly differentiated genes where the fold-change was at least 3
sigGenes <- rownames(deseq2ResDF[deseq2ResDF$padj <= .05 & abs(deseq2ResDF$log2FoldChange) > 3,])

deseq2VST <- deseq2VST[deseq2VST$Gene %in% sigGenes,]

# Convert the VST counts to long format for ggplot2

library(reshape2)

# First compare wide vs long version
deseq2VST_wide <- deseq2VST

deseq2VST_long <- melt(deseq2VST, id.vars=c("Gene"))

head(deseq2VST_wide)

head(deseq2VST_long)

# Now overwrite our original data frame with the long format

deseq2VST <- melt(deseq2VST, id.vars=c("Gene"))

# Make a heatmap

library(ggplot2)
library(scales) 



install.packages("viridis")
library(viridis)

heatmap <- ggplot(deseq2VST, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
heatmap

