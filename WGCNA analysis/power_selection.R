library(WGCNA)
library(reshape2)
library(stringr)
library(DESeq2)
library(magrittr)

options(stringsAsFactors = FALSE)

type = "signed"
corType='bicor'
corFnc = ifelse(corType=="pearson", cor, bicor)
robustY = ifelse(corType=="pearson",T,F)
exprMat <- "TOM"
metadata <- readr::read_tsv("metadata.txt")
dataExpr<-read.table('input_table_clean_2.txt',sep='\t',row.names=1,header=T,quote="", comment="", check.names=F)
dds <- DESeqDataSetFromMatrix(
  countData = dataExpr, 
  colData = metadata, 
  design = ~1 
)
dds_norm <- vst(dds)
normalized_counts <- assay(dds_norm) %>%
  t()




gsg = goodSamplesGenes(dataExpr, verbose = 3)

if (!gsg$allOK){ 
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:",
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:",
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}

#head(dataExpr)

nGenes = ncol(normalized_counts)
nSamples = nrow(normalized_counts)
sampleTree = hclust(dist(normalized_counts), method = "average")
png(file='outliers.png')
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
dev.off()
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(normalized_counts,dataIsExpr = TRUE, powerVector=powers, networkType=type, verbose=5,removeFirst=FALSE,corFnc = corFnc,blockSize=40)

par(mfrow = c(1,2))
cex1 = 0.9
png(file='softthreshold1.png')
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.85,col="red")
dev.off()
png(file='softthreshold2.png')
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,
     cex=cex1, col="red")

dev.off()
power = sft$powerEstimate
power
