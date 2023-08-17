library(WGCNA)
library(reshape2)
library(stringr)
library(DESeq2)
library(magrittr)

options(stringsAsFactors = FALSE)


exprMat <- "TOM"
metadata <- readr::read_tsv("metadata.txt")
dataExpr<-read.table('input_table_clean_2.txt',sep='\t',row.names=1,header=T,quote="", comment="", check.names=F)
dds <- DESeqDataSetFromMatrix(
  countData = dataExpr, # Our prepped data frame with counts
  colData = metadata, # Data frame with annotation for our samples
  design = ~1 # Here we are not specifying a model
)
dds_norm <- vst(dds)
normalized_counts <- assay(dds_norm) %>%
  t()




load('net.wgcna.RData')
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
probes = colnames(normalized_counts)



load(net$TOMFiles[2], verbose=T)
subprobes2=probes[as.numeric(unlist(net$blockGenes[2]))]
submoduleColors2=moduleColors[as.numeric(unlist(net$blockGenes[2]))]
TOM2 <- as.matrix(TOM)
dimnames(TOM2) <- list(subprobes2, subprobes2)
cyt = exportNetworkToCytoscape(TOM2,
             edgeFile = paste(exprMat, ".02.edges.txt", sep=""),
             nodeFile = paste(exprMat, ".02.nodes.txt", sep=""),
             weighted = TRUE, threshold = 0,
             nodeNames = subprobes2, nodeAttr = submoduleColors2)
