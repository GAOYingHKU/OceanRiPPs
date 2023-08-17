library(WGCNA)
library(reshape2)
library(stringr)
library(DESeq2)
library(magrittr)

options(stringsAsFactors = FALSE)
enableWGCNAThreads()

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


rownames(normalized_counts)
nSamples = nrow(normalized_counts)

datTraits <- read.csv("trait1.txt", header=T, as.is=T,sep='\t')
rownames(datTraits)<-datTraits[,1]
datTraits<-datTraits[,-1]

load('net.wgcna.RData')
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
MEs0 <- moduleEigengenes(normalized_counts, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
moduleTraitCor <- cor(MEs, datTraits, method ="pearson",use="complete.obs")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
moduleTraitCor
moduleTraitPvalue
png(filename = "Module-Trait-Relationship1.png", width = 20, height = 120, res=600, unit="cm")
par(mar = c(6, 8.5, 3, 3))
textMatrix <- paste0(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")",sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
labeledHeatmap(
 Matrix = moduleTraitCor,
 xLabels = names(datTraits),
 yLabels = names(MEs),
 ySymbols = names(MEs),
 colorLabels = FALSE,
 colors = blueWhiteRed(50),
 textMatrix = textMatrix,
 setStdMargins = FALSE,
 cex.text = 0.5,
 zlim = c(-1,1),
 main = paste("Module-trait relationships"))
dev.off()
