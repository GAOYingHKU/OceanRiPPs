library(WGCNA)
library(reshape2)
library(stringr)
library(DESeq2)
library(magrittr)

options(stringsAsFactors = FALSE)


type = "signed"
corType='bicor'
corFnc = ifelse(corType=="pearson", cor, bicor)
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
power=9

net = blockwiseModules(normalized_counts, power = power, maxBlockSize = 50000,
                       TOMType = type, minModuleSize = 30,networkType = "signed",
                       numericLabels =TRUE,
                       reassignThreshold = 1e-6, mergeCutHeight = 0.15,
                       saveTOMs=TRUE, corType = corType,
                       maxPOutliers=1, loadTOMs=FALSE,deepSplit = 3,replaceMissingAdjacencies=TRUE,
                       saveTOMFileBase = paste0(exprMat, ".tom"),
                       verbose = 3)
save(net,file='net.wgcna.RData')
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
png(file='dendrogram.png')
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

MEs = net$MEs
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "02-networkConstruction-auto.RData")

MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)
png(file='Eigengene_adjacency_heatmap.png')
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap",
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T,
                      xLabelsAngle = 90)
dev.off()
load(net$TOMFiles[1], verbose=T)
TOM <- as.matrix(TOM)



dissTOM = 1-TOM
plotTOM = dissTOM^7
diag(plotTOM) = NA
png(file='Network_heatmap_plot.png')

probes = colnames(normalized_counts)
subprobes=probes[as.numeric(unlist(net$blockGenes[1]))]
submoduleColors=moduleColors[as.numeric(unlist(net$blockGenes[1]))]
load(net$TOMFiles[1], verbose=T)
TOM <- as.matrix(TOM)
dimnames(TOM) <- list(subprobes, subprobes)
cyt = exportNetworkToCytoscape(TOM,
             edgeFile = paste(exprMat, ".01.edges.txt", sep=""),
             nodeFile = paste(exprMat, ".01.nodes.txt", sep=""),
             weighted = TRUE, threshold = 0,
             nodeNames = subprobes, nodeAttr = submoduleColors)

