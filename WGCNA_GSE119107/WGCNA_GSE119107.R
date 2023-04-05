### Loading thge libraries
library(data.table)
library(GEOquery)
###
setwd("/Users/mokira/MPLAB/R_work/WGCNA_GSE119107")
gset <- getGEO("GSE119107", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL19545", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "000011112222"
sml <- strsplit(gsms, split="")[[1]]

ex <- exprs(gset)
exprs
###

#______________________________________________________________________________#
### WGCNA analysis

library("WGCNA")
#data0=read.table("final_list_DEG_2881_FG.tsv",header = T,sep = "\t")
newex <- t(ex)
sampleTree = hclust(dist(newex), method = "average");


## ggplot for hclust



# plot sample tree
pdf(file = "1-n-sampleClustering.pdf", width = 12, height = 9);
par(cex = 1.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",
     cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
dev.off()

#plot sample tree
tiff("1-n-sampleClustering.tiff", height = 1200, width = 1200, res = 300, 
     pointsize = 12)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",
     hang = 0.1)
dev.off()
# choose a set of soft threshold parameters
powers = c(c(1:30), seq(from = 22, to=40, by=2))
sft = pickSoftThreshold(newex, powerVector = powers, verbose = 5) 
# Scale-free topology fit index as a function of the soft-thresholding power
pdf(file = "2-n-sft.pdf", width = 9, height = 5);
tiff("2-n-sft.tiff", height = 1300, width = 1300,res = 300)
par(mfrow = c(1,2));
#cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red") 
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,col="red")
dev.off()
#===============================================================================
#
#  Turn data expression into topological overlap matrix
#
#===============================================================================

# Turn data expression into topological overlap matrix
power=sft$powerEstimate
power = 30
TOM = TOMsimilarityFromExpr(newex, power = power)
dissTOM = 1-TOM 
# Plot gene tree
geneTree = hclust(as.dist(dissTOM), method = "average");
pdf(file = "3-gene_cluster.pdf", width = 12, height = 9);
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
dev.off()




#=================================Alternate option 1: automatic
cor <- WGCNA::cor
net = blockwiseModules(newex, power = power,
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0,
                       mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = F,
                       saveTOMs = FALSE,
                       verbose = 3)
cor<- stats::cor

sizeGrWindow(12,9)
mergedColors = labels2colors(net$colors)
pdf(file = "4-module_tree_blockwise.pdf", width = 8, height = 6);
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colours")
dev.off()

#===============================================================================
#
#  Construct modules
#
#===============================================================================

# Module identification using dynamic tree cut
# the minClusterSize controls the number of genes presentin final analysis
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 2, 
                            pamRespectsDendro = FALSE,minClusterSize = 10);
table(dynamicMods)
length(table(dynamicMods)) 
# Convert numeric labels into colors
dynamicColors = labels2colors(dynamicMods)
# table(dynamicColors)
# # Plot the dendrogram and colors underneath
# pdf(file = "4-module_tree.pdf", width = 8, height = 6);
# plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",dendroLabels = FALSE,
#                     hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors")
# dev.off()

#===============================================================================
#
#  Merge modules

MEDissThres=0.2
#abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(newex, dynamicColors, cutHeight = MEDissThres, verbose = 3) 
mergedColors = merge$colors  
mergedMEs = merge$newMEs


setwd("/Users/mokira/MPLAB/R_work/WGCNA_GSE119107/Cytoscape2")

for (i in 1:length(merge$newMEs)){
  modules = c(substring(names(merge$newMEs)[i], 3));
  genes = colnames(newex)
  inModule = is.finite(match(dynamicColors,modules))
  modGenes = genes[inModule]
  modTOM=TOM[inModule,inModule]
  dimnames(modTOM)=list(modGenes,modGenes)
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("merge_CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                                 nodeFile = paste("merge_CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                                 weighted = TRUE, threshold = -1, nodeNames = modGenes, nodeAttr = dynamicColors[inModule]);
}

