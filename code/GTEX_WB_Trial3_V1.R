library("data.table")
library(WGCNA)


dir_path="./WGCNA/GTEX/Trial2/"
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
#Read in the female liver data set
femData = fread("Whole-Blood-670-Mean-above1-13k-log2-MinMax_wgcna.csv", header = TRUE, sep = ",");
# Take a quick look at what is in the data set:
dim(femData);  #670 13070 row is sample, column is gene


####transpose the matric##########
datExpr0 = as.data.frame(femData[,-1]) #the first column is the IID
names(datExpr0) = names(femData[,-1])
rownames(datExpr0) = femData$IID  #670 13069


#####check if all samples are good###
gsg = goodSamplesGenes(datExpr0, verbose = 3); # This function iteratively identifies samples and genes with too many missing entries and genes with zero variance.
gsg$allOK
sum(is.na(datExpr0))


#####cluster samples based on the distance
print("Step 1...")
sampleTree = hclust(dist(datExpr0), method = "average");
pdf(file=paste0(dir_path,"1_sampleClustering.pdf"),width=12,height=9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
abline(h = 30, col = "red")
dev.off()   


#####remove all samples below the line
clust = cutreeStatic(sampleTree, cutHeight = 30, minSize = 10)
table(clust)
keepSamples = (clust==1)
datExpr0 = datExpr0[keepSamples, ]
dim(datExpr0)


####2.a.1 choosing soft-thresholding power: analysis of network topology
print("Step 2...")
#powers = c(c(1:10), seq(from = 12, to=60, by=4)) 
powers = c(seq(from = 1, to=80, by=5))
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5) 
pdf(file=paste0(dir_path,"2_scale_independence.pdf"),width=9,height=5)
par(mfrow = c(1,2))
cex1 = 0.8  ##the default is 0.9 but seems bit high for power which equals to 71; cutoff 0.80, power still 41

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
   xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="o",
   main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
   labels=powers,cex=cex1,col="red");
abline(h=0.80,col="red") 

plot(sft$fitIndices[,1], sft$fitIndices[,5],
   xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="o",
   main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

#####2.a.2 Co-expression similarity and adjacency
softPower = 71 
print("Calculating adjacency...")
adjacency = adjacency(datExpr0, power = softPower)  ## gene gene adjacency:  AIJ = (1 + cor(EI ;EJ ))/2

####2.a.3 Topological Overlap Matrix (TOM)
print("Calculating TOM...")
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

###2.a.4 Clustering using TOM
print("Step 3...")
geneTree = hclust(as.dist(dissTOM), method = "average");
pdf(file=paste0(dir_path,"3_gene_clustering.pdf"),width=12,height=9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
   labels = FALSE, hang = 0.04)
dev.off()

#cutting branches, each leaf is a gene. Branches of the dendrogram group together densely interconnected, highly co-expressed genes
print("Step 4...")
minModuleSize = 30 
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                          deepSplit = 2, pamRespectsDendro = FALSE,
                          minClusterSize = minModuleSize);
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
pdf(file=paste0(dir_path,"4_Dynamic_Tree.pdf"),width=8,height=6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                  dendroLabels = FALSE, hang = 0.03,
                  addGuide = TRUE, guideHang = 0.05,
                  main = "Gene dendrogram and module colors")
dev.off()


###2.a.5 merging modules whose expression profile are very similar
print("Step 5...")
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average")
pdf(file=paste0(dir_path,"5_Clustering_module.pdf"),width=7,height=6)
plot(METree, main = "Clustering of module eigengenes",
   xlab = "", sub = "")
MEDissThres = 0.25 
abline(h=MEDissThres, col = "red")
dev.off()


print("Step 6...")
merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
pdf(file=paste0(dir_path,"6_merged_dynamic.pdf"), width = 9, height = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                  c("Dynamic Tree Cut", "Merged dynamic"),
                  dendroLabels = FALSE, hang = 0.03,
                  addGuide = TRUE, guideHang = 0.05)
dev.off()

moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1   ##table(moduleLabels) contains number of genes in each module 
table(moduleLabels)

MEs = mergedMEs
save(MEs, moduleLabels, moduleColors, geneTree, file = "FemaleLiver-02-networkConstruction-stepByStep.RData")
col_names=names(datExpr0)

module_1 = which(moduleLabels==1)  
module_1_col = col_names[module_1]
module_1_gene_exp = datExpr0[,module_1_col]


module_no = which(moduleLabels!=0)  
module_no_col = col_names[module_no]
module_no_gene_exp = datExpr0[,module_no_col]
datExpr1=module_no_gene_exp

