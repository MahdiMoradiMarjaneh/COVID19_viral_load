################################################################
###################   Deconvolution   ##########################

# install CellCODE

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("sva")

install.packages("remotes")
remotes::install_github("mchikina/CellCODE")

## https://github.com/mchikina/CellCODE/blob/master/vignettes/Vignette.Rnw

myFPKM <- as.matrix(read.table("master_FPKM_PTGs.txt", header=TRUE, sep = "\t", row.names = 1, as.is=TRUE))
log2myFPKM <- log2(myFPKM + 0.001)

library(CellCODE)
data("GSE20300")
data("IRIS")
data("DMAP")

irisTag=tagData(IRIS[,c("Bcell-naïve", "CD4Tcell-N0", "CD8Tcell-N0", "Monocyte-Day0", "Neutrophil-Resting", "NKcell-control" )], 2, max=15, ref=log2myFPKM, ref.mean=F)

colnames(irisTag)=c("Bcell", "TcellCD4", "TcellCD8", "Monocyte", "Neutrophil", "NKcell")

Estimate proportions using the tag genes
grp=as.factor(rep(c(1,2), each=48))
SPVs=getAllSPVs(log2myFPKM, grp=grp, irisTag, method="raw", plot=T, mix.par=0.3)

write.table(SPVs, file = "SPVs.txt", sep = "\t")

#####################################################
###################   DE   ##########################

library(tximport)

samples <- read.table("samples.txt", header = TRUE)
files <- file.path("rsem", paste0(samples$sample_id, ".genes.results"))
names(files) <- paste0("sample", 1:21)
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
head(txi.rsem$counts)

txi.rsem$length[txi.rsem$length == 0] <- 1

library(edgeR)

cts <- txi.rsem$counts

normMat <- txi.rsem$length

normMat <- normMat/exp(rowMeans(log(normMat)))
normCts <- cts/normMat

eff.lib <- calcNormFactors(normCts) * colSums(normCts)

normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)

y <- DGEList(cts)
y <- scaleOffset(y, normMat)

keep <- filterByExpr(y)
y <- y[keep, ]

design <- model.matrix(~samples$Z + samples$Bcell + samples$TcellCD4 + samples$TcellCD8 + samples$Monocyte + samples$Neutrophil + samples$NKcell, data=y$samples)
y <- estimateDisp(y, design, robust=TRUE)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef=2)
x <- topTags(qlf, n=200000)
write.table(x, file = "results.txt", sep = "\t")

################################################################

BiocManager::install('grimbough/biomaRt')

library(biomaRt)

packageVersion('biomaRt')

ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

res <- read.table("results.txt", header = TRUE)

ensg <- res$gene_id

# get stable id
ensg.no_version = sapply(strsplit(as.character(ensg),"\\."),"[[",1)
res$gene_id <- ensg.no_version

BM <- getBM(attributes = c('ensembl_gene_id', 'entrezgene_id', 'hgnc_symbol'), filters = 'ensembl_gene_id', values=ensg.no_version, mart=ensembl)
BM <- na.omit(BM)	
BM <- BM[!duplicated(BM$ensembl_gene_id), ]
rownames(BM) <- BM$ensembl_gene_id
res <- res[res$gene_id %in% rownames(BM), ]
res <- unique(res)
res$entrezgene_id <- BM$entrezgene_id[match(res$gene_id, rownames(BM))]
res$hgnc_symbol <- BM$hgnc_symbol[match(res$gene_id, rownames(BM))]

res$entrezgene_id <- paste("hsa", res$entrezgene_id, sep=":")
res$entrezgene_id <- as.factor(res$entrezgene_id)
	
################################################################

fc <- res$logFC[res$PValue <= .05]
names(fc) <- res$entrezgene_id[res$PValue <= .05]
pv <- res$PValue[res$PValue <= .05]
names(pv) <- res$entrezgene_id[res$PValue <= .05]

ref <- as.character(res$entrezgene_id)

require(graph)
.require(ROntoTools)
kpg <- keggPathwayGraphs("hsa", verbose = FALSE)

kpg <- keggPathwayGraphs("hsa", updateCache = TRUE, verbose = TRUE)
kpg <- setEdgeWeights(kpg, edgeTypeAttr = "subtype",
edgeWeightByType = list(activation = 1, inhibition = -1,
expression = 1, repression = -1),
defaultWeight = 0)
kpn <- keggPathwayNames("hsa")
kpg <- setNodeWeights(kpg, weights = alphaMLG(pv), defaultWeight = 1)
head(nodeWeights(kpg[["path:hsa04110"]]))

peRes <- pe(x = fc, graphs = kpg, ref = ref, nboot = 200, verbose = FALSE)

pa_summary <- Summary(peRes, pathNames = kpn, totalAcc = FALSE, totalPert = FALSE,
        pAcc = TRUE, pORA = TRUE, comb.pv =  c("pAcc", "pORA"), comb.pv.func = compute.normalInv, order.by = "pComb")

write.csv(as.data.frame(pa_summary), 
          file="pa_summary_results_1.csv")

plot(peRes, c("pAcc", "pORA"), comb.pv.func = compute.normalInv, threshold = .01)

################################################################
###################   WGCNA   ##################################
### Reference: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/

library(WGCNA)
library(limma)

options(stringsAsFactors = FALSE)

myData <- as.matrix(read.table("master_FPKM_PTGs_SPVs.txt", header=TRUE, sep = "\t", row.names = 1, as.is=TRUE))
myDataT <- t(myData)
myFPKM <- myDataT[, -c(1:7)]
myFPKMlog2 <- log2(myFPKM + 0.001)

myBatch <- as.data.frame(myDataT[, c(1:1)])
mySPV <- as.data.frame(myDataT[, c(2:7)])

myFPKMlog2Adj <- limma::removeBatchEffect(x = t(myFPKMlog2), batch=NULL, batch2=NULL, covariates=mySPV,
                                 design=matrix(1,96,1))
datExpr0 <- t(myFPKMlog2Adj)

gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0) 
     printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
     printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

##### dendrogram and trait heatmap 

sampleTree = hclust(dist(datExpr0), method = "average");
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
abline(h = 280, col = "red");

traitData = read.table("VL_Zscores.txt", header=TRUE, sep = "\t", as.is=TRUE, na.strings = "-")
dim(traitData)
names(traitData)

mySamples = rownames(datExpr0);
traitRows = match(mySamples, traitData$sample_id)
datTraits = traitData[traitRows,];
rownames(datTraits) = traitData[traitRows, 1];
collectGarbage();

datTraits <- datTraits$VL_Zscore

sampleTree = hclust(dist(datExpr0), method = "average")
traitColors = numbers2colors(as.numeric(datTraits), signed = FALSE)

pdf(file="sample_dendrogram_and_trait_heatmap_with.Outlier.pdf", height = 10, width = 20)
plotDendroAndColors(sampleTree, traitColors,
groupLabels = names(datTraits),
main = "Sample dendrogram and trait heatmap")
abline(h = 280, col = "red");
dev.off()

##### dendrogram and trait heatmap (outlier removed)

clust = cutreeStatic(sampleTree, cutHeight = 280, minSize = 10)
table(clust)
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

traitData = read.table("VL_Zscores.txt", header=TRUE, sep = "\t", as.is=TRUE, na.strings = "-")
dim(traitData)
names(traitData)

mySamples = rownames(datExpr);
traitRows = match(mySamples, traitData$sample_id)
datTraits = traitData[traitRows,];
rownames(datTraits) = traitData[traitRows, 1];
collectGarbage();

datTraits <- datTraits$VL_Zscore

sampleTree2 = hclust(dist(datExpr), method = "average")
traitColors = numbers2colors(as.numeric(datTraits), signed = FALSE)

pdf(file="sample_dendrogram_and_trait_heatmap_no.Outlier.pdf", height = 10, width = 20)
plotDendroAndColors(sampleTree2, traitColors,
groupLabels = names(datTraits),
main = "Sample dendrogram and trait heatmap")
dev.off()

##############################

enableWGCNAThreads()

powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

pdf(file="scale_independence.pdf", height = 5, width = 5)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red")
abline(h=0.90,col="red")
dev.off()

pdf(file="mean_connectivity.pdf", height = 5, width = 5)
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

net = blockwiseModules(datExpr, power = 5,
TOMType = "unsigned", minModuleSize = 30,
reassignThreshold = 0, mergeCutHeight = 0.25,
numericLabels = TRUE, pamRespectsDendro = FALSE,
saveTOMs = TRUE,
saveTOMFileBase = "myTOM",
verbose = 3)       

sizeGrWindow(12, 9)

mergedColors = labels2colors(net$colors)

pdf(file="cluster_dendrogram.pdf", height = 10, width = 20)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
"Module colors",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
dev.off()


moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
file = "SVD-networkConstruction-auto.RData")



options(stringsAsFactors = FALSE);

nGenes = ncol(datExpr);
nSamples = nrow(datExpr);

MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
df <- data.frame(moduleTraitCor, moduleTraitPvalue)
colnames(df) <- c("moduleTraitCor", "moduleTraitPvalue")
write.csv(df, file = "moduleTraitCorPvalue.csv")

df.sig <- subset(df, moduleTraitPvalue < 0.01) ## update the threshold in nanostring analysis
moduleTraitCor.sig <- subset(moduleTraitCor, rownames(moduleTraitCor) %in% rownames(df.sig))
moduleTraitPvalue.sig <- subset(moduleTraitPvalue, rownames(moduleTraitPvalue) %in% rownames(df.sig))
MEs.sig <- MEs[, colnames(MEs) %in% rownames(df.sig)]


sizeGrWindow(10,6)

textMatrix = paste(signif(moduleTraitCor.sig, 2), "  (",
signif(moduleTraitPvalue.sig, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor.sig)
par(mar = c(6, 8.5, 3, 3));

pdf(file="Module-trait_relationships.pdf", height = 20, width = 10)
labeledHeatmap(Matrix = moduleTraitCor.sig,
xLabels = c("VL_Zscore"),
yLabels = names(MEs.sig),
ySymbols = names(MEs.sig),
colorLabels = FALSE,
colors = greenWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 2,
zlim = c(-1,1), cex.lab.y = 0.7, 
main = paste("Module-trait relationships"))
dev.off()

############

VL_Zscores = as.data.frame(datTraits);
names(VL_Zscores) = "VL_Zscore"

modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, VL_Zscores, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(VL_Zscores), sep="");
names(GSPvalue) = paste("p.GS.", names(VL_Zscores), sep="")



module = "bisque4"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));

pdf(file="bisque4.pdf", height = 10, width = 10)

verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                           abs(geneTraitSignificance[moduleGenes, 1]),
                           xlab = paste("Module Membership in", module, "module"),
                           ylab = "Gene significance for VL z-score",
                           main = paste("Module membership vs. gene significance\n"),
                           abline = TRUE, 
                           abline.color = "red", pch = 20, cex = 5, col = "black", cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2)
dev.off()


colnames(datExpr)[moduleColors=="bisque4"]

probes = colnames(datExpr)

geneInfo0 = data.frame(substanceBXH = probes,
geneSymbol = colnames(datExpr),
LocusLinkID = colnames(datExpr),
moduleColor = moduleColors,
geneTraitSignificance,
GSPvalue)

modOrder = order(-abs(cor(MEs, VL_Zscores, use = "p")));


for (mod in 1:ncol(geneModuleMembership))
{
oldNames = names(geneInfo0)
geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
MMPvalue[, modOrder[mod]]);
names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

write.csv(geneInfo0, file = "geneInfo0.csv")

hubGenes <- chooseTopHubInEachModule(datExpr, moduleColors, power = 5, type = "unsigned") 
write.csv(hubGenes, file = "hubGenes.csv")

# Visualizing the gene network

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 5);

plotTOM = dissTOM^7;

diag(plotTOM) = NA;

sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")


nSelect = 400

set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];

selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];

sizeGrWindow(9,9)

plotDiss = selectTOM^7;
diag(plotDiss) = NA;

pdf(file="network_heatmap_400genes.pdf", height = 10, width = 10)
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
dev.off()

# Visualizing the network of eigengenes

MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes

VL_Zscore = datTraits; 
names(VL_Zscore) = "VL_Zscore"

MET = orderMEs(cbind(MEs, VL_Zscore))

sizeGrWindow(5,7.5);
par(cex = 0.9)

pdf(file="eigengenes_network.pdf", height = 15, width = 20)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle = 90)
dev.off()

# Visualizing the network of eigengenes (sig modules only)

sig.modules <- c("MEbisque4", "MEthistle2", "MEnavajowhite2", "MEskyblue1", "MElightcyan1", "MEmediumorchid", "MEthistle1", "MEdarkolivegreen", "MEgreenyellow", "MEmagenta", "MEdarkturquoise", "MElightsteelblue1", "MElightcoral", "MEindianred4", "MEivory", "MEsalmon4", "MEdarkgrey", "MEcoral2", "MEantiquewhite4", "MEskyblue3", "MEgrey60", "MEblack", "MEcyan", "MEdarkred")
sig.MEs <- MEs[,colnames(MEs) %in% sig.modules] 
sig.MET = orderMEs(cbind(sig.MEs, VL_Zscore))
pdf(file="sig.eigengenes_network.pdf", height = 10, width = 15)
plotEigengeneNetworks(sig.MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle = 90)
dev.off()

################################################################
###################   BloodGen3Module   ########################
## Reference: https://github.com/Drinchai/BloodGen3Module

library(limma)
library(BloodGen3Module)
library(S4Vectors)

exp_matrix <- as.matrix(read.table("master_FPKM_PTGs.txt", header=TRUE, sep = "\t", row.names = 1, as.is=TRUE))
sample_info <- as.data.frame(read.table("sample_types.txt", header=TRUE, sep = "\t", row.names = 1, as.is=TRUE))

Group_df <- Groupcomparison(exp_matrix,
                            sample_info = sample_info,
                            FC = 0.5,
                            pval = 0.05,
                            FDR = FALSE,
                            Group_column = "sample_type",
                            Test_group = "high",
                            Ref_group = "low")
							
percentages_ttest <- Group_df@assays@data$Percent

write.table(percentages_ttest, file="percentages_ttest.txt", sep="\t", quote=F, col.names=NA)

gridplot(Group_df, 
         cutoff = 10, 
         Ref_group = "low",
         filename= "gridplot_ttest_cutoff10" )

gridplot(Group_df, 
         cutoff = 15, 
         Ref_group = "low",
         filename= "gridplot_ttest_cutoff15")		 
							
							
Group_df <- Groupcomparison(exp_matrix,
                            sample_info = sample_info,
                            FC = 0.5,
                            pval = 0.05,
                            FDR = TRUE,
                            Group_column = "sample_type",
                            Test_group = "high",
                            Ref_group = "low")		


Group_limma <- Groupcomparisonlimma(exp_matrix,
                                    sample_info = sample_info,
                                    FC = 0.5,
                                    pval = 0.05,
                                    FDR = FALSE,
                                    Group_column = "sample_type",
                                    Test_group = "high",
                                    Ref_group = "low")

percentages_limma <- Group_limma@assays@data$Percent
write.table(percentages_limma, file="percentages_limma.txt", sep="\t", quote=F, col.names=NA)

gridplot(Group_limma, 
         cutoff = 10, 
         Ref_group = "low",
         filename= "gridplot_limma_cutoff10" )

gridplot(Group_limma, 
         cutoff = 15, 
         Ref_group = "low",
         filename= "gridplot_limma_cutoff15")	

################################################################
###################   severity plot   ##########################

library(dplyr)
library(ggplot2)

data <- read.table("severity_z-score.txt", header = TRUE, row.names = 1, check.names=FALSE)
pdf(file="severity_z-score_16.samples.pdf", height = 10, width = 10)
ggplot(data) + geom_point(aes(severity, z_score), size = 10) + theme_bw() +
    stat_summary(aes(severity, z_score), fun = mean, geom = "pointrange",fun.max = function(x) mean(x) + sd(x), fun.min = function(x) mean(x) - sd(x), color="red", size=2)
dev.off()

######################################################
############## volcano plot ##########################

library(ggplot2)

de <- read.csv("nasal_nanostring_genes_zscore_corr.csv", header = TRUE, row.names = 1, check.names=FALSE)

ggplot(data=de, aes(x=correlation, y=pvalue)) + geom_point()

p <- ggplot(data=de, aes(x=correlation, y=-log10(pvalue))) + geom_point()

p <- ggplot(data=de, aes(x=correlation, y=-log10(pvalue))) + geom_point() + theme_minimal()

p2 <- p + geom_vline(xintercept=c(-0.5, 0.5), col="red") +
    geom_hline(yintercept=-log10(0.05), col="red")
	
de$correlated <- "No"
de$correlated[de$correlation > 0.5 & de$pvalue < 0.05] <- "Positive"
de$correlated[de$correlation < -0.5 & de$pvalue < 0.05] <- "Negative"

p <- ggplot(data=de, aes(x=correlation, y=-log10(pvalue), col=correlated)) + geom_point() + theme_minimal()

p2 <- p + geom_vline(xintercept=c(-0.5, 0.5), col="red") +
        geom_hline(yintercept=-log10(0.05), col="red")
		

p3 <- p2 + scale_color_manual(values=c("blue", "black", "red"))

mycolors <- c("blue", "red", "black")
names(mycolors) <- c("Negative", "Positive", "No")
p3 <- p2 + scale_colour_manual(values = mycolors)

de$delabel <- NA
de$delabel[de$correlated != "No"] <- de$gene_symbol[de$correlated != "No"]

ggplot(data=de, aes(x=correlation, y=-log10(pvalue), col=correlated, label=delabel)) + 
    geom_point() + 
    theme_minimal() +
    geom_text()
	
library(ggrepel)

pdf(file="volcano.pdf", height = 10, width = 10)

ggplot(data=de, aes(x=correlation, y=-log10(pvalue), col=correlated, label=delabel)) +
        geom_point() + 
        theme_minimal() +
        geom_text_repel() +
        scale_color_manual(values=c("blue", "black", "red")) +
        geom_vline(xintercept=c(-0.5, 0.5), col="red") +
        geom_hline(yintercept=-log10(0.05), col="red") +
		theme_bw()
		
dev.off()

################################################################
##############   compare nanostring with RNAseq data   #########

# convert table to two columns

data <- as.matrix(read.table("blood_NE.txt", header=TRUE, sep = "\t", row.names = 1, as.is=TRUE))
crossdata <- lapply(rownames(data),function(x)sapply(colnames(data),function(y)list(x,y,data[x,y])))
crossdatatmp <- matrix(unlist(crossdata),nrow=3)
crossdatamat <- t(crossdatatmp)
colnames(crossdatamat) <- c("From","To","Value")
crossdatadf <- as.data.frame(crossdatamat,stringsAsFactors=F)
crossdatadf[,3] <- as.numeric(crossdatadf[,3])
crossdatadf
write.table(crossdatadf, file = "blood_NE_converted.txt")


data <- as.matrix(read.table("RNAseq_NE.filtered.txt", header=TRUE, sep = "\t", row.names = 1, as.is=TRUE))
crossdata <- lapply(rownames(data),function(x)sapply(colnames(data),function(y)list(x,y,data[x,y])))
crossdatatmp <- matrix(unlist(crossdata),nrow=3)
crossdatamat <- t(crossdatatmp)
colnames(crossdatamat) <- c("From","To","Value")
crossdatadf <- as.data.frame(crossdatamat,stringsAsFactors=F)
crossdatadf[,3] <- as.numeric(crossdatadf[,3])
crossdatadf
write.table(crossdatadf, file = "RNAseq_NE.filtered_converted.txt")