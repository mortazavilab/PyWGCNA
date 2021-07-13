library(WGCNA)
library(RColorBrewer)
library(cowplot)
library(Hmisc)
library("ComplexHeatmap")
ht_opt$message = FALSE
library(stringr)
options(stringsAsFactors = FALSE);

###---------------------- 
# Read data
print("Read Data")
expressionList = read.table('expressionList', header = TRUE);

## Prepare and clean data
#Remove rows with less than 1 TPM
expressionList = expressionList[expressionList[,ncol(expressionList)]>1,]

datExpr0 = as.data.frame(t(expressionList[,-c(1)]));
names(datExpr0) = expressionList$gene_id
rownames(datExpr0) = names(expressionList)[-c(1)];

# Check that all genes and samples have sufficiently low numbers of missing values.
gsg = goodSamplesGenes(datExpr0, verbose = 3);
#if not okay 
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

## Clustering
sampleTree = hclust(dist(datExpr0), method = "average");
# The user should change the dimensions if the window is too large or too small.
pdf(file = "../sampleClusteringCleaning.pdf", width = 25, height = 6);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
# Plot a line to show the cut
abline(h = 50000, col = "red");
dev.off();

# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 50000, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust!=0)

datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


#convert trancript ID to gene ID
datExpr = as.data.frame(t(datExpr))
for (i in c(1:dim(datExpr)[1])) {
  row.names(datExpr)[i] = strsplit(row.names(datExpr)[i], "\\.")[[1]][1]
}
datExpr = as.data.frame(t(datExpr))

collectGarbage();

save(datExpr, file = "../data/data_input.RData")

###---------------------- 
## Modules construction
print("WGCNA")

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = "signed")
# Plot the results:
pdf(file = "../summarypower.pdf", width = 10, height = 5);
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.9,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off();

## Set Power
softPower = 13;
adjacency = adjacency(datExpr, power = softPower, type = "signed");

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency, TOMType = "signed");
dissTOM = 1-TOM
save(TOM, file = "../result_power5_signed/data/TOM.RData")

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
pdf(file = "../dendrogram.pdf", width = 18, height = 6);
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
dev.off();


# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 50;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
# Plot the dendrogram and colors underneath
pdf(file = "../DynamicTreeCut.pdf", width = 18, height = 6);
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off();


# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
pdf(file = "../eigenesgenes.pdf", width = 24, height = 5);
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.2
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off();
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;
pdf(file = "../geneDendro-3.pdf", wi = 18, he = 9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off();

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, dynamicColors, file = "../result_power5_signed/data/Data-networkConstruction.RData")

rm(list = ls())

###---------------------- 
## ANALYSING
print("Analysing")
setwd("../result_power5_signed//")
load(file = "data/data_input.RData");
datTraits = read.csv('../data/datTraits', stringsAsFactors = FALSE, header = TRUE, row.names = 1);
load(file = "data/TOM.RData");
load(file = "data/Data-networkConstruction.RData");
annot = read.table('../data/geneList', header = TRUE, stringsAsFactors = FALSE);

datTraits = datTraits[-c(1,2),]

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels 
datME = moduleEigengenes(datExpr, moduleColors)$eigengenes
datME$MEgrey = NULL
MEs = orderMEs(datME)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

## names
xlabels = c("Age(day)", "Sex")
ylabels = c()
for (label in names(MEs)) {
  ylabels = c(ylabels, substr(label, 3, 10000))
}
ylabels = capitalize(ylabels)

## Quantifying module-trait associations
pdf(file = "Module-traitRelationships.pdf", width = 4, height = 6);
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(5, 7, 1, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = xlabels,
               yLabels = ylabels,
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               cex.lab = 1,
               font.lab.x = 2,
               font.lab.y = 2,
               main = "")
dev.off();


## Quantifying module-trait associations
pdf(file = "Module-traitRelationships_t.pdf", width = 10, height = 3);
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(5, 5, 1, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = t(moduleTraitCor),
               xLabels = ylabels,
               yLabels = xlabels,
               xSymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = t(textMatrix),
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               cex.lab = 1,
               font.lab.x = 2,
               font.lab.y = 2,
               main = "")
dev.off();


## Plot the dendrogram
pdf(file = "EigengeneDendrogram.pdf", width = 6, height = 4);
par(cex = 1.0)
plotEigengeneNetworks(datME, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
dev.off();

pdf(file = "EigengeneDendrogramHeatmap.pdf", width = 6, height = 4);
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off();

## Network visuallization
# Select module probes
modules = names(table(moduleColors))
for (module in modules) {
  # Select module probes
  probes = names(datExpr)
  inModule = is.finite(match(moduleColors, module));
  modProbes = probes[inModule];
  # Select the corresponding Topological Overlap
  modTOM = TOM[inModule, inModule];
  
  cytB = exportNetworkToCytoscape(modTOM,nodeFile = paste("data/", paste(module, "Node.txt",sep="" ),sep="" ), weighted=TRUE, threshold = 0,
                                  nodeNames = modProbes,nodeAttr = moduleColors[inModule])
  
}


## convert gene ID to gene Name for each module
modules = names(table(moduleColors))
for (module in modules) {
  node = read.csv(paste0("data/", module, "Node.txt"), sep = "\t")
  colnames(node) = c("Gene ID", "Gene symbol", "Module class")
  for (i in c(1:dim(node)[1])) {
    node[i,1] = strsplit(node[i,1], "\\.")[[1]][1]
    if (length(which(annot$gene_id == node[i,1])) != 0) {
      node[i,2] = annot$gene_name[which(annot$gene_id == node[i,1])]
    }
    else {
      node[i,2] = node[i,1]
    }
  }
  write.csv(node, paste0("data/", module, "Node.txt"), row.names=F, quote = F)
}

## find size of each module
modules = names(table(moduleColors))
len = c()
for (module in modules) {
  node = read.csv(paste("data/", paste(module, "Node.txt",sep="" ), sep="" ))
  len = c(len, dim(node)[1])
}

len = cbind(modules, len)
len = as.data.frame(len)
colnames(len) = c("moduleColor", "size")
write.csv(len, 'data/moduleSize.txt', row.names = FALSE, quote = F)


len = read.csv('data/moduleSize.txt')
pdf(file = "sizeModule.pdf");
barplot(as.integer(len$size),
        #names.arg = len$moduleColor,
        col = len$moduleColor,
        main = "Size of each module",
        ylab = "number of genes",
        xlab = "modules")
dev.off()


## Diagnostics: displaying module heatmap and the eigengene
## Add color bar
day = rep(0, dim(datExpr)[1])
sex = rep("Male", dim(datExpr)[1])
for (i in c(1:dim(datExpr)[1])) {
  if (datTraits$age.days.[i] == 540) {
    day[i] = "18-20 months"
  } else if (datTraits$age.days.[i] == 60) {
    day[i] = "2 months"
  } else {
    day[i] = paste0(datTraits$age.days.[i], " days")
  }
  
  if (datTraits$sex[i] == 0){
    sex[i] = "Female"
  }
}


## Plot module heatmap eigen gene
modules = names(table(moduleColors))
modules = modules[modules != "grey"]
for (which.module in modules) {
  # Recalculate MEs with color labels
  label = t(scale(datExpr[, moduleColors==which.module]))
  colnames(label) = paste0("sample_", seq(ncol(label)))
  ME=datME[, paste("ME",which.module, sep="")]
  
  for (i in c(1:length(row.names(label)))) {
    index = which(row.names(label)[i] == annot$Gene.ID)
    if (length(index) != 0) {
      if (annot$Gene.Symbol[index] %in% c("Cst7", "Gfap", "Apoe", "Trem2", "Thy1", "Gm4924", "Gas5", "Tyrobp", "Ctsz", "Ccl6", "Ctss")) {
        row.names(label)[i] = annot$Gene.Symbol[index]
      }
      else {
        row.names(label)[i] = ""
      }
    } else {
      row.names(label)[i] = ""
    }
  }
  
  pdf(file = paste0("ModuleHeatmapEigengene_", which.module, ".pdf"), width = 15, height = 8);
  column_ha = HeatmapAnnotation(Timepoint = day,
                                Sex = sex,
                                col = list(Timepoint = c("4 days" = "blue", "10 days" = "darkblue",
                                                         "14 days" = "deeppink", "25 days" = "darkviolet",
                                                         "36 days" = "cyan", "2 months" = "deepskyblue",
                                                         "18-20 months" = "palegreen"),
                                           Sex = c("Female" = "green", "Male" = "yellow")),
                                eigengeneExp = anno_barplot(ME, baseline = 0, gp = gpar(fill = which.module)),
                                show_annotation_name = c(Timepoint = F, Sex = F, eigengeneExpression = T),
                                gap = unit(2, "points"), show_legend = TRUE, annotation_height = c(1,1,5), height = unit(6, "cm"))
  heatmap = Heatmap(label, column_title = which.module, cluster_rows = TRUE, cluster_columns = FALSE, show_row_dend = FALSE,
                    show_row_names = TRUE, show_column_names = FALSE, top_annotation = column_ha, 
                    column_title_gp = gpar(fontsize = 20, fontface = "bold"), show_heatmap_legend = FALSE)
  print(heatmap)
  
  dev.off();
}


## Plot avg exp
modules = names(table(moduleColors))
modules = modules[modules != "grey"]

len = read.csv('data/moduleSize.txt')

barplots = vector(mode = "list", length = length(modules))
for (which.module in modules) {
  # Recalculate MEs with color labels
  label = t(scale(datExpr[, moduleColors==which.module]))
  colnames(label) = paste0("sample_", seq(ncol(label)))
  ME=datME[, paste("ME",which.module, sep="")]
  row.names(label) = NULL
  
  barplots[[which(which.module == modules)]] = ggplot(as.data.frame(ME), aes(x=c(1:length(ME)), y=ME)) + 
    geom_bar(color = "black", fill = which.module, stat = "identity") +
    ylab(paste0(which.module, " module (", len$size[len$moduleColor == which.module], " genes)")) +
    ylim(-1, 1) + 
    xlim(0.5, length(ME)+0.5) + 
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_text(colour="black", size = 18, face = "bold"),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
}

pdf(file = "ModuleHeatmapEigengene.pdf", width = 15, height = 2+4*length(modules));
par(oma=c(3,3,3,5))
annotDay = ggplot(as.data.frame(day)) +
  geom_bar(mapping = aes(x = c(1:length(day)), y = 0.1, fill = day), stat = "identity", width = 1) +
  labs(fill = "Time point") +
  scale_fill_manual(values=c("4 days" = "blue", "10 days" = "darkblue",
                             "14 days" = "deeppink", "25 days" = "darkviolet",
                             "36 days" = "cyan", "2 months" = "deepskyblue",
                             "18-20 months" = "palegreen")) +
  xlim(0.5, length(day)+0.5) + 
  theme_void() +
  theme(legend.text=element_text(colour="black", size = 16),
        legend.title=element_text(colour="black", size = 16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank())

annotSex = ggplot(as.data.frame(sex)) +
  geom_bar(mapping = aes(x = c(1:length(sex)), y = 0.1, fill = sex), stat = "identity", width = 1) +
  labs(fill = "Sex") +
  scale_fill_manual(values=c("Female" = "green", "Male" = "yellow")) +
  xlim(0.5, length(sex)+0.5) + 
  theme_void() +
  theme(legend.text=element_text(colour="black", size = 16),
        legend.title=element_text(colour="black", size = 16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank())


legend = plot_grid(get_legend(annotDay), get_legend(annotSex), ncol = 1)
annotDay = annotDay + theme(legend.position = "none")
annotSex = annotSex + theme(legend.position = "none")
annot = plot_grid(annotDay, annotSex, align = "v", ncol = 1, axis = "tb", rel_heights = c(1, 1))

plot = plot_grid(annot, barplots[[1]], barplots[[2]], barplots[[3]], barplots[[4]], barplots[[5]], barplots[[6]], 
                 barplots[[7]], barplots[[8]], barplots[[9]], barplots[[10]], barplots[[11]], 
                 align = "v", ncol = 1, axis = "tb", rel_heights = c(1, rep(4, length(modules))))
plot_grid(plot, legend, nrow = 1, rel_widths = c(10, 2))

dev.off()



pdf(file = "ModuleHeatmapEigengeneV2.pdf", width = 50, height = 2+2*length(modules));
par(oma=c(3,3,3,5))
annotDay = ggplot(as.data.frame(day)) +
  geom_bar(mapping = aes(x = c(1:length(day)), y = 0.1, fill = day), stat = "identity", width = 1) +
  labs(fill = "Time point") +
  scale_fill_manual(values=c("4 days" = "blue", "10 days" = "darkblue",
                             "14 days" = "deeppink", "25 days" = "darkviolet",
                             "36 days" = "cyan", "2 months" = "deepskyblue",
                             "18-20 months" = "palegreen")) +
  xlim(0.5, length(day)+0.5) + 
  theme_void() +
  theme(legend.text=element_text(colour="black", size = 16),
        legend.title=element_text(colour="black", size = 16, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank())

annotSex = ggplot(as.data.frame(sex)) +
  geom_bar(mapping = aes(x = c(1:length(sex)), y = 0.1, fill = sex), stat = "identity", width = 1) +
  labs(fill = "Sex") +
  scale_fill_manual(values=c("Female" = "green", "Male" = "yellow")) +
  xlim(0.5, length(sex)+0.5) + 
  theme_void() +
  theme(legend.text=element_text(colour="black", size = 16),
        legend.title=element_text(colour="black", size = 16, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank())


legend = plot_grid(get_legend(annotDay), get_legend(annotSex), ncol = 1)
annotDay = annotDay + theme(legend.position = "none")
annotSex = annotSex + theme(legend.position = "none")
annot = plot_grid(annotDay, annotSex, align = "v", ncol = 1, axis = "tb", rel_heights = c(1, 1))

plot = plot_grid(annot, annot, annot, barplots[[1]], barplots[[2]], barplots[[3]], barplots[[4]], barplots[[5]], barplots[[6]], 
                 barplots[[7]], barplots[[8]], barplots[[9]], barplots[[10]], barplots[[11]],
                 align = "v", ncol = 3, axis = "tb", rel_heights = c(1, rep(4, length(modules))))
plot_grid(plot, legend, nrow = 1, rel_widths = c(10, 1))

dev.off()


