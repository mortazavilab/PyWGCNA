library(WGCNA)
library(RColorBrewer)
library(cowplot)
library(Hmisc)
library("ComplexHeatmap")
ht_opt$message = FALSE
library(stringr)
options(stringsAsFactors = FALSE);
library(analyze.stuff)
library(dynamicTreeCut)

.interpolate <- function (data, index) 
{
  i = round(index)
  n = length(data)
  if (i < 1) 
    return(data[1])
  if (i >= n) 
    return(data[n])
  r = index - i
  data[i] * (1 - r) + data[i + 1] * r
}

CoreSize <- function (BranchSize, minClusterSize) 
{
  BaseCoreSize = minClusterSize/2 + 1
  if (BaseCoreSize < BranchSize) {
    CoreSize = as.integer(BaseCoreSize + sqrt(BranchSize - 
                                                BaseCoreSize))
  }
  else CoreSize = BranchSize
  CoreSize
}

.chunkSize = 100

minClusterSize = 50
dendro = geneTree
distM = dissTOM
deepSplit = 2
pamRespectsDend = FALSE

cutHeight = NULL
method = "hybrid"
maxCoreScatter = NULL
minGap = NULL
maxAbsCoreScatter = NULL
minAbsGap = NULL
minSplitHeight = NULL
minAbsSplitHeight = NULL
externalBranchSplitFnc = NULL
minExternalSplit = NULL
externalSplitOptions = list()
externalSplitFncNeedsDistance = NULL
assumeSimpleExternalSpecification = TRUE
pamStage = TRUE
pamRespectsDendro = TRUE
useMedoids = FALSE
maxDistToLabel = NULL
maxPamDist = cutHeight
respectSmallClusters = TRUE
verbose = 2
indent = 0

dendro = dendro
distM = distM
cutHeight = cutHeight
minClusterSize = minClusterSize
deepSplit = deepSplit
maxCoreScatter = maxCoreScatter
minGap = minGap
maxAbsCoreScatter = maxAbsCoreScatter
minAbsGap = minAbsGap
minSplitHeight = minSplitHeight
minAbsSplitHeight = minAbsSplitHeight
externalBranchSplitFnc = externalBranchSplitFnc
minExternalSplit = minExternalSplit
externalSplitOptions = externalSplitOptions
externalSplitFncNeedsDistance = externalSplitFncNeedsDistance
assumeSimpleExternalSpecification = assumeSimpleExternalSpecification
pamStage = pamStage
amRespectsDendro = pamRespectsDendro
useMedoids = useMedoids
maxPamDist = maxPamDist
respectSmallClusters = respectSmallClusters
verbose = verbose
indent = indent

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
MEs$MEgrey = NULL
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
save(MEs, moduleLabels, moduleColors, geneTree, dynamicColors, file = "../data/Data-networkConstruction.RData")


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


##-----------
.chunkSize = 100

minClusterSize = 50
dendro = geneTree
distM = dissTOM
deepSplit = 2
pamRespectsDend = FALSE

cutHeight = NULL
method = "hybrid"
maxCoreScatter = NULL
minGap = NULL
maxAbsCoreScatter = NULL
minAbsGap = NULL
minSplitHeight = NULL
minAbsSplitHeight = NULL
externalBranchSplitFnc = NULL
minExternalSplit = NULL
externalSplitOptions = list()
externalSplitFncNeedsDistance = NULL
assumeSimpleExternalSpecification = TRUE
pamStage = TRUE
pamRespectsDendro = TRUE
useMedoids = FALSE
maxDistToLabel = NULL
maxPamDist = cutHeight
respectSmallClusters = TRUE
verbose = 2
indent = 0

dendro = dendro
distM = distM
cutHeight = cutHeight
minClusterSize = minClusterSize
deepSplit = deepSplit
maxCoreScatter = maxCoreScatter
minGap = minGap
maxAbsCoreScatter = maxAbsCoreScatter
minAbsGap = minAbsGap
minSplitHeight = minSplitHeight
minAbsSplitHeight = minAbsSplitHeight
externalBranchSplitFnc = externalBranchSplitFnc
minExternalSplit = minExternalSplit
externalSplitOptions = externalSplitOptions
externalSplitFncNeedsDistance = externalSplitFncNeedsDistance
assumeSimpleExternalSpecification = assumeSimpleExternalSpecification
pamStage = pamStage
amRespectsDendro = pamRespectsDendro
useMedoids = useMedoids
maxPamDist = maxPamDist
respectSmallClusters = respectSmallClusters
verbose = verbose
indent = indent

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 50;

for (merge in 1:16) if (dendro$height[merge] <= cutHeight) {
  print(merge)
  if (dendro$merge[merge, 1] < 0 & dendro$merge[merge, 
                                                2] < 0) {
    print("A")
    nBranches = nBranches + 1
    branch.isBasic[nBranches] = TRUE
    branch.isTopBasic[nBranches] = TRUE
    branch.singletons[[nBranches]] = c(-dendro$merge[merge, 
    ], extender)
    branch.basicClusters[[nBranches]] = extender
    branch.mergingHeights[[nBranches]] = c(rep(dendro$height[merge], 
                                               2), extender)
    branch.singletonHeights[[nBranches]] = c(rep(dendro$height[merge], 
                                                 2), extender)
    IndMergeToBranch[merge] = nBranches
    RootBranch = nBranches

  }
  else if (sign(dendro$merge[merge, 1]) * sign(dendro$merge[merge, 
                                                            2]) < 0) {
    print("C")
    clust = IndMergeToBranch[max(dendro$merge[merge, 
    ])]
    if (clust == 0) 
      stop("Internal error: a previous merge has no associated cluster. Sorry!")
    gene = -min(dendro$merge[merge, ])
    ns = branch.nSingletons[clust] + 1
    nm = branch.nMerge[clust] + 1
    if (branch.isBasic[clust]) {
      if (ns > length(branch.singletons[[clust]])) {
        branch.singletons[[clust]] = c(branch.singletons[[clust]], 
                                       extender)
        branch.singletonHeights[[clust]] = c(branch.singletonHeights[[clust]], 
                                             extender)
      }
      branch.singletons[[clust]][ns] = gene
      branch.singletonHeights[[clust]][ns] = dendro$height[merge]
    }
    else {
      onBranch[gene] = clust
    }
    if (nm >= length(branch.mergingHeights[[clust]])) 
      branch.mergingHeights[[clust]] = c(branch.mergingHeights[[clust]], 
                                         extender)
    branch.mergingHeights[[clust]][nm] = dendro$height[merge]
    branch.size[clust] = branch.size[clust] + 1
    branch.nMerge[clust] = nm
    branch.nSingletons[clust] = ns
    IndMergeToBranch[merge] = clust
    RootBranch = clust
  }
  else {
    print("B")
    clusts = IndMergeToBranch[dendro$merge[merge, ]]
    sizes = branch.size[clusts]
    rnk = rank(sizes, ties.method = "first")
    small = clusts[rnk[1]]
    large = clusts[rnk[2]]
    sizes = sizes[rnk]
    branch1 = branch.singletons[[large]][1:sizes[2]]
    branch2 = branch.singletons[[small]][1:sizes[1]]
    spyMatch = FALSE
    if (!is.null(spyIndex)) {
      n1 = length(intersect(branch1, spyIndex))
      if ((n1/length(branch1) > 0.99 && n1/length(spyIndex) > 
           0.99)) {
        printFlush(paste("Found spy match for branch 1 on merge", 
                         merge))
        spyMatch = TRUE
      }
      n2 = length(intersect(branch2, spyIndex))
      if ((n2/length(branch1) > 0.99 && n2/length(spyIndex) > 
           0.99)) {
        printFlush(paste("Found spy match for branch 2 on merge", 
                         merge))
        spyMatch = TRUE
      }
    }
    if (branch.isBasic[small]) {
      coresize = CoreSize(branch.nSingletons[small], 
                           minClusterSize)
      Core = branch.singletons[[small]][c(1:coresize)]
      SmAveDist = mean(colSums(distM[Core, Core, drop = FALSE])/(coresize - 
                                                                   1))
    }
    else {
      SmAveDist = 0
    }
    if (branch.isBasic[large]) {
      coresize = CoreSize(branch.nSingletons[large], 
                           minClusterSize)
      Core = branch.singletons[[large]][c(1:coresize)]
      LgAveDist = mean(colSums(distM[Core, Core])/(coresize - 
                                                     1))
    }
    else {
      LgAveDist = 0
    }
    mergeDiagnostics[merge, ] = c(small, branch.size[small], 
                                  SmAveDist, dendro$height[merge] - SmAveDist, 
                                  large, branch.size[large], LgAveDist, dendro$height[merge] - 
                                    LgAveDist, NA)
    SmallerScores = c(branch.isBasic[small], branch.size[small] < 
                        minClusterSize, SmAveDist > maxAbsCoreScatter, 
                      dendro$height[merge] - SmAveDist < minAbsGap, 
                      dendro$height[merge] < minAbsSplitHeight)
    if (SmallerScores[1] * sum(SmallerScores[-1]) > 
        0) {
      DoMerge = TRUE
      SmallerFailSize = !(SmallerScores[3] | SmallerScores[4])
    }
    else {
      LargerScores = c(branch.isBasic[large], branch.size[large] < 
                         minClusterSize, LgAveDist > maxAbsCoreScatter, 
                       dendro$height[merge] - LgAveDist < minAbsGap, 
                       dendro$height[merge] < minAbsSplitHeight)
      if (LargerScores[1] * sum(LargerScores[-1]) > 
          0) {
        DoMerge = TRUE
        SmallerFailSize = !(LargerScores[3] | LargerScores[4])
        x = small
        small = large
        large = x
        sizes = rev(sizes)
      }
      else {
        DoMerge = FALSE
      }
    }
    if (DoMerge) {
      mergeDiagnostics$merged[merge] = 1
    }
    if (!DoMerge && (nExternalSplits > 0) && branch.isBasic[small] && 
        branch.isBasic[large]) {
      if (verbose > 4) 
        printFlush(paste0("Entering external split code on merge ", 
                          merge))
      branch1 = branch.singletons[[large]][1:sizes[2]]
      branch2 = branch.singletons[[small]][1:sizes[1]]
      if (verbose > 4 | spyMatch) 
        printFlush(paste0("  ..branch lengths: ", 
                          sizes[1], ", ", sizes[2]))
      es = 0
      while (es < nExternalSplits && !DoMerge) {
        es = es + 1
        args = externalSplitOptions[[es]]
        args = c(args, list(branch1 = branch1, branch2 = branch2))
        extSplit = do.call(externalBranchSplitFnc[[es]], 
                           args)
        if (spyMatch) 
          printFlush(" .. external criterion ", es, 
                     ": ", extSplit)
        DoMerge = extSplit < minExternalSplit[es]
        externalMergeDiags[merge, es] = extSplit
        mergeDiagnostics$merged[merge] = if (DoMerge) 
          2
        else 0
      }
    }
    if (DoMerge) {
      print("HHHHH")
      branch.failSize[[small]] = SmallerFailSize
      branch.mergedInto[small] = large
      branch.attachHeight[small] = dendro$height[merge]
      branch.isTopBasic[small] = FALSE
      nss = branch.nSingletons[small]
      nsl = branch.nSingletons[large]
      ns = nss + nsl
      if (branch.isBasic[large]) {
        nExt = ceiling((ns - length(branch.singletons[[large]]))/.chunkSize)
        if (nExt > 0) {
          if (verbose > 5) 
            printFlush(paste("Extending singletons for branch", 
                             large, "by", nExt, " extenders."))
          branch.singletons[[large]] = c(branch.singletons[[large]], 
                                         rep(extender, nExt))
          branch.singletonHeights[[large]] = c(branch.singletonHeights[[large]], 
                                               rep(extender, nExt))
        }
        branch.singletons[[large]][(nsl + 1):ns] = branch.singletons[[small]][1:nss]
        branch.singletonHeights[[large]][(nsl + 1):ns] = branch.singletonHeights[[small]][1:nss]
        branch.nSingletons[large] = ns
      }
      else {
        if (!branch.isBasic[small]) 
          stop("Internal error: merging two composite clusters. Sorry!")
        onBranch[branch.singletons[[small]]] = large
      }
      nm = branch.nMerge[large] + 1
      if (nm > length(branch.mergingHeights[[large]])) 
        branch.mergingHeights[[large]] = c(branch.mergingHeights[[large]], 
                                           extender)
      branch.mergingHeights[[large]][nm] = dendro$height[merge]
      branch.nMerge[large] = nm
      branch.size[large] = branch.size[small] + branch.size[large]
      IndMergeToBranch[merge] = large
      RootBranch = large
    }
    else {
      if (branch.isBasic[large] & !branch.isBasic[small]) {
        x = large
        large = small
        small = x
        sizes = rev(sizes)
      }
      if (branch.isBasic[large] | (pamStage & pamRespectsDendro)) {
        nBranches = nBranches + 1
        branch.attachHeight[c(large, small)] = dendro$height[merge]
        branch.mergedInto[c(large, small)] = nBranches
        if (branch.isBasic[small]) {
          addBasicClusters = small
        }
        else addBasicClusters = branch.basicClusters[[small]]
        if (branch.isBasic[large]) {
          addBasicClusters = c(addBasicClusters, large)
        }
        else addBasicClusters = c(addBasicClusters, 
                                  branch.basicClusters[[large]])
        branch.isBasic[nBranches] = FALSE
        branch.isTopBasic[nBranches] = FALSE
        branch.basicClusters[[nBranches]] = addBasicClusters
        branch.mergingHeights[[nBranches]] = c(rep(dendro$height[merge], 
                                                   2), extender)
        branch.nMerge[nBranches] = 2
        branch.size[nBranches] = sum(sizes)
        branch.nBasicClusters[nBranches] = length(addBasicClusters)
        IndMergeToBranch[merge] = nBranches
        RootBranch = nBranches
      }
      else {
        addBasicClusters = if (branch.isBasic[small]) 
          small
        else branch.basicClusters[[small]]
        nbl = branch.nBasicClusters[large]
        nb = branch.nBasicClusters[large] + length(addBasicClusters)
        if (nb > length(branch.basicClusters[[large]])) {
          nExt = ceiling((nb - length(branch.basicClusters[[large]]))/.chunkSize)
          branch.basicClusters[[large]] = c(branch.basicClusters[[large]], 
                                            rep(extender, nExt))
        }
        branch.basicClusters[[large]][(nbl + 1):nb] = addBasicClusters
        branch.nBasicClusters[large] = nb
        branch.size[large] = branch.size[large] + 
          branch.size[small]
        nm = branch.nMerge[large] + 1
        if (nm > length(branch.mergingHeights[[large]])) 
          branch.mergingHeights[[large]] = c(branch.mergingHeights[[large]], 
                                             extender)
        branch.mergingHeights[[large]][nm] = dendro$height[merge]
        branch.nMerge[large] = nm
        branch.attachHeight[small] = dendro$height[merge]
        branch.mergedInto[small] = large
        IndMergeToBranch[merge] = large
        RootBranch = large
      }
    }
  }
  if (verbose > 2) 
    pind = .updateProgInd(merge/nMerge, pind)
}



print(clust)
#print(nBranches)
print(MxBranches)
print(branch.isBasic[clusts])
print(branch.isTopBasic[nBranches])
print(branch.failSize)
print(branch.rootHeight)
print(branch.size[clusts])
print(branch.nMerge[clusts])
print(branch.nSingletons[clusts])
print(branch.nBasicClusters)
print(branch.mergedInto)
print(branch.attachHeight)
print(branch.singletons[[2]])
#print(branch.basicClusters[[merge]])
print(branch.mergingHeights[[clust]])
print(branch.singletonHeights[[clust]])
