### WGCNA on Astrangia poculata ###

library(DESeq2)
library(WGCNA)
library(tidyverse)


# Load dds object from previous analyses
load("../../DESeq2/Astrangia/temperature_DDS.RData")

# Variance stabilizing transformation, it's important to make sure blind = TRUE 
vst = vst(dds, blind = TRUE)

# Let WGCNA use multi threads (do only once)
allowWGCNAThreads()  

# transpose data frame
input = t(as.data.frame(assay(vst)))

# filter genes that are useful in WGCNA
gsg = goodSamplesGenes(input)
gsg$allOK #FALSE

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(input)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(input)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  input = input[gsg$goodSamples, gsg$goodGenes]
}
gsg = goodSamplesGenes(input)
gsg$allOK #TRUE

# read in trait spreadsheet
traits = read.delim("astrangia_expDesign.csv", sep = ",", row.names = 1)

# cluster samples
sampleTree = hclust(dist(input), method = "average")

# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(traits, signed = TRUE);

# Plot the sample dendrogram and the colors underneath.
pdf("figures/Astrangia_dendro_and_traits.pdf")

plotDendroAndColors(sampleTree,
                    traitColors,
                    groupLabels = names(traits),
                    main = "Sample dendrogram and trait heatmap")

dev.off()

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(input, powerVector = powers, verbose = 5, networkType = "signed")


#    Power SFT.R.sq  slope truncated.R.sq mean.k. median.k. max.k.
# 1      1    0.751 -4.210          0.921  4920.0  4.69e+03   8500
# 2      2    0.951 -3.390          0.961  1180.0  1.02e+03   3500
# 3      3    0.960 -2.850          0.965   388.0  3.06e+02   1890
# 4      4    0.950 -2.560          0.967   161.0  1.13e+02   1180
# 5      5    0.931 -2.360          0.980    81.1  4.61e+01    812
# 6      6    0.883 -2.220          0.981    47.9  2.04e+01    590
# 7      7    0.820 -2.070          0.965    32.3  9.67e+00    446
# 8      8    0.752 -1.890          0.925    24.2  4.82e+00    346
# 9      9    0.680 -1.630          0.803    19.6  2.52e+00    275
# 10    10    0.940 -1.110          0.925    16.8  1.37e+00    223
# 11    12    0.917 -0.969          0.908    13.8  4.44e-01    198
# 12    14    0.877 -0.893          0.885    12.4  1.60e-01    191
# 13    16    0.858 -0.843          0.879    11.7  6.30e-02    188
# 14    18    0.797 -0.808          0.826    11.2  2.60e-02    188
# 15    20    0.788 -0.786          0.821    11.0  1.13e-02    188

# Plot the soft threshold scale independence
pdf("figures/scale_independence.pdf")
cex1 = 0.9;
plot(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",
     type="n",
     main = paste("Scale independence"));

text(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,
     cex=cex1,
     col="red");
abline(h=0.90,col="red") # this line corresponds to using an R^2 cut-off of h

dev.off()


# Mean connectivity as a function of the soft-thresholding power

pdf("figures/Mean_connectivity.pdf")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


# automatic network construction and module detection, choosing a soft threshold power of 14
net = blockwiseModules(input,
                       power = 14,
                       networkType = "signed",
                       minModuleSize = 600,
                       numericLabels = FALSE,
                       mergeCutHeight = .25,
                       saveTOMFileBase = "Astrangia_TOM",
                       verbose = 3)

# number of genes associated with each module
table(net$colors)


# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]],
                    mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05)


moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "Astrangia_networkConstruction.RData")


MEList = moduleEigengenes(input, colors = net$unmergedColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
plot(METree,
     main = "Clustering of module eigengenes",
     xlab = "",
     sub = "")



# Define numbers of genes and samples
nGenes = ncol(input);
nSamples = nrow(input);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(input, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
pdf("figures/Module_trait_relationship.pdf")
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(traits),
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



# Define variable weight containing the weight column of datTrait
weight = as.data.frame(traits$Mapped_Sym);
names(weight) = "Brown Phenotype"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(input, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(input, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");

# intramodular analysis, pulling out eigengenes of interest

module = "brown"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for brown phenotype",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

names(t(input))
names(t(input))[moduleColors=="brown"]


test_list = rownames(t(input))
annot = read.delim(file = "astrangia_iso2go.tab", sep = "\t") %>%
  filter(Protein_ID %in% test_list)

dim(annot)
names(annot)
probes = rownames(t(input))
probes2annot = match(probes, annot$Protein_ID)

# The following is the number or probes without annotation:
sum(is.na(probes2annot))

# Create the starting data frame
geneInfo0 = data.frame(Protein_ID = probes,
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue) %>%
  left_join(annot, by = "Protein_ID")
# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, weight, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Brown.Phenotype));
geneInfo = geneInfo0[geneOrder, ]

write.csv(geneInfo, file = "test.csv")



# Read in the probe annotation
annot = read.csv(file = "GeneAnnotation.csv");
# Match probes in the data set to the probe IDs in the annotation file
probes = names(datExpr)
probes2annot = match(probes, annot$substanceBXH)
# Get the corresponding Locuis Link IDs
allLLIDs = annot$LocusLinkID[probes2annot];
# $ Choose interesting modules
intModules = c("brown", "red", "salmon")
for (module in intModules)
{
  # Select module probes
  modGenes = (moduleColors==module)
  # Get their entrez ID codes
  modLLIDs = allLLIDs[modGenes];
  # Write them into a file
  fileName = paste("LocusLinkIDs-", module, ".txt", sep="");
  write.table(as.data.frame(modLLIDs), file = fileName,
              row.names = FALSE, col.names = FALSE)
}
# As background in the enrichment analysis, we will use all probes in the analysis.
fileName = paste("LocusLinkIDs-all.txt", sep="");
write.table(as.data.frame(allLLIDs), file = fileNam
