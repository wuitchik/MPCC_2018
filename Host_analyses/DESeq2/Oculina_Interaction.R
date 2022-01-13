## Differential expression in Oculina arbuscula
## Interaction between temperature and symbiont state
## Likelihood ratio test

library(DESeq2)
library(tidyverse)

counts = read.csv("Oculina/outlier_removed_host_counts.csv", row.names = 1)
expDesign = read.csv("Oculina/outlier_removed_expDesign.csv") 

# Interaction between treatment and symbiotic state

dds = DESeqDataSetFromMatrix(countData = counts,
                             colData = expDesign,
                             design = ~ Treatment + Sym.Status + Treatment:Sym.Status)
dds = DESeq(dds)

dds = DESeq(dds,
            test="LRT",
            reduced = ~ Treatment + Sym.Status) 

tempXsym = results(dds) 

summary(tempXsym)

# out of 20488 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 0, 0%
# LFC < 0 (down)     : 0, 0%
# outliers [1]       : 17, 0.083%
# low counts [2]     : 0, 0%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

tempXsym = as.data.frame(tempXsym) %>%
  arrange(padj)

write.csv(tempXsym, "Oculina/tempXsym_interaction.csv")

save.image("Oculina/interaction_DDS.RData")
