## Differential expression in Oculina arbuscula
## Cold results
## Wald test

library(DESeq2)
library(tidyverse)

counts = read.csv("Oculina/outlier_removed_host_counts.csv", row.names = 1)
expDesign = read.csv("Oculina/outlier_removed_expDesign.csv") 

dds = DESeqDataSetFromMatrix(countData = counts,
                             colData = expDesign,
                             design = ~ Treatment + Sym.Status)
dds = DESeq(dds)

## Cold 
cold_results = results(dds, alpha = 0.05, contrast = c("Treatment", "Cold", "Control"))
cold_summary = summary(cold_results)

# out of 20488 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 1204, 5.9%
# LFC < 0 (down)     : 1663, 8.1%
# outliers [1]       : 22, 0.11%
# low counts [2]     : 9367, 46%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

cold_results = as.data.frame(cold_results) %>%
  arrange(padj)

write.csv(cold_results, "Oculina/cold_results.csv")

## Heat 
heat_results = results(dds, alpha = 0.05, contrast = c("Treatment", "Heat", "Control"))
heat_summary = summary(heat_results)

# out of 20488 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 125, 0.61%
# LFC < 0 (down)     : 207, 1%
# outliers [1]       : 22, 0.11%
# low counts [2]     : 15601, 76%
# (mean count < 5)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

heat_results = as.data.frame(heat_results) %>%
  arrange(padj)

write.csv(heat_results, "Oculina/heat_results.csv")

save.image("Oculina/temperature_DDS.RData")
