## Differential expression in Astrangia poculata
## Symbiont state
## Wald

library(DESeq2)
library(tidyverse)

counts = read.csv("Astrangia/outlier_removed_host_counts.csv", row.names = 1)
expDesign = read.csv("Astrangia/outlier_removed_expDesign.csv") %>%
  select(-X)

# Symbiont state across all treatments

dds = DESeqDataSetFromMatrix(countData = counts,
                             colData = expDesign,
                             design = ~ Treatment + Sym.Status)
dds = DESeq(dds)

sym_results = results(dds, alpha = 0.05, contrast = c("Sym.Status", "Brown", "White"))
sym_summary = summary(sym_results)

# out of 33650 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 79, 0.23%
# LFC < 0 (down)     : 63, 0.19%
# outliers [1]       : 20, 0.059%
# low counts [2]     : 22426, 67%
# (mean count < 2)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results


sym_results = as.data.frame(sym_results) %>%
  arrange(padj)

write.csv(sym_results, "Astrangia/sym_allTreatments_results.csv")



as_expDesign_ORv2 = as_expDesign_OR %>%
  mutate(Treatment_by_Sym.Status = paste(Treatment,Sym.Status, sep = "_"))

as_dds.sym = DESeqDataSetFromMatrix(countData = as_counts_OR, 
                                    colData = as_expDesign_ORv2,
                                    design = ~ Treatment_by_Sym.Status)

as_dds.sym = DESeq(as_dds.sym)
as_cold.sym_results = results(as_dds.sym, alpha = 0.05, contrast = c("Treatment_by_Sym.Status", "Cold_Brown", "Cold_White"))
as_cold.sym_summary = summary(as_cold.sym_results)