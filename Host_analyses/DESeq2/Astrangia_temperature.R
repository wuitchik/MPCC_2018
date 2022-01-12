## Differential expression in Astrangia poculata
## Cold results
## Wald test

library(DESeq2)
library(tidyverse)

counts = read.csv("Astrangia/outlier_removed_host_counts.csv", row.names = 1)
expDesign = read.csv("Astrangia/outlier_removed_expDesign.csv") %>%
  select(-X)

dds = DESeqDataSetFromMatrix(countData = as_counts_OR,
                             colData = as_expDesign_OR,
                             design = ~ Treatment + Sym.Status)
as_dds = DESeq(as_dds)

as_cold_results = results(as_dds, alpha = 0.05, contrast = c("Treatment", "Cold", "Control"))
as_cold_summary = summary(as_cold_results)
write.csv(as_cold_results, "Astrangia/cold_results.csv")