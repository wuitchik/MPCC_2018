## Differential expression in Astrangia poculata
## Symbiont state

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


## Symbiont state in control

expDesign = expDesign %>%
  mutate(Treatment_by_Sym.Status = paste(Treatment,Sym.Status, sep = "_"))

dds.sym = DESeqDataSetFromMatrix(countData = counts,
                                 colData = expDesign,
                                 design = ~ Treatment_by_Sym.Status)
dds.sym = DESeq(dds.sym)
control_sym_results = results(dds.sym,
                              alpha = 0.05,
                              contrast = c("Treatment_by_Sym.Status", "Control_Brown", "Control_White"))
control_sym_summary = summary(control_sym_results)

# out of 33652 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 0, 0%
# LFC < 0 (down)     : 0, 0%
# outliers [1]       : 12, 0.036%
# low counts [2]     : 2, 0.0059%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

control_sym_results = as.data.frame(control_sym_results) %>%
  arrange(padj)

write.csv(control_sym_results, "Astrangia/sym_control_results.csv")

# Symbiont state in heat

heat_sym_results = results(dds.sym,
                           alpha = 0.05,
                           contrast = c("Treatment_by_Sym.Status", "Heat_Brown", "Heat_White"))
heat_sym_summary = summary(heat_sym_results)

# out of 33652 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 91, 0.27%
# LFC < 0 (down)     : 8, 0.024%
# outliers [1]       : 12, 0.036%
# low counts [2]     : 23711, 70%
# (mean count < 3)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

heat_sym_results = as.data.frame(heat_sym_results) %>%
  arrange(padj)

write.csv(heat_sym_results, "Astrangia/heat_sym_results.csv")

# Symbiont state in cold

cold_sym_results = results(dds.sym,
                           alpha = 0.05,
                           contrast = c("Treatment_by_Sym.Status", "Cold_Brown", "Cold_White"))
cold_sym_summary = summary(cold_sym_results)

# out of 33652 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 1, 0.003%
# LFC < 0 (down)     : 0, 0%
# outliers [1]       : 12, 0.036%
# low counts [2]     : 2, 0.0059%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

cold_sym_results = as.data.frame(cold_sym_results) %>%
  arrange(padj)

write.csv(cold_sym_results, "Astrangia/cold_sym_results.csv")

save.image(file = "Astrangia/Astrangia_sym_state.RData")
