## Differential expression in Oculina arbuscula
## Symbiont state

library(DESeq2)
library(tidyverse)

counts = read.csv("Oculina/outlier_removed_host_counts.csv", row.names = 1)
expDesign = read.csv("Oculina/outlier_removed_expDesign.csv")

# Symbiont state across all treatments

dds = DESeqDataSetFromMatrix(countData = counts,
                             colData = expDesign,
                             design = ~ Treatment + Sym.Status)
dds = DESeq(dds)

sym_results = results(dds, alpha = 0.05, contrast = c("Sym.Status", "Brown", "White"))
sym_summary = summary(sym_results)

# out of 20488 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 58, 0.28%
# LFC < 0 (down)     : 56, 0.27%
# outliers [1]       : 22, 0.11%
# low counts [2]     : 10535, 51%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results


sym_results = as.data.frame(sym_results) %>%
  arrange(padj)

write.csv(sym_results, "Oculina/sym_allTreatments_results.csv")

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

# out of 20488 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 0, 0%
# LFC < 0 (down)     : 1, 0.0049%
# outliers [1]       : 16, 0.078%
# low counts [2]     : 0, 0%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results


control_sym_results = as.data.frame(control_sym_results) %>%
  arrange(padj)

write.csv(control_sym_results, "Oculina/sym_control_results.csv")

# Symbiont state in heat

heat_sym_results = results(dds.sym,
                           alpha = 0.05,
                           contrast = c("Treatment_by_Sym.Status", "Heat_Brown", "Heat_White"))
heat_sym_summary = summary(heat_sym_results)

# out of 20488 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 0, 0%
# LFC < 0 (down)     : 1, 0.0049%
# outliers [1]       : 16, 0.078%
# low counts [2]     : 0, 0%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

heat_sym_results = as.data.frame(heat_sym_results) %>%
  arrange(padj)

write.csv(heat_sym_results, "Oculina/heat_sym_results.csv")

# Symbiont state in cold

cold_sym_results = results(dds.sym,
                           alpha = 0.05,
                           contrast = c("Treatment_by_Sym.Status", "Cold_Brown", "Cold_White"))
cold_sym_summary = summary(cold_sym_results)

# out of 20488 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 2, 0.0098%
# LFC < 0 (down)     : 1, 0.0049%
# outliers [1]       : 16, 0.078%
# low counts [2]     : 0, 0%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

cold_sym_results = as.data.frame(cold_sym_results) %>%
  arrange(padj)

write.csv(cold_sym_results, "Oculina/cold_sym_results.csv")

save.image(file = "Oculina/Oculina_sym_state.RData")

