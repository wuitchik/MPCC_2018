## Comparing to Dixon et al. stress modules (2020)

# libraries
library(tidyverse)
library(cowplot)

# Read in Dixon data
Dixon_MF = read.table("../Dixon_MF.csv", header = T)
Dixon_CC = read.table("../Dixon_CC.csv", header = T)
Dixon_BP = read.table("../Dixon_BP.csv", header = T)

# Heat x sym module
heat_xsym_BP = read.table("../../GO_DEGs/MWU_BP_sym_heat_results.csv", header = T)

# Terms in both sets
BP_goods=intersect(heat_xsym_BP$term,Dixon_BP$term)
data1=heat_xsym_BP[heat_xsym_BP$term %in% BP_goods,]
data2=Dixon_BP[Dixon_BP$term %in% BP_goods,]

# Combine them
bp_plot=merge(data1,data2,by="term")

ggplot(bp_plot, aes(delta.rank.x, delta.rank.y, label = name.y)) +
  scale_fill_distiller(palette = "Reds", direction = 1) +
  geom_hex() +
  labs( x = "Brown vs White in Heat Treatment",
        y = "Dixon") +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.75) +
  geom_smooth(method=lm, color="black", se =FALSE) +
  theme_cowplot()

