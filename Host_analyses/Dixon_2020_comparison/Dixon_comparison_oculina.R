### Comparing the GO Delta Ranks to the stress module by Dixon et al. 2020
library(tidyverse)
library(patchwork)
library(cowplot)

# Read in Dixon's red module
Dixon_BP = read.table("Dixon_BP.csv", header = T)
Dixon_MF = read.table("Dixon_MF.csv", header = T)
Dixon_CC = read.table("Dixon_CC.csv", header = T)

# Read in Ocurangia GO, Symbiont comparisons
Ocu_Heat_Sym_BP = read.table("../GO_DEGs/Oculina/MWU_BP_sym_control_vs_hot_results_modified_pvalues.csv", header = T)
Ocu_Heat_Sym_MF = read.table("../GO_DEGs/Oculina/MWU_MF_sym_control_vs_hot_results_modified_pvalues.csv", header = T)
Ocu_Heat_Sym_CC = read.table("../GO_DEGs/Oculina/MWU_CC_sym_control_vs_hot_results_modified_pvalues.csv", header = T)

Ocu_Heat_apo_BP = read.table("../GO_DEGs/Oculina/MWU_BP_apo_control_vs_hot_results_modified_pvalues.csv", header = T)
Ocu_Heat_apo_MF = read.table("../GO_DEGs/Oculina/MWU_MF_apo_control_vs_hot_results_modified_pvalues.csv", header = T)
Ocu_Heat_apo_CC = read.table("../GO_DEGs/Oculina/MWU_CC_apo_control_vs_hot_results_modified_pvalues.csv", header = T)


#### Heat Treatment ####


#### Brown Phenotype 
# Biological Processes
heat_sym_goods = intersect(Ocu_Heat_Sym_BP$term, Dixon_BP$term)
heat_sym = Ocu_Heat_Sym_BP[Ocu_Heat_Sym_BP$term %in% heat_sym_goods,]
Dixon_BP_set = Dixon_BP[Dixon_BP$term %in% heat_sym_goods,]

# Combine them
heat_sym_BP_data = merge(heat_sym, Dixon_BP_set, by="term")

heat_sym_BP_plot = ggplot(heat_sym_BP_data, aes(delta.rank.x, delta.rank.y, label = name.y)) +
  scale_fill_distiller(palette = "Reds", direction = 1) +
  geom_hex() +
  labs( x = "Brown Phenotype",
        y = "Dixon") +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.75) +
  geom_smooth(method=lm, color="black", se =FALSE) +
  theme_cowplot()

# Molecular Function
heat_sym_goods_MF = intersect(Ocu_Heat_Sym_MF$term, Dixon_MF$term)
heat_sym_MF = Ocu_Heat_Sym_MF[Ocu_Heat_Sym_MF$term %in% heat_sym_goods_MF,]
Dixon_MF_set = Dixon_MF[Dixon_MF$term %in% heat_sym_goods_MF,]

# Combine them
heat_sym_MF_data = merge(heat_sym_MF, Dixon_MF_set, by="term")

heat_sym_MF_plot = ggplot(heat_sym_MF_data, aes(delta.rank.x, delta.rank.y, label = name.y)) +
  scale_fill_distiller(palette = "Reds", direction = 1) +
  geom_hex() +
  labs( x = "Brown Phenotype",
        y = "Dixon") +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.75) +
  geom_smooth(method=lm, color="black", se =FALSE) +
  theme_cowplot()

# Cellular Component
heat_sym_goods_CC = intersect(Ocu_Heat_Sym_CC$term, Dixon_CC$term)
heat_sym_CC = Ocu_Heat_Sym_CC[Ocu_Heat_Sym_CC$term %in% heat_sym_goods_CC,]
Dixon_CC_set = Dixon_CC[Dixon_CC$term %in% heat_sym_goods_CC,]

# Combine them
heat_sym_CC_data = merge(heat_sym_CC, Dixon_CC_set, by="term")

heat_sym_CC_plot = ggplot(heat_sym_CC_data, aes(delta.rank.x, delta.rank.y, label = name.y)) +
  scale_fill_distiller(palette = "Reds", direction = 1) +
  geom_hex() +
  labs( x = "Brown Phenotype",
        y = "Dixon") +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.75) +
  geom_smooth(method=lm, color="black", se =FALSE) +
  theme_cowplot()


#### White Phenotype 
# Biological Processes
heat_apo_goods = intersect(Ocu_Heat_apo_BP$term, Dixon_BP$term)
heat_apo = Ocu_Heat_apo_BP[Ocu_Heat_apo_BP$term %in% heat_apo_goods,]
Dixon_BP_set = Dixon_BP[Dixon_BP$term %in% heat_apo_goods,]

# Combine them
heat_apo_BP_data = merge(heat_apo, Dixon_BP_set, by="term")

heat_apo_BP_plot = ggplot(heat_apo_BP_data, aes(delta.rank.x, delta.rank.y, label = name.y)) +
  scale_fill_distiller(palette = "Reds", direction = 1) +
  geom_hex() +
  labs( x = "White Phenotype",
        y = "Dixon") +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.75) +
  geom_smooth(method=lm, color="black", se =FALSE) +
  theme_cowplot()

# Molecular Function
heat_apo_goods_MF = intersect(Ocu_Heat_apo_MF$term, Dixon_MF$term)
heat_apo_MF = Ocu_Heat_apo_MF[Ocu_Heat_apo_MF$term %in% heat_apo_goods_MF,]
Dixon_MF_set = Dixon_MF[Dixon_MF$term %in% heat_apo_goods_MF,]

# Combine them
heat_apo_MF_data = merge(heat_apo_MF, Dixon_MF_set, by="term")

heat_apo_MF_plot = ggplot(heat_apo_MF_data, aes(delta.rank.x, delta.rank.y, label = name.y)) +
  scale_fill_distiller(palette = "Reds", direction = 1) +
  geom_hex() +
  labs( x = "White Phenotype",
        y = "Dixon") +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.75) +
  geom_smooth(method=lm, color="black", se =FALSE) +
  theme_cowplot()

# Cellular Component
heat_apo_goods_CC = intersect(Ocu_Heat_apo_CC$term, Dixon_CC$term)
heat_apo_CC = Ocu_Heat_apo_CC[Ocu_Heat_apo_CC$term %in% heat_apo_goods_CC,]
Dixon_CC_set = Dixon_CC[Dixon_CC$term %in% heat_apo_goods_CC,]

# Combine them
heat_apo_CC_data = merge(heat_apo_CC, Dixon_CC_set, by="term")

heat_apo_CC_plot = ggplot(heat_apo_CC_data, aes(delta.rank.x, delta.rank.y, label = name.y)) +
  scale_fill_distiller(palette = "Reds", direction = 1) +
  geom_hex() +
  labs( x = "White Phenotype",
        y = "Dixon") +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.75) +
  geom_smooth(method=lm, color="black", se =FALSE) +
  theme_cowplot()

((heat_sym_BP_plot + heat_sym_MF_plot + heat_sym_CC_plot) /
    (heat_apo_BP_plot + heat_apo_MF_plot + heat_apo_CC_plot))

ggsave("Oculina_heat_within_sym.pdf", 
       last_plot(),
       width = 10,
       height = 7,
       units = "in")

#### Cold Treatment ####
# Read in Ocurangia GO, Symbiont comparisons
Ocu_Cold_Sym_BP = read.table("../GO_DEGs/Oculina/MWU_BP_sym_control_vs_cold_results_modified_pvalues.csv", header = T)
Ocu_Cold_Sym_MF = read.table("../GO_DEGs/Oculina/MWU_BP_sym_control_vs_cold_results_modified_pvalues.csv", header = T)
Ocu_Cold_Sym_CC = read.table("../GO_DEGs/Oculina/MWU_BP_sym_control_vs_cold_results_modified_pvalues.csv", header = T)

Ocu_Cold_apo_BP = read.table("../GO_DEGs/Oculina/MWU_BP_apo_control_vs_cold_results_modified_pvalues.csv", header = T)
Ocu_Cold_apo_MF = read.table("../GO_DEGs/Oculina/MWU_BP_apo_control_vs_cold_results_modified_pvalues.csv", header = T)
Ocu_Cold_apo_CC = read.table("../GO_DEGs/Oculina/MWU_BP_apo_control_vs_cold_results_modified_pvalues.csv", header = T)

#### Brown Phenotype 

# Biological Processes
Cold_sym_goods = intersect(Ocu_Cold_Sym_BP$term, Dixon_BP$term)
Cold_sym = Ocu_Cold_Sym_BP[Ocu_Cold_Sym_BP$term %in% Cold_sym_goods,]
Dixon_BP_set = Dixon_BP[Dixon_BP$term %in% Cold_sym_goods,]

# Combine them
Cold_sym_BP_data = merge(Cold_sym, Dixon_BP_set, by="term")

Cold_sym_BP_plot = ggplot(Cold_sym_BP_data, aes(delta.rank.x, delta.rank.y, label = name.y)) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  geom_hex() +
  labs( x = "Sym State",
        y = "Dixon") +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.75) +
  geom_smooth(method=lm, color="black", se =FALSE) +
  theme_cowplot()

# Molecular Function
Cold_sym_goods_MF = intersect(Ocu_Cold_Sym_MF$term, Dixon_MF$term)
Cold_sym_MF = Ocu_Cold_Sym_MF[Ocu_Cold_Sym_MF$term %in% Cold_sym_goods_MF,]
Dixon_MF_set = Dixon_MF[Dixon_MF$term %in% Cold_sym_goods_MF,]

# Combine them
Cold_sym_MF_data = merge(Cold_sym_MF, Dixon_MF_set, by="term")

Cold_sym_MF_plot = ggplot(Cold_sym_MF_data, aes(delta.rank.x, delta.rank.y, label = name.y)) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  geom_hex() +
  labs( x = "Sym State",
        y = "Dixon") +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.75) +
  geom_smooth(method=lm, color="black", se =FALSE) +
  theme_cowplot()

# Cellular Component
Cold_sym_goods_CC = intersect(Ocu_Cold_Sym_CC$term, Dixon_CC$term)
Cold_sym_CC = Ocu_Cold_Sym_CC[Ocu_Cold_Sym_CC$term %in% Cold_sym_goods_CC,]
Dixon_CC_set = Dixon_CC[Dixon_CC$term %in% Cold_sym_goods_CC,]

# Combine them
Cold_sym_CC_data = merge(Cold_sym_CC, Dixon_CC_set, by="term")

Cold_sym_CC_plot = ggplot(Cold_sym_CC_data, aes(delta.rank.x, delta.rank.y, label = name.y)) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  geom_hex() +
  labs( x = "Sym State",
        y = "Dixon") +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.75) +
  geom_smooth(method=lm, color="black", se =FALSE) +
  theme_cowplot()

#### White phenotype

# Biological Processes
Cold_apo_goods = intersect(Ocu_Cold_apo_BP$term, Dixon_BP$term)
Cold_apo = Ocu_Cold_apo_BP[Ocu_Cold_apo_BP$term %in% Cold_apo_goods,]
Dixon_BP_set = Dixon_BP[Dixon_BP$term %in% Cold_apo_goods,]

# Combine them
Cold_apo_BP_data = merge(Cold_apo, Dixon_BP_set, by="term")

Cold_apo_BP_plot = ggplot(Cold_apo_BP_data, aes(delta.rank.x, delta.rank.y, label = name.y)) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  geom_hex() +
  labs( x = "apo State",
        y = "Dixon") +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.75) +
  geom_smooth(method=lm, color="black", se =FALSE) +
  theme_cowplot()

# Molecular Function
Cold_apo_goods_MF = intersect(Ocu_Cold_apo_MF$term, Dixon_MF$term)
Cold_apo_MF = Ocu_Cold_apo_MF[Ocu_Cold_apo_MF$term %in% Cold_apo_goods_MF,]
Dixon_MF_set = Dixon_MF[Dixon_MF$term %in% Cold_apo_goods_MF,]

# Combine them
Cold_apo_MF_data = merge(Cold_apo_MF, Dixon_MF_set, by="term")

Cold_apo_MF_plot = ggplot(Cold_apo_MF_data, aes(delta.rank.x, delta.rank.y, label = name.y)) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  geom_hex() +
  labs( x = "apo State",
        y = "Dixon") +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.75) +
  geom_smooth(method=lm, color="black", se =FALSE) +
  theme_cowplot()

# Cellular Component
Cold_apo_goods_CC = intersect(Ocu_Cold_apo_CC$term, Dixon_CC$term)
Cold_apo_CC = Ocu_Cold_apo_CC[Ocu_Cold_apo_CC$term %in% Cold_apo_goods_CC,]
Dixon_CC_set = Dixon_CC[Dixon_CC$term %in% Cold_apo_goods_CC,]

# Combine them
Cold_apo_CC_data = merge(Cold_apo_CC, Dixon_CC_set, by="term")

Cold_apo_CC_plot = ggplot(Cold_apo_CC_data, aes(delta.rank.x, delta.rank.y, label = name.y)) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  geom_hex() +
  labs( x = "apo State",
        y = "Dixon") +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.75) +
  geom_smooth(method=lm, color="black", se =FALSE) +
  theme_cowplot()

((Cold_sym_BP_plot + Cold_sym_MF_plot + Cold_sym_CC_plot) /
    (Cold_apo_BP_plot + Cold_apo_MF_plot + Cold_apo_CC_plot))

ggsave("Oculina_cold_within_sym.pdf", 
       last_plot(),
       width = 10,
       height = 7,
       units = "in")
