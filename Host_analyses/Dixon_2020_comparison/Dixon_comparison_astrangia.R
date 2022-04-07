### Comparing the GO Delta Ranks to the stress module by Dixon et al. 2020
library(tidyverse)
library(patchwork)
library(cowplot)

# Read in Dixon's red module
Dixon_BP = read.table("Dixon_BP.csv", header = T)
Dixon_MF = read.table("Dixon_MF.csv", header = T)
Dixon_CC = read.table("Dixon_CC.csv", header = T)

# Read in Astrangia GO, Symbiont comparisons
Ast_Heat_Sym_BP = read.table("../GO_DEGs/Astrangia/GO_MWU_outputs/MWU_BP_heat_sym_results_modified_pvalues.csv", header = T)
Ast_Heat_Sym_MF = read.table("../GO_DEGs/Astrangia/GO_MWU_outputs/MWU_MF_heat_sym_results_modified_pvalues.csv", header = T)
Ast_Heat_Sym_CC = read.table("../GO_DEGs/Astrangia/GO_MWU_outputs/MWU_CC_heat_sym_results_modified_pvalues.csv", header = T)

#### Heat Treatment ####

# Biological Processes
heat_sym_goods = intersect(Ast_Heat_Sym_BP$term, Dixon_BP$term)
heat_sym = Ast_Heat_Sym_BP[Ast_Heat_Sym_BP$term %in% heat_sym_goods,]
Dixon_BP_set = Dixon_BP[Dixon_BP$term %in% heat_sym_goods,]

# Combine them
heat_sym_BP_data = merge(heat_sym, Dixon_BP_set, by="term")

heat_sym_BP_plot = ggplot(heat_sym_BP_data, aes(delta.rank.x, delta.rank.y, label = name.y)) +
  scale_fill_distiller(palette = "Reds", direction = 1) +
  geom_hex() +
  labs( x = "Sym State",
        y = "Dixon") +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.75) +
  geom_smooth(method=lm, color="black", se =FALSE) +
  theme_cowplot()

# Molecular Function
heat_sym_goods_MF = intersect(Ast_Heat_Sym_MF$term, Dixon_MF$term)
heat_sym_MF = Ast_Heat_Sym_MF[Ast_Heat_Sym_MF$term %in% heat_sym_goods_MF,]
Dixon_MF_set = Dixon_MF[Dixon_MF$term %in% heat_sym_goods_MF,]

# Combine them
heat_sym_MF_data = merge(heat_sym_MF, Dixon_MF_set, by="term")

heat_sym_MF_plot = ggplot(heat_sym_MF_data, aes(delta.rank.x, delta.rank.y, label = name.y)) +
  scale_fill_distiller(palette = "Reds", direction = 1) +
  geom_hex() +
  labs( x = "Sym State",
        y = "Dixon") +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.75) +
  geom_smooth(method=lm, color="black", se =FALSE) +
  theme_cowplot()

# Cellular Component
heat_sym_goods_CC = intersect(Ast_Heat_Sym_CC$term, Dixon_CC$term)
heat_sym_CC = Ast_Heat_Sym_CC[Ast_Heat_Sym_CC$term %in% heat_sym_goods_CC,]
Dixon_CC_set = Dixon_CC[Dixon_CC$term %in% heat_sym_goods_CC,]

# Combine them
heat_sym_CC_data = merge(heat_sym_CC, Dixon_CC_set, by="term")

heat_sym_CC_plot = ggplot(heat_sym_CC_data, aes(delta.rank.x, delta.rank.y, label = name.y)) +
  scale_fill_distiller(palette = "Reds", direction = 1) +
  geom_hex() +
  labs( x = "Sym State",
        y = "Dixon") +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.75) +
  geom_smooth(method=lm, color="black", se =FALSE) +
  theme_cowplot()


#### Control Treatment ####
# Read in Astrangia GO, Symbiont comparisons
Ast_control_Sym_BP = read.table("../GO_DEGs/Astrangia/GO_MWU_outputs/MWU_BP_sym_control_results_modified_pvalues.csv", header = T)
Ast_control_Sym_MF = read.table("../GO_DEGs/Astrangia/GO_MWU_outputs/MWU_MF_sym_control_results_modified_pvalues.csv", header = T)
Ast_control_Sym_CC = read.table("../GO_DEGs/Astrangia/GO_MWU_outputs/MWU_CC_sym_control_results_modified_pvalues.csv", header = T)

# Biological Processes
control_sym_goods = intersect(Ast_control_Sym_BP$term, Dixon_BP$term)
control_sym = Ast_control_Sym_BP[Ast_control_Sym_BP$term %in% control_sym_goods,]
Dixon_BP_set = Dixon_BP[Dixon_BP$term %in% control_sym_goods,]

# Combine them
control_sym_BP_data = merge(control_sym, Dixon_BP_set, by="term")

control_sym_BP_plot = ggplot(control_sym_BP_data, aes(delta.rank.x, delta.rank.y, label = name.y)) +
  scale_fill_distiller(palette = "Greys", direction = 1) +
  geom_hex() +
  labs( x = "Sym State",
        y = "Dixon") +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.75) +
  geom_smooth(method=lm, color="black", se =FALSE) +
  theme_cowplot()

# Molecular Function
control_sym_goods_MF = intersect(Ast_control_Sym_MF$term, Dixon_MF$term)
control_sym_MF = Ast_control_Sym_MF[Ast_control_Sym_MF$term %in% control_sym_goods_MF,]
Dixon_MF_set = Dixon_MF[Dixon_MF$term %in% control_sym_goods_MF,]

# Combine them
control_sym_MF_data = merge(control_sym_MF, Dixon_MF_set, by="term")

control_sym_MF_plot = ggplot(control_sym_MF_data, aes(delta.rank.x, delta.rank.y, label = name.y)) +
  scale_fill_distiller(palette = "Greys", direction = 1) +
  geom_hex() +
  labs( x = "Sym State",
        y = "Dixon") +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.75) +
  geom_smooth(method=lm, color="black", se =FALSE) +
  theme_cowplot()

# Cellular Component
control_sym_goods_CC = intersect(Ast_control_Sym_CC$term, Dixon_CC$term)
control_sym_CC = Ast_control_Sym_CC[Ast_control_Sym_CC$term %in% control_sym_goods_CC,]
Dixon_CC_set = Dixon_CC[Dixon_CC$term %in% control_sym_goods_CC,]

# Combine them
control_sym_CC_data = merge(control_sym_CC, Dixon_CC_set, by="term")

control_sym_CC_plot = ggplot(control_sym_CC_data, aes(delta.rank.x, delta.rank.y, label = name.y)) +
  scale_fill_distiller(palette = "Greys", direction = 1) +
  geom_hex() +
  labs( x = "Sym State",
        y = "Dixon") +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.75) +
  geom_smooth(method=lm, color="black", se =FALSE) +
  theme_cowplot()


#### Cold Treatment ####
# Read in Astrangia GO, Symbiont comparisons
Ast_Cold_Sym_BP = read.table("../GO_DEGs/Astrangia/GO_MWU_outputs/MWU_BP_cold_sym_results_modified_pvalues.csv", header = T)
Ast_Cold_Sym_MF = read.table("../GO_DEGs/Astrangia/GO_MWU_outputs/MWU_MF_cold_sym_results_modified_pvalues.csv", header = T)
Ast_Cold_Sym_CC = read.table("../GO_DEGs/Astrangia/GO_MWU_outputs/MWU_CC_cold_sym_results_modified_pvalues.csv", header = T)

# Biological Processes
Cold_sym_goods = intersect(Ast_Cold_Sym_BP$term, Dixon_BP$term)
Cold_sym = Ast_Cold_Sym_BP[Ast_Cold_Sym_BP$term %in% Cold_sym_goods,]
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
Cold_sym_goods_MF = intersect(Ast_Cold_Sym_MF$term, Dixon_MF$term)
Cold_sym_MF = Ast_Cold_Sym_MF[Ast_Cold_Sym_MF$term %in% Cold_sym_goods_MF,]
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
Cold_sym_goods_CC = intersect(Ast_Cold_Sym_CC$term, Dixon_CC$term)
Cold_sym_CC = Ast_Cold_Sym_CC[Ast_Cold_Sym_CC$term %in% Cold_sym_goods_CC,]
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



((heat_sym_BP_plot + heat_sym_MF_plot + heat_sym_CC_plot) /
    (control_sym_BP_plot + control_sym_MF_plot+ control_sym_CC_plot) /
    (Cold_sym_BP_plot + Cold_sym_MF_plot + Cold_sym_CC_plot) )

ggsave("Astrangia_Sym_state.pdf", 
       last_plot(),
       width = 10,
       height = 7,
       units = "in")
