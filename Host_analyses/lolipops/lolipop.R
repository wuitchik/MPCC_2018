# lolipop script

library(tidyverse)
library(GOfuncR)

# first trying to get to parent GO terms
Heat_Sym_BP = read.table("../GO_DEGs/Astrangia/MWU_BP_sym_control_vs_hot_results_modified_pvalues.csv", header = T) %>%
  filter(p.adj < 0.05)
Heat_Sym_MF = read.table("../GO_DEGs/Astrangia/MWU_MF_sym_control_vs_hot_results_modified_pvalues.csv", header = T)
Heat_Sym_CC = read.table("../GO_DEGs/Astrangia/MWU_CC_sym_control_vs_hot_results_modified_pvalues.csv", header = T)

Heat_apo_BP = read.table("../GO_DEGs/Astrangia/MWU_BP_apo_control_vs_hot_results_modified_pvalues.csv", header = T)%>%
  filter(p.adj < 0.05)
Heat_apo_MF = read.table("../GO_DEGs/Astrangia/MWU_MF_apo_control_vs_hot_results_modified_pvalues.csv", header = T)
Heat_apo_CC = read.table("../GO_DEGs/Astrangia/MWU_CC_apo_control_vs_hot_results_modified_pvalues.csv", header = T)

# testing

Heat_Sym_BP_terms = Heat_Sym_BP$term

Heat_Sym_BP_parents = get_parent_nodes(Heat_Sym_BP_terms)


Heat_apo_BP_terms = Heat_apo_BP$term

Heat_apo_BP_parents = get_parent_nodes(Heat_apo_BP_terms)


write.csv(Heat_apo_BP_parents, "apo_test.csv")
write.csv(Heat_Sym_BP_parents, "sym_test.csv")
