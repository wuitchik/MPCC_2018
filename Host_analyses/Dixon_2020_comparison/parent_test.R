## playing with parent GO terms


library(GO.db)
library(tidyverse)
library(GOSim)



# Read in Dixon's red module
Dixon_BP = read.table("Dixon_BP.csv", header = T)
Dixon_MF = read.table("Dixon_MF.csv", header = T)
Dixon_CC = read.table("Dixon_CC.csv", header = T)

# Read in Astrangia GO, Symbiont comparisons
Heat_Sym_BP = read.table("../GO_DEGs/Astrangia/MWU_BP_sym_control_vs_hot_results_modified_pvalues.csv", header = T)
Heat_Sym_MF = read.table("../GO_DEGs/Astrangia/MWU_MF_sym_control_vs_hot_results_modified_pvalues.csv", header = T)
Heat_Sym_CC = read.table("../GO_DEGs/Astrangia/MWU_CC_sym_control_vs_hot_results_modified_pvalues.csv", header = T)

Heat_apo_BP = read.table("../GO_DEGs/Astrangia/MWU_BP_apo_control_vs_hot_results_modified_pvalues.csv", header = T)
Heat_apo_MF = read.table("../GO_DEGs/Astrangia/MWU_MF_apo_control_vs_hot_results_modified_pvalues.csv", header = T)
Heat_apo_CC = read.table("../GO_DEGs/Astrangia/MWU_CC_apo_control_vs_hot_results_modified_pvalues.csv", header = T)


terms = Heat_Sym_BP$term

terms <- c("GO:0008022","GO:0001071")

parents = as.list(GOMFANCESTOR[terms])

## Bimap interface:
# Convert the object to a list
xx <- as.list(GOBPPARENTS)
# Remove GO IDs that do not have any parent
xx <- xx[!is.na(xx)]
if(length(xx) > 0){
  # Get the children GO IDs for the first elents of xx
  goids <- xx[[1]]
  # Find out the GO terms for the first parent goid
  GOID(GOTERM[[goids[1]]])
  Term(GOTERM[[goids[1]]])
  Synonym(GOTERM[[goids[1]]])
  Secondary(GOTERM[[goids[1]]])
  Definition(GOTERM[[goids[1]]])
  Ontology(GOTERM[[goids[1]]])
}

