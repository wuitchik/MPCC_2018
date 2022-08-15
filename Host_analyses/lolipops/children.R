# Looking at the oxidative stress response
library(tidyverse)
library(GOfuncR)

# GO:0006979 -- Oxidative stress

children_oxidative = get_child_nodes("GO:0006979") 
children_oxidative_list = children_oxidative$child_go_id

# GO MWU
Heat_Sym_BP = read.table("../GO_DEGs/Astrangia/MWU_BP_sym_control_vs_hot_results_modified_pvalues.csv", header = T) %>%
  filter(term %in% children_oxidative_list) %>%
  mutate(Treatment = "Heat",
         Phenotype = "Brown")

Heat_apo_BP = read.table("../GO_DEGs/Astrangia/MWU_BP_apo_control_vs_hot_results_modified_pvalues.csv", header = T)%>%
  filter(term %in% children_oxidative_list)%>%
  mutate(Treatment = "Heat",
         Phenotype = "White")

Cold_Sym_BP = read.table("../GO_DEGs/Astrangia/MWU_BP_sym_control_vs_cold_results_modified_pvalues.csv", header = T) %>%
  filter(term %in% children_oxidative_list)%>%
  mutate(Treatment = "Cold",
         Phenotype = "Brown")

Cold_apo_BP = read.table("../GO_DEGs/Astrangia/MWU_BP_apo_control_vs_cold_results_modified_pvalues.csv", header = T)%>%
  filter(term %in% children_oxidative_list)%>%
  mutate(Treatment = "Cold",
         Phenotype = "White")

# Merged 

all = rbind(Heat_Sym_BP, Heat_apo_BP, Cold_Sym_BP, Cold_apo_BP) %>%
  mutate(Phenotype = as.factor(Phenotype)) %>%
  mutate(Treatment_phenotype = paste(Treatment, Phenotype, sep = " ")) %>%
  mutate(type=ifelse(Phenotype=="Brown","Highlighted","Normal")) 

test = aov(delta.rank ~ Treatment + Phenotype + Treatment:Phenotype, all)
summary(test)

TukeyHSD(test, conf.level = 0.95)

# plot

# Make colour pallet

# Points
shapes = c("White" = 1, "Brown" = 19)

# Colours
fill_cols = c("Cold White" = "white", "Cold Brown" = "#96bcfa","Heat Brown" = "#de835b","Heat White" = "white")
colour_cols = c("Cold White" = "#68a2ff", "Cold Brown" = "#68a2ff","Heat Brown" = "#ea6227","Heat White" = "#ea6227")

oxidative_plot = 
ggplot(data = all, aes(Treatment_phenotype, delta.rank)) +
  geom_point(aes(colour = Treatment_phenotype, shape = Phenotype), size = 3, stroke = 1.25, position=position_jitterdodge()) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0, colour = "black") +
  stat_summary(fun = "mean", size = 0.5, colour = "black") +
  #geom_boxplot(aes(fill = Treatment_phenotype, colour = Treatment_phenotype), lwd = 1, alpha =.5) +
  scale_fill_manual(values = fill_cols) +
  scale_colour_manual(values = colour_cols) +
  scale_shape_manual(values = shapes) +
  ylab("'Oxidative stress' (GO:0006979) \n Delta-Rank") +
  xlab("") +
  theme_classic() +
  theme(legend.position = "none") 


##### GO:0006950 -- response to stress

children_stress = get_child_nodes("GO:0006950") 
children_stress_list = children_stress$child_go_id

# GO MWU
Heat_Sym_BP = read.table("../GO_DEGs/Astrangia/MWU_BP_sym_control_vs_hot_results_modified_pvalues.csv", header = T) %>%
  filter(term %in% children_stress_list) %>%
  mutate(Treatment = "Heat",
         Phenotype = "Brown")

Heat_apo_BP = read.table("../GO_DEGs/Astrangia/MWU_BP_apo_control_vs_hot_results_modified_pvalues.csv", header = T)%>%
  filter(term %in% children_stress_list)%>%
  mutate(Treatment = "Heat",
         Phenotype = "White")

Cold_Sym_BP = read.table("../GO_DEGs/Astrangia/MWU_BP_sym_control_vs_cold_results_modified_pvalues.csv", header = T) %>%
  filter(term %in% children_stress_list)%>%
  mutate(Treatment = "Cold",
         Phenotype = "Brown")

Cold_apo_BP = read.table("../GO_DEGs/Astrangia/MWU_BP_apo_control_vs_cold_results_modified_pvalues.csv", header = T)%>%
  filter(term %in% children_stress_list)%>%
  mutate(Treatment = "Cold",
         Phenotype = "White")

# Merged 

all = rbind(Heat_Sym_BP, Heat_apo_BP, Cold_Sym_BP, Cold_apo_BP) %>%
  mutate(Phenotype = as.factor(Phenotype)) %>%
  mutate(Treatment_phenotype = paste(Treatment, Phenotype, sep = " ")) 

test = aov(delta.rank ~ Treatment + Phenotype + Treatment:Phenotype, all)
summary(test)

TukeyHSD(test, conf.level = 0.95)

# plot

# Make colour pallet

# Points
shapes = c("White" = 1, "Brown" = 19)

# Colours
fill_cols = c("Cold White" = "white", "Cold Brown" = "#96bcfa","Heat Brown" = "#de835b","Heat White" = "white")
colour_cols = c("Cold White" = "#68a2ff", "Cold Brown" = "#68a2ff","Heat Brown" = "#ea6227","Heat White" = "#ea6227")


ggplot(data = all, aes(Treatment, delta.rank)) +
  geom_point(aes(colour = Treatment_phenotype, shape = Phenotype), size = 3, stroke = 1.25, position=position_jitterdodge()) +
  geom_boxplot(aes(fill = Treatment_phenotype, colour = Treatment_phenotype), lwd = 1, alpha =.5) +
  scale_fill_manual(values = fill_cols) +
  scale_colour_manual(values = colour_cols) +
  scale_shape_manual(values = shapes) +
  ylab("GO Term Delta-Rank") +
  theme_classic() +
  theme(legend.position = "none") 




