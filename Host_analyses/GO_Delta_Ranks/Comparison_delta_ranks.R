## Comparing GO Delta Ranks

# load libraries
library(tidyverse)

# CHOOSE SET OF INPUT FILES TO LOOP RUN -----------------------------------
astrangia_data_path = "Astrangia/MWU_outputs/"   # path to the astrangia data
oculina_data_path = "Oculina/MWU_outputs/" # Oculina path

astrangia_results_files <- dir(astrangia_data_path, pattern = ".csv*") # get file names
oculina_results_files <- dir(oculina_data_path, pattern = ".csv*") # get file names
astrangia_names <- sub("\\.csv.*", "", astrangia_results_files)
oculina_names <- sub("\\.csv.*", "", oculina_results_files)

# Load all result files from folder
for(i in astrangia_names){
  filepath = file.path(astrangia_data_path,paste(i,".csv", sep="")) # sets up .csv files
  assign(i, read.table(filepath, header = T))
}

for(i in oculina_names){
  filepath = file.path(oculina_data_path,paste(i,".csv", sep="")) # sets up .csv files
  assign(i, read.table(filepath, header = T))
}

#### Biological Processes and sym state

BP_heat_sym_goods=intersect(astrangia_BP_heat_sym$term,oculina_BP_heat_sym$term)
astrangia = astrangia_BP_heat_sym[astrangia_BP_heat_sym$term %in% BP_heat_sym_goods,]
oculina = oculina_BP_heat_sym[oculina_BP_heat_sym$term %in% BP_heat_sym_goods,]

# Combine them
merged_data = merge(astrangia,oculina,by="term")


# This is to manually look for interesting go terms, and you can play with it in excel
BP_intersect = merged_data %>%
  filter(p.adj.x <0.1) %>%
  filter(p.adj.y < 0.1)

write.csv(BP_intersect, "BP_heat_sym_intersect.csv")

# Read back in your manipulated csv for those that you want to use as labels
BP_intersect = read.csv("BP_heat_sym_intersect.csv")

# Here is the actual plot, lots of it is redundant 

ggplot(BP_intersect, aes(delta.rank.x, delta.rank.y, label = name.y)) +
  geom_point() + 
  labs( x = "Astrangia",
        y = "Oculina") +
  labs(title = "BP Heat Sym State") +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.75) +
  theme_cowplot()

