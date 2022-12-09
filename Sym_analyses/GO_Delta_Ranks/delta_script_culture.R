## Delta ranks between cold and hot analyses for symbionts in host.

#### Libraries and files ####
library(tidyverse)
library(plotly)
library(viridis)
library(scales)

setwd("/Users/hannahaichelman/Documents/BU/Host_Buffering/MPCC_2018/Sym_analyses/GO_Delta_Ranks/")

# Load in both hot and cold comparisons for symbionts in host
SymInHost_Heat_BP = read.table("/Users/hannahaichelman/Documents/BU/Host_Buffering/MPCC_2018/Sym_analyses/GO_Analyses/MWU_BP_oculina_hot_syminhost_GO.csv", header = T)
SymInHost_Cold_BP = read.table("/Users/hannahaichelman/Documents/BU/Host_Buffering/MPCC_2018/Sym_analyses/GO_Analyses/MWU_BP_oculina_cold_syminhost_GO.csv", header = T)

SymInHost_Heat_MF = read.table("/Users/hannahaichelman/Documents/BU/Host_Buffering/MPCC_2018/Sym_analyses/GO_Analyses/MWU_MF_oculina_hot_syminhost_GO.csv", header = T)
SymInHost_Cold_MF = read.table("/Users/hannahaichelman/Documents/BU/Host_Buffering/MPCC_2018/Sym_analyses/GO_Analyses/MWU_MF_oculina_cold_syminhost_GO.csv", header = T)

SymInHost_Heat_CC = read.table("/Users/hannahaichelman/Documents/BU/Host_Buffering/MPCC_2018/Sym_analyses/GO_Analyses/MWU_CC_oculina_hot_syminhost_GO.csv", header = T)
SymInHost_Cold_CC = read.table("/Users/hannahaichelman/Documents/BU/Host_Buffering/MPCC_2018/Sym_analyses/GO_Analyses/MWU_CC_oculina_cold_syminhost_GO.csv", header = T)

# Load in both hot and cold comparisons for symbionts in culture
Culture_Heat_BP = read.table("/Users/hannahaichelman/Documents/BU/Host_Buffering/MPCC_2018/Sym_analyses/GO_Analyses/MWU_BP_culture_hot_GO.csv", header = T)
Culture_Cold_BP = read.table("/Users/hannahaichelman/Documents/BU/Host_Buffering/MPCC_2018/Sym_analyses/GO_Analyses/MWU_BP_culture_cold_GO.csv", header = T)

Culture_Heat_MF = read.table("/Users/hannahaichelman/Documents/BU/Host_Buffering/MPCC_2018/Sym_analyses/GO_Analyses/MWU_MF_culture_hot_GO.csv", header = T)
Culture_Cold_MF = read.table("/Users/hannahaichelman/Documents/BU/Host_Buffering/MPCC_2018/Sym_analyses/GO_Analyses/MWU_MF_culture_cold_GO.csv", header = T)

Culture_Heat_CC = read.table("/Users/hannahaichelman/Documents/BU/Host_Buffering/MPCC_2018/Sym_analyses/GO_Analyses/MWU_CC_culture_hot_GO.csv", header = T)
Culture_Cold_CC = read.table("/Users/hannahaichelman/Documents/BU/Host_Buffering/MPCC_2018/Sym_analyses/GO_Analyses/MWU_CC_culture_cold_GO.csv", header = T)

#### Symbiont in Host ####

# merge into one data set

Hot_v_Cold_BP = merge(SymInHost_Heat_BP, SymInHost_Cold_BP, by = "term") %>%
  mutate(divergence = delta.rank.x - delta.rank.y) %>%
  mutate(divergence = abs(divergence)) %>%
  arrange(desc(divergence))

Hot_v_Cold_MF = merge(SymInHost_Heat_MF, SymInHost_Cold_MF, by = "term") %>%
  mutate(divergence = delta.rank.x - delta.rank.y) %>%
  mutate(divergence = abs(divergence)) %>%
  arrange(desc(divergence))

Hot_v_Cold_CC = merge(SymInHost_Heat_CC, SymInHost_Cold_CC, by = "term") %>%
  mutate(divergence = delta.rank.x - delta.rank.y) %>%
  mutate(divergence = abs(divergence)) %>%
  arrange(desc(divergence))

## experimental plot
library(wesanderson)
scale_colour_gradientn(colours = c("darkred", "orange", "yellow", "white"))
pal <- wes_palette("Zissou1", 20, type = "continuous")

# BP
HotvCold_BP_plot = 
  ggplot(data = Hot_v_Cold_BP, aes(x = delta.rank.x, y = delta.rank.y, colour = divergence, label = name.y)) +
  geom_point() +
  scale_fill_continuous() +
  labs(x = "Heat Delta-Rank",
       y = "Cold Delta-Rank",
       title = "Sym in Host Biological Process") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  #scale_color_viridis(option = "magma", direction = -1) +
  scale_color_gradientn(colours = pal) +
  theme_bw()

ggplotly(HotvCold_BP_plot)

# MF
HotvCold_MF_plot = 
  ggplot(data = Hot_v_Cold_MF, aes(x = delta.rank.x, y = delta.rank.y, colour = divergence, label = name.y)) +
  geom_point() +
  scale_fill_continuous() +
  labs(x = "Heat Delta-Rank",
       y = "Cold Delta-Rank",
       title = "Sym in Host Molecular Function") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_color_viridis(option = "G", direction = -1) +
  theme_bw()

ggplotly(HotvCold_MF_plot)


# CC
HotvCold_CC_plot = 
  ggplot(data = Hot_v_Cold_CC, aes(x = delta.rank.x, y = delta.rank.y, colour = divergence, label = name.y)) +
  geom_point() +
  scale_fill_continuous() +
  labs(x = "Heat Delta-Rank",
       y = "Cold Delta-Rank",
       title = "Sym in Host Cellular Component") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_color_viridis(option = "G", direction = -1) +
  theme_bw()

ggplotly(HotvCold_CC_plot)


#### Symbiont in Culture ####

# merge into one data set

Hot_v_Cold_BP_culture = merge(Culture_Heat_BP, Culture_Cold_BP, by = "term") %>%
  mutate(divergence = delta.rank.x - delta.rank.y) %>%
  mutate(divergence = abs(divergence)) %>%
  arrange(desc(divergence))

Hot_v_Cold_MF_culture = merge(Culture_Heat_MF, Culture_Cold_MF, by = "term") %>%
  mutate(divergence = delta.rank.x - delta.rank.y) %>%
  mutate(divergence = abs(divergence)) %>%
  arrange(desc(divergence))

Hot_v_Cold_CC_culture = merge(Culture_Heat_CC, Culture_Cold_CC, by = "term") %>%
  mutate(divergence = delta.rank.x - delta.rank.y) %>%
  mutate(divergence = abs(divergence)) %>%
  arrange(desc(divergence))


# BP
HotvCold_BP_culture_plot = 
  ggplot(data = Hot_v_Cold_BP_culture, aes(x = delta.rank.x, y = delta.rank.y, colour = divergence, label = name.y)) +
  geom_point() +
  scale_fill_continuous() +
  labs(x = "Heat Delta-Rank",
       y = "Cold Delta-Rank",
       title = "Culture Biological Process") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_color_viridis(option = "G", direction = -1) +
  theme_bw()

ggplotly(HotvCold_BP_culture_plot)

# MF
HotvCold_MF_culture_plot = 
  ggplot(data = Hot_v_Cold_MF_culture, aes(x = delta.rank.x, y = delta.rank.y, colour = divergence, label = name.y)) +
  geom_point() +
  scale_fill_continuous() +
  labs(x = "Heat Delta-Rank",
       y = "Cold Delta-Rank",
       title = "Culture Molecular Function") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_color_viridis(option = "G", direction = -1) +
  theme_bw()

ggplotly(HotvCold_MF_culture_plot)


# CC
HotvCold_CC_culture_plot = 
  ggplot(data = Hot_v_Cold_CC_culture, aes(x = delta.rank.x, y = delta.rank.y, colour = divergence, label = name.y)) +
  geom_point() +
  scale_fill_continuous() +
  labs(x = "Heat Delta-Rank",
       y = "Cold Delta-Rank",
       title = "Culture Cellular Component") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_color_viridis(option = "G", direction = -1) +
  theme_bw()

ggplotly(HotvCold_CC_culture_plot)


#### Heat Maps of Interesting GOs ####
## Symbionts in Host ##
# read in go.obo database to figure out the go category id's that we want to pull.

go.obo = read.delim("/Users/hannahaichelman/Documents/BU/Host_Buffering/MPCC_2018/Sym_analyses/Oculina_GO_Analyses/go.obo")

# these terms are visually pulled from the Cold GO tree outputs, all terms related to photosynthesis 
inds = c(which(go.obo$format.version..1.2 == "name: tetrapyrrole binding"), # MF
         which(go.obo$format.version..1.2 == "name: chlorophyll binding"), # MF
         which(go.obo$format.version..1.2 == "name: photosystem"), # CC
         which(go.obo$format.version..1.2 == "name: thylakoid membrane"), # CC
         which(go.obo$format.version..1.2 == "name: photosynthesis, light harvesting"), # BP
         which(go.obo$format.version..1.2 == "name: protein-chromophore linkage") # BP
         ) 

want.go.obo = go.obo %>% 
  slice(sort(unique(c(inds - 1)))) %>%
  mutate_at("format.version..1.2", str_replace, "id: ", "")

##pulling all of the immunity terms and making a data frame
rldpval = read.delim("/Users/hannahaichelman/Documents/BU/Host_Buffering/MPCC_2018/Sym_analyses/tables/oculina_cold_syminhost_results.txt", header = TRUE) %>%
  rownames_to_column(var = "gene")

iso2go_sym = read.delim("/Users/hannahaichelman/Documents/BU/Host_Buffering/MPCC_2018/Sym_analyses/Oculina_GO_Analyses/B_psygmophilum_isogroup_to_GOterm.tab", sep = "\t", header = FALSE) %>%
  dplyr::rename("gene" = "V1") %>%
  dplyr::rename("GO_ID" = "V2")
  
photo_genes = iso2go_sym %>%
  filter(grepl("GO:0009521|GO:0009765|GO:0016168|GO:0018298|GO:0042651|GO:0046906", GO_ID)) %>%
  left_join(rldpval)
head(photo_genes)
