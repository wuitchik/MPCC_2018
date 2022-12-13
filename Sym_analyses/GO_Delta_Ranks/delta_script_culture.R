## Delta ranks between cold and hot analyses for symbionts in host.

#### Libraries and files ####
library(tidyverse)
library(plotly)
library(viridis)
library(scales)

setwd("/Users/hannahaichelman/Documents/BU/Host_Buffering/MPCC_2018/Sym_analyses/GO_Delta_Ranks/")

# Load in both hot and cold comparisons for symbionts in host
SymInHost_Heat_BP = read.table("/Users/hannahaichelman/Documents/BU/Host_Buffering/MPCC_2018/Sym_analyses/Oculina_GO_Analyses/MWU_BP_oculina_hot_syminhost_GO.csv", header = T)
SymInHost_Cold_BP = read.table("/Users/hannahaichelman/Documents/BU/Host_Buffering/MPCC_2018/Sym_analyses/Oculina_GO_Analyses/MWU_BP_oculina_cold_syminhost_GO.csv", header = T)

SymInHost_Heat_MF = read.table("/Users/hannahaichelman/Documents/BU/Host_Buffering/MPCC_2018/Sym_analyses/Oculina_GO_Analyses/MWU_MF_oculina_hot_syminhost_GO.csv", header = T)
SymInHost_Cold_MF = read.table("/Users/hannahaichelman/Documents/BU/Host_Buffering/MPCC_2018/Sym_analyses/Oculina_GO_Analyses/MWU_MF_oculina_cold_syminhost_GO.csv", header = T)

SymInHost_Heat_CC = read.table("/Users/hannahaichelman/Documents/BU/Host_Buffering/MPCC_2018/Sym_analyses/Oculina_GO_Analyses/MWU_CC_oculina_hot_syminhost_GO.csv", header = T)
SymInHost_Cold_CC = read.table("/Users/hannahaichelman/Documents/BU/Host_Buffering/MPCC_2018/Sym_analyses/Oculina_GO_Analyses/MWU_CC_oculina_cold_syminhost_GO.csv", header = T)

# Load in both hot and cold comparisons for symbionts in culture
Culture_Heat_BP = read.table("/Users/hannahaichelman/Documents/BU/Host_Buffering/MPCC_2018/Sym_analyses/Oculina_GO_Analyses/MWU_BP_culture_hot_GO.csv", header = T)
Culture_Cold_BP = read.table("/Users/hannahaichelman/Documents/BU/Host_Buffering/MPCC_2018/Sym_analyses/Oculina_GO_Analyses/MWU_BP_culture_cold_GO.csv", header = T)

Culture_Heat_MF = read.table("/Users/hannahaichelman/Documents/BU/Host_Buffering/MPCC_2018/Sym_analyses/Oculina_GO_Analyses/MWU_MF_culture_hot_GO.csv", header = T)
Culture_Cold_MF = read.table("/Users/hannahaichelman/Documents/BU/Host_Buffering/MPCC_2018/Sym_analyses/Oculina_GO_Analyses/MWU_MF_culture_cold_GO.csv", header = T)

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
  dplyr::slice(sort(unique(c(inds - 1)))) %>%
  mutate_at("format.version..1.2", str_replace, "id: ", "")

##pulling all of the photosynthesis terms and making a data frame
rldpval = read.csv("/Users/hannahaichelman/Documents/BU/Host_Buffering/MPCC_2018/Sym_analyses/tables/SymInHost_RLDandPVALS.csv", header = TRUE) %>%
  dplyr::rename("gene" = "X")

# read in iso2go file for the algae
iso2go_sym = read.delim("/Users/hannahaichelman/Documents/BU/Host_Buffering/MPCC_2018/Sym_analyses/Oculina_GO_Analyses/B_psygmophilum_isogroup_to_GOterm.tab", sep = "\t", header = FALSE) %>%
  dplyr::rename("gene" = "V1") %>%
  dplyr::rename("GO_ID" = "V2")
 
#pulling all of the photosynthesis terms from the iso2go 
photo_genes = iso2go_sym %>%
  filter(grepl("GO:0009521|GO:0009765|GO:0016168|GO:0018298|GO:0042651|GO:0046906", GO_ID)) %>%
  left_join(rldpval)
head(photo_genes)
str(photo_genes)

# add gene name to this df - trim down names and add to photo_genes
iso2gene = read.delim("/Users/hannahaichelman/Documents/BU/Host_Buffering/MPCC_TagSeq/References/B_psygmophilum_transcriptome/B_psygmophilum_isogroup_to_genename.tab", sep = "\t", header = FALSE) %>%
  dplyr::rename("gene" = "V1") %>%
  dplyr::rename("gene_name" = "V2")

iso2gene$gene_name = gsub("OS=.*", "", iso2gene$gene_name)

#iso2gene$gene_name = substr(iso2gene$gene_name, 0,40)
head(iso2gene)
str(iso2gene)

photo_genes_anno = photo_genes %>%
  left_join(iso2gene)
head(photo_genes_anno)


p.val = 0.10 # raw p-value for GO enriched

# filter based on p-values from deseq results
conds=photo_genes_anno[photo_genes_anno$pval.cold.syminhost<=p.val & !is.na(photo_genes_anno$pval.cold.syminhost),]
length(conds[,1])
#6
head(conds)

# add annotation as row names and remove p-values and identifying info from dataframe
row.names(conds)=conds$gene_name
exp = conds[, c(3:23)]
head(exp)

means=apply(exp,1,mean) # calculate means of rows
explc=exp-means # subtracting them
head(explc)
# this explc object is what we can use to make our heatmap

# now make heatmap
# set color palette
col0=colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)

# make treatment data frame
treatment = as.factor(sapply(strsplit(colnames(explc), split = "_"), "[[", 2)) %>%
  revalue(c("C" = "control", "F" = "cold", "H" = "heat"))
expDesign = data.frame(colnames(explc), treatment)
expDesign = expDesign %>%
  column_to_rownames(var = "colnames.explc.")

my_colour = list(treatment = c(cold = "#74c476", heat = "#fd8d3c", control = "#a6611a"))

#heat map of all photosynthesis genes
library(pheatmap)

pdf("/Users/hannahaichelman/Documents/BU/Host_Buffering/MPCC_2018/Sym_analyses/plots/syminhost_heatmap_photosynthesis.pdf", height = 6, width = 9, onefile = F)
pheatmap(explc, cluster_cols = TRUE, scale = "row", color = col0, annotation_col = expDesign, annotation_colors = my_colour, show_rownames = TRUE, show_colnames = FALSE, border_color = "NA")
dev.off()