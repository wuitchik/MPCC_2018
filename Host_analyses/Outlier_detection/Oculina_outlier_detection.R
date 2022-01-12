### Here, we explore for potential outlier samples in Oculina arbuscula that result from poor sequencing or other such errors

library(DESeq)
library(arrayQualityMetrics)
library(tidyverse)

counts = read.csv("Oculina/oculina_counts.csv", row.names = 1) %>%
  rename_all(funs(str_replace_all(.,".counts.txt", "")))

expDesign = read.csv("Oculina/oculina_expDesign.csv") %>%
  mutate(Treatment = case_when(
    substr(Tank.ID,1,7) == "control" ~ "Control",
    substr(Tank.ID,1,4) == "cold" ~ "Cold", 
    substr(Tank.ID,1,4) == "heat" ~ "Heat")) %>%
  mutate(Sym.Status = case_when(
    Symbiotic.status == "Sym" ~ "Brown",
    Symbiotic.status == "Apo" ~ "White")) %>%
  dplyr::rename(Sample = Coral.ID) %>%
  dplyr::select(-Tank.ID, -Symbiotic.status, -X)

# Here I randomly select one clonal sample from genotypes X,Y,Z balanced by treatment

set.seed(1) # setting seed makes it so I consistently use the same random samples

x.genotype = expDesign %>%
  filter(Genotype == "X") %>%
  group_by(Treatment) %>%
  slice_sample(n=1)

y.genotype = expDesign %>%
  filter(Genotype == "Y") %>%
  group_by(Treatment) %>%
  slice_sample(n=1)

z.genotype = expDesign %>%
  filter(Genotype == "Z") %>%
  group_by(Treatment) %>%
  slice_sample(n=2)

excluded_list = c(x.genotype$Sample, y.genotype$Sample, z.genotype$Sample)

expDesign = expDesign %>%
  filter(!Sample %in% excluded_list)

write.csv(expDesign, "Oculina/clones_removed_expDesign.csv", row.names = F)

counts = counts %>%
  select(expDesign$Sample) %>%
  na.omit()

write.csv(counts, "Oculina/clones_removed_counts.csv", row.names = T)

# Using arrayQualityMetrics to look for outliers

real = newCountDataSet(counts, expDesign)
real = estimateSizeFactors(real)

cds = estimateDispersions(real,method="blind")
vsdBlind = varianceStabilizingTransformation(cds)
arrayQualityMetrics(vsdBlind,intgroup=c("Treatment"), force=TRUE, outdir = "Oculina/arrayQualityMetrics") 

# Visual inspection of Oculina/index.html reveals one outlier (R9) which was removed

counts = counts %>%
  select(-R9)

expDesign = expDesign %>%
  filter(Sample %in% colnames(counts)) 

write.csv(counts, "Oculina/outlier_removed_host_counts.csv")
write.csv(expDesign, "Oculina/outlier_removed_expDesign.csv")
