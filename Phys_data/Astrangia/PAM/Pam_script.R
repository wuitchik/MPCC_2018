###  PAM analyses 
# Script written by Daniel Wuitchik on November 2018

# First things first, load libraries
library(dplyr)
library(ggplot2)

# next up, we set the working directory
setwd("~/Documents/Classes/marine_phys/data/PAM/")

# load the data
pam = read.csv("PAM_practice.csv", header = TRUE) %>%
  as_tibble()

# organize my data frame
pam.organized = pam %>% 
  select(coral, genotype, Tank, Sym.Status, Treatment,PAM.Day.2.Average, PAM.Day.5.Average, PAM.Day.8.Average )

# I then organized this in excel to better read
pam258 = read.csv("PAM_258.csv", header = T) %>% as_tibble()
pam258$Day = as.factor(pam258$Day)
pam258 = as.data.frame(pam258)

pdf("pam_treatment.pdf")
ggplot(pam258, aes(x = Day, y = PAM, fill = Treatment)) +
  geom_boxplot() +
  facet_grid(. ~ Sym.Status) +
  scale_fill_manual(values=c("steelblue3", "seagreen3", "coral"))+
  theme_classic()
dev.off()









