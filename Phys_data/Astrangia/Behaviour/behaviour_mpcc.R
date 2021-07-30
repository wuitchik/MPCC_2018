### Behaviour analysis on MPCC astrangia 

library(dplyr)
library(reshape2)
library(ggplot2)

setwd("~/Documents/Experiments/Freezing_hot/Behaviour/")


behaviour = read.csv("behaviour_temperature.csv", header = T) %>% as_tibble() 
cols = c("hot" = "coral", "freezing" = "turquoise3")

behaviour$temperature = as.factor(behaviour$temperature)
behaviour$behaviour = as.factor(behaviour$behaviour)
behaviour$treatment_control = as.factor(b)

freezing = behaviour %>%
  filter(treatment == "freezing")

hot = behaviour %>%
  filter(experiment == "hot")

freezing = behaviour %>%
  filter(experiment == "freezing")
freezing$Day = as.factor(freezing$Day) 

freezing5 = freezing %>%
  filter(Day == c("Day1", "Day4", "Day7", "Day11", "Day15"))

ggplot(data = hot, aes(x = temperature, y = behaviour)) +
  geom_jitter(aes(colour = treatment), alpha = .9) +
  geom_smooth(aes(group = treatment, color = treatment, fill = treatment),method = lm) + 
  scale_color_manual(values = cols)+
  scale_fill_manual(values= cols)+  
  theme_cowplot()
ggsave("hot_plot.pdf")

library(RColorBrewer)

ggplot(data = hot, aes(x = Day)) +
  geom_bar(position = "fill", aes(fill = behaviour)) +
  facet_grid( . ~ treatment_control) + 
  scale_fill_brewer(palette = "Oranges")

ggplot(data = freezing, aes(x = reorder(Day))) +
  geom_bar(position = "fill", aes(fill = behaviour)) +
  facet_grid( . ~ treatment_control) + 
  scale_fill_brewer(palette = "GnBu")

ggplot(data = freezing5, aes(x = reorder(Day))) +
  geom_bar(position = "fill", aes(fill = behaviour)) +
  facet_grid( . ~ treatment_control) + 
  scale_fill_brewer(palette = "GnBu")

ggplot(data = freezing, aes(x = temperature, y = behaviour)) +
  geom_jitter(aes(colour = treatment), alpha = .9) +
  geom_smooth(aes(group = treatment, color = treatment, fill = treatment),method = lm) + 
  scale_color_manual(values = cols)+
  scale_fill_manual(values= cols)+  
  theme_cowplot()
ggsave("freezing_plot.pdf")


ggplot(data = behaviour, aes(x = temperature, y = behaviour)) +
  geom_jitter(aes(colour = treatment), alpha = .9) +
  geom_smooth(aes(group = treatment, color = treatment, fill = treatment)) + 
  scale_color_manual(values = cols)+
  scale_fill_manual(values= cols)+  
  theme_cowplot()
ggsave("behaviour_sym.jpeg")







