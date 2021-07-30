
### Colour Analyses 


## Libraries
library(cowplot)
library(dplyr)
library(readr)
library(reshape2)
library(ggplot2)
library(lavaan)
library(plyr)
library(readr)
library(caTools)
library(bitops)


packages <- c("ggplot2", "dplyr", "lavaan", "plyr", "cowplot", "rmarkdown", 
              "readr", "caTools", "bitops")

## WD
setwd("~/Documents/Classes/marine_phys/data/Colour_analyses/")

## Data manipulation

colour = read.csv("Colour_Analysis.csv", header = T) %>% as_tibble() %>%
  select(-Green, -Blue, - Day, -X, -X.1)

colour$coral_id = colour$Genotype %>% 
  substring(1,2) %>%
  substring(2,2)

m.colour = melt(colour) %>% as_tibble() %>%
  select(-variable, -Genotype) %>%
  rename(value = "Red_intensity") %>%
  rename(coral_id = "Genotype")

write.csv(m.colour, "mod_colour.csv", row.names = F)
colour$temperature = as.factor(colour$temperature)



write.csv(m.colour, "melted_colour.csv", row.names = F)
temp_colour = read.csv("Colour_Analysis_temperature.csv", header = T) %>%
  as_tibble() %>%
  select(-Day, -Green, -Blue) %>%
  rename(Red = "Red_channel_intensity")
temp_colour$Temperature = as.factor(temp_colour$Temperature)
temp_colour$Genotype = temp_colour$Genotype %>% 
  substring(1,2) %>%
  substring(2,2)




## Data visulization

ggplot(temp_colour, aes(x = Temperature, y = Red_channel_intensity, fill = Sym.State)) +
  geom_boxplot() +
 # geom_smooth(aes(group = Sym.State, color = Sym.State, fill = Sym.State)) + 
  scale_fill_manual(values=c("sienna3", "grey47"))+  
 # scale_color_manual(values=c("sienna3", "grey47"))+
  theme_cowplot()

ggplot(data = temp_colour, aes(x = Temperature, y = Red_channel_intensity)) +
  geom_jitter(aes(colour = Sym.State), alpha = 0.3) +
  geom_smooth(aes(group = Sym.State, color = Sym.State, fill = Sym.State)) + 
  scale_color_manual(values=c("sienna3", "grey47"))+
  scale_fill_manual(values=c("sienna3", "grey47"))+  
  theme_cowplot()


c1 = aov(Red_channel_intensity ~ Temperature*Sym.State, data = temp_colour)
summary(c1)
plot(c1)


#                       Df Sum Sq Mean Sq F value   Pr(>F)    
#Temperature             6   9833    1639   5.126 5.27e-05 ***
 # Sym.State               1 167496  167496 523.872  < 2e-16 ***
  #Temperature:Sym.State   6   3846     641   2.005   0.0653 .  
#Residuals             268  85687     320                     
---
TukeyHSD(c2)

c2 = lmer(Red_channel_intensity ~ Temperature*Sym.State + (1 | Genotype), data = temp_colour)


bwstat2=aov(Growth~Treatment, data=bw_final)
summary(bwstat2)
par(mfrow=c(2,2))
plot(bwstat2)



