###BW Analysis 
library(dplyr)
library(tidyverse)
library(ggplot2)
library(reshape2)

Sample Treatment Genotype
1     AA2   control        A
2     AB3   control        B
3     AB4   control        B
4     AC2   control        C
5     AD2   control        D
6     AD4   control        D
7     AF8   control        F
8     AQ2   control        Q
9     AQ4   control        Q
10    AO4   control        O
11    AC3   control        C
12    AQ5   control        Q
13    AA1       hot        A
14    AB2       hot        B
15    AB5       hot        B
16    AB6       hot        B
17    AD6       hot        D
18    AK4       hot        K
19    AN1       hot        N
20    AO1       hot        O
21    AO2       hot        O
22    AQ3       hot        Q
23    AQ6       hot        Q


setwd("~/Documents/Boston University/Marine Semester/MPCC")

#load data
bw = read.csv("BW_calc1.csv", header = TRUE) %>% 
  as_tibble()
head(bw)

m.bw=melt(bw)
head(m.bw)

colnames(m.bw)= c("Nubbin","Treatment", "Sym", "Reglue", "Day_Comparison","Growth")

head(m.bw)

bw_final = m.bw %>% as_tibble() %>%
  select(Nubbin, Treatment, Sym, Reglue, Day_Comparison, Growth) %>%
  filter(Reglue == "N")%>%
  filter(Day_Comparison == "FmI")

bw_final$genotype = bw_final$Nubbin %>% 
  substring(1,2) %>%
  substring(2,2)

colnames(bw_final)= c("Nubbin","Treatment", "Sym", "Reglue", "Day_Comparison","Growth", "Genotype")

head(bw_final)

bw_final$Genotype=as.factor(bw_final$Genotype)

head(bw_final)

ggplot(bw_final, aes(x = Day_Comparison, y = Growth, fill = Treatment)) +
  geom_boxplot() +
  #facet_grid(genotype ~ .) +
  scale_fill_manual(values=c("steelblue3", "seagreen3", "coral"))+
  theme_bw()+
  xlab("Treatment")+ylab("% Change in Calcification Rate")

ggplot(bw_final, aes(x = Sym, y = Growth, fill = Treatment)) +
  geom_boxplot() +
  #facet_grid(genotype ~ .) +
  scale_fill_manual(values=c("steelblue3", "seagreen3", "coral"))+
  theme_bw()+
  xlab("Symbiotic State")+ylab("% Change in Calcification Rate")

bw_final1 = m.bw %>% as_tibble() %>%
  select(Nubbin, Treatment, Sym, Reglue, Day_Comparison, Growth) %>%
  filter(Reglue == "N")%>%
  filter(Day_Comparison == "FmI")

ggplot(bw_final, aes(x = Genotype, y = Growth, fill = Treatment)) +
  geom_point() +
  #facet_grid(genotype ~ .) +
  scale_fill_manual(values=c("steelblue3", "seagreen3", "coral"))+
  theme_bw()+
  xlab("Symbiotic State")+ylab("% Change in Calcification Rate")

bwstat=aov(Growth~Day_Comparison, data=bw_final)
summary(bwstat)
par(mfrow=c(2,2))
plot(bwstat)

TukeyHSD(bwstat)


bwstat2=aov(Growth~Treatment, data=bw_final)
summary(bwstat2)
par(mfrow=c(2,2))
plot(bwstat2)

TukeyHSD(bwstat2)

bwstat3=aov(Growth~Genotype, data=bw_final)
summary(bwstat3)
par(mfrow=c(2,2))
plot(bwstat3)



bwstat3=aov(Growth~Treatment*Sym*Genotype, data=bw_final)
summary(bwstat3)
par(mfrow=c(2,2))
plot(bwstat3)

bwstat4=aov(Growth~Sym, data=bw_final)
summary(bwstat4)
par(mfrow=c(2,2))
plot(bwstat4)


bwstat5=aov(Growth~Genotype, data=bw_final)
summary(bwstat5)
par(mfrow=c(2,2))
plot(bwstat5)

bwstat5=aov(Growth~Genotype, data=bw_final)
summary(bwstat5)
par(mfrow=c(2,2))
plot(bwstat5)

bwstat6=aov(Growth~Treatment*Sym, data=bw_final)
summary(bwstat6)
par(mfrow=c(2,2))
plot(bwstat6)

TukeyHSD(bwstat5)

library(lme4)
test = glmer(formula = growth ~ Treatment + sym + ( 1 | genotype ), family = gaussian, data = bw_)
