# This script analyzes the temperature data from the Oculina arbuscula holobiont experiment and the Breviolum psygmophilum culture experiment. 
# Author: Hannah E Aichelman
# Last Updated: January 1, 2023

#### Load in and format data ####
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(Rmisc)

setwd("/Users/hannahaichelman/Documents/BU/Host_Buffering/")

host_exp = read.csv("/Users/hannahaichelman/Documents/BU/Host_Buffering/MPCC_TagSeq/Oculina_temperature.csv") %>%
  select(-X) %>%
  mutate(treatment = factor(treatment))
head(host_exp)
str(host_exp)

culture_exp = read.csv("/Users/hannahaichelman/Documents/BU/Host_Buffering/Sym_TagSeq/All_Temps.csv") %>%
  mutate(Incubator = factor(Incubator))
head(culture_exp)
str(culture_exp)

# work with date time format to subset date
culture_exp$DateTime<-strptime(culture_exp$DateTime, format="%m/%d/%y %H:%M")
culture_exp$datetime_ct <- as.POSIXct(culture_exp$DateTime, format="%Y-%m-%dT%H:%M:%S")

culture_exp$Date<-format(culture_exp$datetime_ct,"%D")
culture_exp$Date<-as.POSIXct(culture_exp$Day, format="%m/%d/%y")
culture_exp$Day = as.factor(culture_exp$Day)

#### Make plots ####

#cols = c("control" = "#a6611a", "cold" = "#74c476", "hot" = "#fd8d3c")
cols_host = c("control" = "#a6611a", "cold" = "#74c476", "hot" = "#fd8d3c")

host.plot = host_exp %>%
  ggplot(aes(x = day, y = temp_C, color = treatment))+
  theme_bw() +
  geom_point(aes(color = treatment), size=2, alpha = 0.9)+
  #geom_line() +
  scale_color_manual(name = "Treatment", values = cols_host)+
  xlab("Day") +
  scale_y_continuous(name = "Temperature (°C)", breaks = seq(2,32,4)) +
  scale_x_continuous(breaks = seq(1,15,1))+
  geom_hline(yintercept = 32, linetype = "solid", color = "#fd8d3c", size = 1) +
  geom_hline(yintercept = 6, linetype = "solid", color = "#74c476", size = 1) +
  ggtitle("A. Holobiont Experiment") +
  theme(legend.position = c(.1, .2), legend.background = element_rect(color = "black"))
host.plot
ggsave(host.plot, filename = "/Users/hannahaichelman/Documents/BU/Host_Buffering/MPCC_2018/Sym_analyses/plots/Host_Temp_Plot.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)

# now culture plot
cols_sym = c("control" = "#a6611a", "cool" = "#74c476", "heat" = "#fd8d3c")

culture.plot = culture_exp %>%
  ggplot(aes(x = Date, y = TempC, color = Incubator))+
  theme_bw() +
  geom_point(aes(color = Incubator), size=2, alpha = 0.9, shape = 21)+
  scale_color_manual(name = "Treatment", values = cols_sym)+
  xlab("Day") +
  scale_y_continuous(name = "Temperature (°C)", breaks = seq(2,32,4)) +
  scale_x_discrete(labels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15"))+
  geom_hline(yintercept = 32, linetype = "dashed", color = "#fd8d3c", size = 1) +
  geom_hline(yintercept = 6, linetype = "dashed", color = "#74c476", size = 1) +
  ggtitle("B. Culture Experiment") +
  theme(legend.position = c(.1, .2), legend.background = element_rect(color = "black"))
culture.plot
ggsave(culture.plot, filename = "/Users/hannahaichelman/Documents/BU/Host_Buffering/MPCC_2018/Sym_analyses/plots/Culture_Temp_Plot.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)


# combine host and symbiont plots
figs.combined = ggarrange(host.plot, culture.plot, ncol = 2, nrow = 1)
ggsave(figs.combined, filename = "/Users/hannahaichelman/Documents/BU/Host_Buffering/MPCC_2018/Sym_analyses/plots/Combined_Temp_Plot.pdf", width=8, height=4, units=c("in"), useDingbats=FALSE)

#### Summarize temperatures ####
str(host_exp)
host_temps = summarySE(data = host_exp, measurevar = "temp_C", groupvars = c("treatment","day"))
host_temps_fullavg = summarySE(data = host_exp, measurevar = "temp_C", groupvars = c("treatment"))

str(culture_exp)
culture_temps = summarySE(data = culture_exp, measurevar = "TempC", groupvars = c("Incubator","Date"))
culture_temps_fullavg = summarySE(data = culture_exp, measurevar = "TempC", groupvars = c("Incubator"))

