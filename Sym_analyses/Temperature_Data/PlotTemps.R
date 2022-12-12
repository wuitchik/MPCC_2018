library(tidyverse)

setwd("/Users/hannahaichelman/Documents/BU/Host_Buffering/")

host_exp = read.csv("/Users/hannahaichelman/Documents/BU/Host_Buffering/MPCC_TagSeq/Oculina_temperature.csv") %>%
  select(-X) %>%
  mutate(treatment = factor(treatment))
head(host_exp)
str(host_exp)

sym_exp = read.csv("/Users/hannahaichelman/Documents/BU/Host_Buffering/Sym_TagSeq/All_Temps.csv") %>%
  mutate(Incubator = factor(Incubator))
head(sym_exp)
str(sym_exp)


cols = c("control" = "#a6611a", "cold" = "#74c476", "hot" = "#fd8d3c")
cols_host = c("control" = "#a6611a", "cold" = "#74c476", "hot" = "#fd8d3c")

host.plot.blues = host_exp %>%
  ggplot(aes(x = day, y = temp_C, color = treatment))+
  geom_point(aes(color = treatment), size=2, alpha = 0.9)+
  #geom_line() +
  scale_color_manual(name = "Treatment", values = cols_host)+
  xlab("Day") +
  scale_y_continuous(name = "Temperature (°C)", breaks = seq(2,32,4)) +
  scale_x_continuous(breaks = seq(1,16,1))+
  theme_bw() 
host.plot.blues
ggsave(host.plot, filename = "/Users/hannahaichelman/Documents/BU/Host_Buffering/Host_Temp_Plot.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)

cols_sym = c("control" = "#a6611a", "cool" = "#74c476", "heat" = "#fd8d3c")
sym.plot = sym_exp %>%
  ggplot(aes(x = Day, y = TempC, color = Incubator))+
  geom_point(aes(color = Incubator), size=2, alpha = 0.9, shape = 21)+
  scale_color_manual(name = "Treatment", values = cols_sym)+
  xlab("Day") +
  scale_y_continuous(name = "Temperature (°C)", breaks = seq(2,32,4)) +
  scale_x_continuous(breaks = seq(1,16,1))+
  theme_bw() 
sym.plot
ggsave(sym.plot, filename = "/Users/hannahaichelman/Documents/BU/Host_Buffering/Culture_Temp_Plot.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)


