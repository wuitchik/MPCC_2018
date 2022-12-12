# Analyze Oculina arbuscula physiology data

library(tidyverse)
library(Rmisc)

# Read in Fv/Fm data

pam = read.csv('/Users/hannahaichelman/Documents/BU/Host_Buffering/MPCC_2018/Phys_data/Oculina/Oculina_PAM.csv')

str(pam)
pam$day = as.factor(pam$day)
pam$sym_state = as.factor(pam$sym_state)
pam$genet = as.factor(str_sub(pam$coral_id, start = 1, end = 1))

pam2 = pam %>%
  mutate(avgfvfm = rowMeans(select(., fvfm1, fvfm2, fvfm3))) %>%
  filter(complete.cases(avgfvfm))
  

pam.sym = pam2 %>%
  subset(sym_state == 'Sym')

pam.apo = pam2 %>%
  subset(sym_state == 'Apo')

pam.sym.summary = summarySE(data = pam.sym, measurevar = "avgfvfm", groupvars = c("day", "treatment"))

pam.sym.plot = ggplot(pam.sym.summary, aes(x = day, y = avgfvfm, color = treatment))+
  geom_point(size = 3)+
  geom_errorbar()+
  
  
  
  