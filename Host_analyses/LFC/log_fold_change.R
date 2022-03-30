### Log fold change analysis 

# libraries
library(tidyverse)
library(DESeq2)
library(patchwork)
library(ggpubr) 

# load astrangia data
load("../DESeq2/Astrangia/Astrangia_sym_state.RData")

white_astrangia_cold = as.data.frame(results(dds.sym, alpha = 0.05, contrast = c("Treatment_by_Sym.Status", "Cold_White", "Control_White")))
white_astrangia_heat = as.data.frame(results(dds.sym, alpha = 0.05, contrast = c("Treatment_by_Sym.Status", "Heat_White", "Control_White")))
brown_astrangia_cold = as.data.frame(results(dds.sym, alpha = 0.05, contrast = c("Treatment_by_Sym.Status", "Cold_Brown", "Control_Brown")))
brown_astrangia_heat = as.data.frame(results(dds.sym, alpha = 0.05, contrast = c("Treatment_by_Sym.Status", "Heat_Brown", "Control_Brown")))
                                    
# overwrite with oculina
load("../DESeq2/Oculina/Oculina_sym_state.RData")

white_oculina_cold = as.data.frame(results(dds.sym, alpha = 0.05, contrast = c("Treatment_by_Sym.Status", "Cold_White", "Control_White")))
white_oculina_heat = as.data.frame(results(dds.sym, alpha = 0.05, contrast = c("Treatment_by_Sym.Status", "Heat_White", "Control_White")))
brown_oculina_cold = as.data.frame(results(dds.sym, alpha = 0.05, contrast = c("Treatment_by_Sym.Status", "Cold_Brown", "Control_Brown")))
brown_oculina_heat = as.data.frame(results(dds.sym, alpha = 0.05, contrast = c("Treatment_by_Sym.Status", "Heat_Brown", "Control_Brown")))

# merge data sets
cold_white_data = white_astrangia_cold %>%
  merge(white_oculina_cold, by=0) %>%
  dplyr::rename(GeneID = "Row.names",
                Astrangia = "log2FoldChange.x",
                Oculina = "log2FoldChange.y")

cold_brown_data = brown_astrangia_cold %>%
  merge(brown_oculina_cold, by=0) %>%
  dplyr::rename(GeneID = "Row.names",
                Astrangia = "log2FoldChange.x",
                Oculina = "log2FoldChange.y")

heat_white_data = white_astrangia_heat %>%
  merge(white_oculina_heat, by=0) %>%
  dplyr::rename(GeneID = "Row.names",
                Astrangia = "log2FoldChange.x",
                Oculina = "log2FoldChange.y")

heat_brown_data = brown_astrangia_heat %>%
  merge(brown_oculina_heat, by=0) %>%
  dplyr::rename(GeneID = "Row.names",
                Astrangia = "log2FoldChange.x",
                Oculina = "log2FoldChange.y")


# make into a figure

cold_white_LFC = ggplot(cold_white_data, aes(Astrangia, Oculina)) +
  geom_point(color = "#0D47A1", shape = 1) +
  geom_smooth(method = "lm", color = "black") +
  xlim(-8, 8) +
  ylim(-8, 8) +
  stat_regline_equation(label.y = 5.5, aes(label = ..eq.label..)) +
  stat_regline_equation(aes(label = ..rr.label..)) +
  theme_classic()

cold_brown_LFC = ggplot(cold_brown_data, aes(Astrangia, Oculina)) +
  geom_point(color = "#0D47A1") +
  geom_smooth(method = "lm", color = "black") +
  xlim(-8, 8) +
  ylim(-8, 8) +
  stat_regline_equation(label.y = 5.5, aes(label = ..eq.label..)) +
  stat_regline_equation(aes(label = ..rr.label..)) +
  theme_classic()

heat_white_LFC = ggplot(heat_white_data, aes(Astrangia, Oculina)) +
  geom_point(color = "#C62828", shape = 1) +
  geom_smooth(method = "lm", color = "black") +
  stat_regline_equation(label.y = 5.5, aes(label = ..eq.label..)) +
  stat_regline_equation(aes(label = ..rr.label..)) +
  xlim(-8, 8) +
  ylim(-8, 8) +
  theme_classic()

heat_brown_LFC = ggplot(heat_brown_data, aes(Astrangia, Oculina)) +
  geom_point(color = "#C62828") +
  geom_smooth(method = "lm", color = "black") +
  stat_regline_equation(label.y = 5.5, aes(label = ..eq.label..)) +
  stat_regline_equation(aes(label = ..rr.label..)) +
  xlim(-8, 8) +
  ylim(-8, 8) +
  theme_classic()

# make full figure
(cold_white_LFC + cold_brown_LFC) / (heat_white_LFC + heat_brown_LFC) + plot_annotation(tag_levels = "A")

# test if significantly different


