library(cowplot)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(gridExtra)
library(latex2exp)

custom_plot_theme <-
  list(
    theme_set(theme_bw(base_size = 16) +
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      axis.ticks.length = unit(-0.25, 'cm'),
                      axis.text.x = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), 'cm'), size = 16),
                      axis.text.y = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), 'cm'), size = 16),
                      plot.title = element_text(hjust = 0.5, size = 16),
                      legend.title = element_blank(),
                      legend.position = 'top',
                      legend.text = element_text(size = 16))),
    scale_color_manual(values = c('#02679a', '#f89821', '#3ac6f3', '#4d4d4f', '#e75480')),
    scale_fill_manual(values = c('#02679a', '#f89821', '#3ac6f3', '#4d4d4f', '#e75480')))

plot_colors <- c('#02679a', '#f89821', '#3ac6f3', '#4d4d4f', '#e75480')
