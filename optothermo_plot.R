#plot optothermocycler data
library(tidyverse)
library(ggplot2)
library(rstatix)
library(plotrix)
library(lemon)
library(multcompView)

#set path to this folder
path <- "set/path/to/tmorita_2023_data/"
source(paste0(path, "/figure_group_colors.R"))

#########
#function---------------------------------------------------------------
#########
plot_opto_behavior <- function(dat,color_param,ylim){
  match_group <- color_param$breaks %in% unique(dat$group)
  dat$group <- factor(dat$group, 
                         levels = color_param$breaks[match_group])
  
  dat %>% 
  group_by(group,window) %>% 
  summarise(pmos_mean = mean(pmos),
            pmos_sd = sd(pmos),
            pmos_se = std.error(pmos)) %>%
  ggplot(aes(x=window,
             y=pmos_mean,
             color=group,
             fill=group)) +
  geom_line(aes(group = group,
                linetype = group),
            position = position_dodge(0.05)) +
  geom_pointrange(aes(ymin = pmos_mean - pmos_se,
                      ymax = pmos_mean + pmos_se,
                      linetype = group),
                  shape=21,
                  fatten=5,
                  position = position_dodge(0.05)) +
  scale_fill_manual(values = color_param$fill[match_group]) +
  scale_color_manual(values = color_param$col[match_group]) +
  labs(y = '% mosquitoes responding') +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        legend.justification ="top",
        legend.title = element_blank(),
        legend.text.align = 0) +
  theme(axis.line = element_line(color="black"),
        axis.text = element_text(size=8,color="black"),
        axis.ticks = element_line(color="black"),
        axis.title = element_text(size=9),
        legend.text = element_text(size=9)) +
  coord_capped_cart(ylim=ylim,bottom='none', left='both')
}

#######################################
# plot antenna tip cut optothermocycler --------------------------------
#######################################
f <- paste0(path,"/figure4/optothermocycler_antennatipcut.tsv")
dat <- read_tsv(file = f)
plot_opto_tipcut <- plot_opto_behavior(dat,
                                        figcolors$antcut,
                                        c(0,50))
opto_tipcut_stats <- t.test(pmos ~ group, dat)

########################################
# plot antenna full cut optothermocycler--------------------------------
########################################
f <- paste0(path,"/figure4/optothermocycler_antennafullcut.tsv")
dat <- read_tsv(file = f)
plot_opto_fullcut <- plot_opto_behavior(dat,
                                        figcolors$antcut,
                                        c(0,50))
opto_fullcut_stats <- t.test(pmos ~ group, dat)

  