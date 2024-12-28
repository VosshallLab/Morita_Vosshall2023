#plot orco heat-seek analysis for Figure 2

library(tidyverse)
library(multcompView) #load for letter statistics
library(rstatix) #load for statistical testing
library(patchwork) #load for organizing panels
library(ggbeeswarm)
library(lemon)

#set path to this folder
path <- "set/path/to/tmorita_2025_data/"
source(paste0(path, "scripts/figure_group_colors.R"))

f <- paste0(path,"figure2/orco_dwell_quant.tsv")
dat <- read_tsv(file=f)

##########
#functions ----------------------------------------------------
##########
#plot data
plot_dwell_experiment <- function(dat,color_param,ylim){
  groups <- color_param$breaks %in% dat$group
  dat$group <- factor(dat$group, 
                         levels = color_param$breaks)
  
  dat %>% ggplot(aes(x = group, 
                     y = y,
                     fill=group, 
                     color=group)) +
    geom_boxplot(alpha=0.4,
                 outlier.shape = NA,
                 coef = 0) +
    geom_beeswarm(shape=1,cex=2) +
    scale_color_manual(values=color_param$col[groups]) +
    scale_fill_manual(values=color_param$col[groups]) +
    scale_x_discrete(breaks=color_param$breaks,
                     labels=color_param$labels) +
    theme_classic() +
    labs(y = "Landing frequency") +
    theme(legend.position = "none",
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank()) +
    theme(axis.line = element_line(color="black"),
          axis.text = element_text(size=8,color="black"),
          axis.ticks = element_line(color="black"),
          axis.title = element_text(size=9),
          legend.text = element_text(size=9)) +
    coord_capped_cart(ylim=ylim,left='both')
}  
#do stats
stats_anova <- function(dat){
  #ANOVA, one-way
  anova.rr <- aov(y ~ group, data=dat)
  summary(anova.rr)
  #Tukey posthoc test
  tukey.rr <- TukeyHSD(anova.rr)
  #Compact letter display
  cld.rr <-multcompLetters4(anova.rr, tukey.rr)
  #summary table
  exp_summary <- group_by(dat, group) %>%
    summarise(mean=mean(y), sd=sd(y), n=n()) %>%
    arrange(desc(mean)) %>%
    merge(data.frame(letters=cld.rr$group$Letters), 
          by.x = "group", 
          by.y = "row.names")
  stat_out <- list(anova.rr=anova.rr,
                   anova_summary = summary(anova.rr),
                   tukey.rr=tukey.rr,
                   exp_summary=exp_summary)
  return(stat_out)
}

#######################
#plot landing frequency ---------------------------------------
#######################
plot_orco_n_landing <- dat %>% 
  mutate(y = landing_n/n_mos/t_sampling) %>%
  plot_dwell_experiment(figcolors$orco,c(0,0.1))
#landing frequency stats
stats_orco_n_landing <- dat %>% 
  mutate(y = landing_n/n_mos/t_sampling) %>%
  stats_anova()

#######################
#plot takeoff frequency ---------------------------------------
#######################
plot_orco_n_takeoff <- dat %>% 
  mutate(y = takeoff_n/n_mos/t_sampling) %>%
  plot_dwell_experiment(figcolors$orco,c(0,0.1))
#takeoff frequency stats
stats_orco_n_takeoff <- dat %>% 
  mutate(y = takeoff_n/n_mos/t_sampling) %>%
  stats_anova()

################
#plot dwell time ----------------------------------------------
################
plot_orco_dwelltime <-dat %>% 
  mutate(y = time) %>%
  plot_dwell_experiment(figcolors$orco,c(0,10))
#swell time stats
stats_orco_dwelltime <- dat %>% 
  mutate(y = time) %>%
  stats_anova()
