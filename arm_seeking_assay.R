#plot arm seeking experiments
library(tidyverse)
library(ggbeeswarm2)
library(lemon)
library(multcompView)
library(plotrix) #for stderr

#set path to this folder
path <- "set/path/to/tmorita_2023_data/"
source(paste0(path, "/figure_group_colors.R"))

##########
#Functions ----------------------------------------------------------------
##########
#plot attraction assay 
plot_arm_seeking_line <- function(d,color_param,ylimits){

  d$group <- factor(d$group, 
                    levels = color_param$breaks)

  d <- d %>% filter(time=='0' |
                      time=='60' | 
                      time=='120' |
                      time=='180' |
                      time=='240' |
                      time=='300')
  summary <- d %>% 
    group_by(group,time) %>% 
    summarise(armseek_mean = mean(mos_p),
              armseek_sd = sd(mos_p),
              armseek_se = std.error(mos_p))
  
  summary %>% ggplot(aes(time,armseek_mean,fill=group,color=group)) +
    stat_summary(geom = "line",fun = mean,na.rm = TRUE) +
    geom_pointrange(aes(ymin = armseek_mean - armseek_se, ymax = armseek_mean + armseek_se)) +
    theme_classic() +
    labs(y="% next to arm",
        x="time (min)") +
    scale_fill_manual(values = color_param$col,
                      labels = color_param$labels) +
    scale_color_manual(values = color_param$col,
                       labels = color_param$labels) +
    scale_x_continuous(breaks=seq(0,300,60),
                       labels=c("0", "1", "2","3","4","5")) +
    theme(axis.line = element_line(color="black"),
          axis.text = element_text(size=8,color="black"),
          axis.ticks = element_line(color="black"),
          axis.title = element_text(size=9),
          legend.text = element_text(size=9)) +
    coord_capped_cart(ylim = ylimits,left='both',bottom = 'both')
}
#stats for attraction assay
stats_armseek <- function(d,t){
  d <- d %>% filter(time == t)
  #ANOVA, one-way
  anova.rr <- aov(mos_p ~ group, data=d)
  summary(anova.rr)
  #Tukey posthoc test
  tukey.rr <- TukeyHSD(anova.rr)
  #Compact letter display
  cld.rr <-multcompLetters4(anova.rr, tukey.rr)
  #summary table
  exp_summary <- group_by(d, group) %>%
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

###########
#worn nylon ---------------------------------------------------------------
###########
f <- paste0(path,"figure3/nylon_seeking_worn.tsv")
dat <- read_tsv(file = f)
wornnylon_plot_line <- plot_arm_seeking_line(dat,figcolors$orco,c(0,50))
stats_armseek(dat,60)

#############
#unworn nylon -------------------------------------------------------------
#############
f <- paste0(path,"figure3/nylon_seeking_unworn.tsv")
dat <- read_tsv(file = f)
unwornnylon_plot_line <- plot_arm_seeking_line(dat,figcolors$orco,c(0,50))
stats_armseek(dat,60)

#########
#arm orco ----------------------------------------------------------------
#########
f <- paste0(path,"figure3/arm_seeking_arm.tsv")
dat <- read_tsv(file = f)
arm_plot_line <- plot_arm_seeking_line(dat,figcolors$orco,c(0,50))
stats_armseek(dat,180)

##############
#arm+deet orco -----------------------------------------------------------
##############s
f <- paste0(path,"figure3/arm_seeking_deet.tsv")
dat <- read_tsv(file = f)
armDEET_plot_line <- plot_arm_seeking_line(dat,figcolors$orco,c(0,50))
stats_armseek(dat,0)



