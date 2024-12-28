#plot feeding behavior for Figure 4

library(tidyverse)
library(ggbeeswarm2)
library(lemon)
library(multcompView)

#set path to this folder
path <- "set/path/to/tmorita_2025_data/"
source(paste0(path, "scripts/figure_group_colors.R"))

##########
#functions ---------------------------------------------------------------
##########
#plot feeding data
plot_feeding <- function(dat,color_param){
  
  groups <- color_param$breaks %in% dat$group
  dat$group <- factor(dat$group,levels=color_param$breaks)
  
  plt <- dat %>%
    ggplot(aes(x = group, y = fed/total*100, fill=group, color = group)) +
    geom_boxplot(alpha=0.4,
                 outlier.shape = NA,
                 coef = 0) +
    geom_beeswarm(shape=1,cex=2) +
    scale_color_manual(values=color_param$col[groups])+
    scale_fill_manual(values=color_param$fill[groups]) +
    theme_classic() +
    ylab("% blood-fed")+
    scale_x_discrete(breaks=color_param$breaks,
                     labels=color_param$labels) +
    theme(legend.position = "none",
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank()) +
    theme(axis.line = element_line(color="black"),
          axis.text = element_text(size=8,color="black"),
          axis.ticks = element_line(color="black"),
          axis.title = element_text(size=9),
          legend.text = element_text(size=9)) +
    coord_capped_cart(ylim=c(0,100),left='both')
  return(plt)
}

#do stats on feeding data
stat_feeding <- function(dat){
  #ANOVA, one-way
  dat <- dat %>% 
    mutate(percent_fed = fed/total*100)
  anova.rr <- aov(percent_fed ~ group, data=dat)
  summary(anova.rr)
  #Tukey posthoc test
  tukey.rr <- TukeyHSD(anova.rr)
  #Compact letter display
  cld.rr <-multcompLetters4(anova.rr, tukey.rr)
  #summary table
  exp_summary <- group_by(dat, group) %>%
    summarise(mean=mean(percent_fed), sd=sd(percent_fed), n=n()) %>%
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

##################################
#antenna tipcut glytube experiment ---------------------------------------
##################################
f <- paste0(path,"figure4/glytube_antennatipcut.tsv")
dat <- read_tsv(file = f)
antcut_glytubefeeding_plot <- plot_feeding(dat,figcolors$antcut)
antcut_glytube_stat <- stat_feeding(dat)

######################################
#antenna tipcut arm feeding experiment -----------------------------------
######################################
f <- paste0(path,"figure4/armfeed_antennatipcut.tsv")
dat <- read_tsv(file = f)
antcut_armfeeding_plot <- plot_feeding(dat,figcolors$antcut)
antcut_armfeeding_stat <- stat_feeding(dat)

#############################
#tarsi cut glytube experiment --------------------------------------------
#############################
f <- paste0(path,"figure4/glytube_tarsicut.tsv")
dat <- read_tsv(file = f)
tarsicut_glytubefeeding_plot <- plot_feeding(dat,figcolors$tarsicut)
tarsicut_glytube_stat <- stat_feeding(dat)

#################################
#tarsi cut arm feeding experiment -----------------------------------------
#################################
f <- paste0(path,"figure4/armfeed_tarsicut.tsv")
dat <- read_tsv(file = f)
tarsicut_armfeeding_plot <- plot_feeding(dat,figcolors$tarsicut)
tarsicut_armfeeding_stat <- stat_feeding(dat)



