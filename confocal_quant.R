#Plot confocal imagind data
library(tidyverse)
library(rstatix) #load for statistical testing
library(ggbeeswarm2) #for scatter plots
library(lemon) #load for axis style
library(multcompView) #load for letter statistics
library(plotrix) #for stderr

#set path to this folder
path <- "path/to/this/folter/tmorita_2023_data/"

##########
#Funcitons -------------------------------------------------------------
##########
#plot scatter plot
plot_count <- function(dat,ylimits){
  plt <- dat %>% 
    ggplot(aes(x=legs, y=count, color=driver)) +
    geom_boxplot(alpha=0.4,
                 outlier.shape = NA,
                 coef = 0)+
    geom_beeswarm(shape=1,cex=2)+
    facet_wrap(~driver,
               nrow = 1,
               scale="free_x") +
    theme_classic() +
    scale_color_manual(breaks = c("brp","ir25a","ir76b","gr4"),
                       values=c("black","dodgerblue3","orchid3","green")) +
    theme(panel.spacing.x = unit(1, "lines"),
          strip.background = element_blank(),
          legend.position='none',
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          #strip.text.x = element_blank(),
          axis.title.x = element_blank())+
    theme(axis.ticks = element_line(color = "black"),
          axis.text=element_text(size=8,color="black"),
          axis.title=element_text(size=9,color="black"),
          legend.text=element_text(size=9)) +
    ylab('Normalized read counts') +
    coord_capped_cart(ylim = ylimits,left='both')
  return(plt)
}
#plot cranial caudal counts
plot_cranialcaudal_count <- function(dat,ylimits){
  plt <- dat %>% 
    ggplot(aes(x=side, y=count, color=driver)) +
    geom_boxplot(alpha=0.4,
                 outlier.shape = NA,
                 coef = 0)+
    geom_beeswarm(shape=1,cex=2)+
    facet_wrap(~driver,
               nrow = 1,
               scale="free_x") +
    theme_classic() +
    scale_color_manual(breaks = c("brp","ir25a","ir76b","gr4"),
                       values=c("black","dodgerblue3","orchid3","green")) +
    theme(panel.spacing.x = unit(1, "lines"),
          strip.background = element_blank(),
          legend.position='none',
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          #strip.text.x = element_blank(),
          axis.title.x = element_blank())+
    theme(axis.ticks = element_line(color = "black"),
          axis.text=element_text(size=8,color="black"),
          axis.title=element_text(size=9,color="black"),
          legend.text=element_text(size=9)) +
    ylab('Normalized read counts') +
    coord_capped_cart(ylim = ylimits,left='both')
  return(plt)
}
#stat
stat_confocal <- function(dat){
  #ANOVA, one-way
  anova.rr <- aov(count ~ legs, data=dat)
  summary(anova.rr)
  #Tukey posthoc test
  tukey.rr <- TukeyHSD(anova.rr)
  #Compact letter display
  cld.rr <-multcompLetters4(anova.rr, tukey.rr)
  #summary table
  exp_summary <- group_by(dat, legs) %>%
    summarise(mean=mean(count), sd=sd(count), n=n()) %>%
    arrange(desc(mean)) %>%
    merge(data.frame(letters=cld.rr$legs$Letters), 
          by.x = "legs", 
          by.y = "row.names")
  stat_out <- list(anova.rr=anova.rr,
                   anova_summary = summary(anova.rr),
                   tukey.rr=tukey.rr,
                   exp_summary=exp_summary)
  return(stat_out)
}

######################
#Plot cell body counts -------------------------------------------------
######################
#cell body count
f <-paste0(path,"figure5/tarsi_cell_count.tsv")
dat <- read_tsv(file = f)
dat$legs <- factor(dat$legs, levels=c("Foreleg","Midleg","Hindleg"))
dat$driver <- factor(dat$driver, levels=c("brp","ir25a","ir76b","gr4"))
cell_count_plot <- plot_count(dat, y=c(0,100))

#####################
#Plot sensilla counts --------------------------------------------------
#####################
#sensilla count
f <-paste0(path,"figure5/tarsi_sensilla_count.tsv")
dat <- read_tsv(file = f)
dat$legs <- factor(dat$legs, levels=c("Foreleg","Midleg","Hindleg"))
dat$driver <- factor(dat$driver, levels=c("brp","ir25a","ir76b","gr4"))
sensilla_count_plot <- plot_count(dat, y=c(0,25))

###########################
#Plot cranial-caudal counts --------------------------------------------
###########################
#cranial-caudal count
f <-paste0(path,"figure5/tarsi_cranialcaudal_count.tsv")
dat <- read_tsv(file = f)
dat$side <- factor(dat$side, levels=c("cranial","caudal"))
dat$driver <- factor(dat$driver, levels=c("brp","ir25a","ir76b","gr4"))
#plot forelegs
dat %>% filter(legs == 'foreleg') %>% plot_cranialcaudal_count(y=c(0,60))
#plot midlegs
dat %>% filter(legs == 'midleg') %>% plot_cranialcaudal_count(y=c(0,60))