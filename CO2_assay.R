#plot CO2 experiments
library(tidyverse)
library(ggbeeswarm2)
library(lemon)
library(multcompView)
library(plotrix) #for stderr
library(patchwork)

#set path to this folder
path <- "set/path/to/tmorita_2023_data/"
source(paste0(path, "/figure_group_colors.R"))

#########
#Function ---------------------------------------------------------------
#########
#plot CO2 activation assay
plot_CO2_line <- function(df,color_param,ylimits){
  df$group <- factor(df$group, 
                     levels = color_param$breaks)
  
  df <- df %>% filter(time=='0' |
                        time=='60' | 
                        time=='120' |
                        time=='180' |
                        time=='240' |
                        time=='300')
  summary <- df %>% 
    group_by(group,time) %>% 
    summarise(activation_mean = mean(pFlying),
              actvation_sd = sd(pFlying),
              activation_se = std.error(pFlying))
  
  summary %>% ggplot(aes(time,activation_mean,fill=group,color=group)) +
    stat_summary(geom = "line",fun = mean,na.rm = TRUE) +
    geom_pointrange(aes(ymin = activation_mean - activation_se, ymax = activation_mean + activation_se)) +
    theme_classic() +
    labs(y="% mosquitoes flying",
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
#plot pre post quantification
plot_CO2_quantification <- function(df,color_param,ylimits){
  df$group <- factor(df$group, 
                     levels = color_param$breaks)
  df %>% 
    mutate_at('time', as.factor) %>%
    filter(time=='0' | time=='60') %>% 
    mutate(trialn = paste0(session,'_',trialName),
           condition = case_when(time == 0 ~ "Pre", time == 60 ~ "Post")) %>%
    mutate(condition = factor(condition, levels = c("Pre","Post"))) %>%
    ggplot(aes(x=condition,y=pFlying,color=group,fill=group)) +
    geom_beeswarm(color=alpha('black',0.1),shape=16)+
    geom_line(aes(group = trialn),color=alpha('black',0.1),na.rm = TRUE) + 
    theme_classic()+
    stat_summary(fun = median, geom = "line", aes(group = 1),na.rm = TRUE,cex=1) +
    stat_summary(fun = median, geom = "point", aes(group = 1),na.rm = TRUE, cex=1) +
    theme_classic()+
    labs(y='% mosquitoes flying')+
    scale_color_manual(values=color_param$col)+
    scale_fill_manual(values=color_param$fill) +
    facet_wrap(~group, scale="free_x",nrow = 1) +
    theme(panel.spacing.x = unit(.5, "lines"),
          legend.justification ="top",
          legend.title = element_blank(),
          legend.position = "none",
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank()) +
    theme(axis.line = element_line(color="black"),
          axis.text = element_text(size=8,color="black"),
          axis.ticks = element_line(color="black"),
          axis.title = element_text(size=9),
          legend.text = element_text(size=9)) +
    coord_capped_cart(ylim = ylimits, left='both')
}
#CO2 activation assay stats
stats_CO2 <- function(df,t){
  d <- df %>% filter(time == t)
  #ANOVA, one-wayd
  anova.rr <- aov(pFlying ~ group, data=d)
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

##################
#Orco CO2 response ------------------------------------------------------
##################
f <- paste0(path,"figure3/CO2_orco.tsv")
dat <- read_tsv(file = f)
orco_CO2_plot_line <- plot_CO2_line(dat,figcolors$orco,c(0,100))
orco_CO2_quantification <- plot_CO2_quantification(dat,figcolors$orco,c(0,100))
stats_CO2(dat,60)

###################
#Ir25a CO2 response -----------------------------------------------------
###################
f <- paste0(path,"figure1/CO2_ir25a.tsv")
dat <- read_tsv(file = f)
ir25a_CO2_plot_line <- plot_CO2_line(dat,figcolors$ir25a,c(0,100))
ir25a_CO2_quantification <- plot_CO2_quantification(dat,figcolors$ir25a,c(0,100))
stats_CO2(dat,300)

##############
#orco scitraks ----------------------------------------------------------
##############
f <- paste0(path,"figure2/orco_scitraks.tsv")
dat <- read_tsv(file = f)
dat$group <- factor(dat$group,
                                  levels=c('orl_wt',
                                           'orco_16_het',
                                           'orco_516_ko',
                                           'orco_25_ko'))
dat$stim <- factor(dat$stim,
                   levels=c('pre',
                            'post'))

figcolors$col <- c("black",
                   "gray30",
                   "firebrick3",
                   "firebrick3")

figcolors$fill <- c("black",
                    "gray30",
                    "firebrick3",
                    "firebrick3")

figcolors$breaks = c("orl_wt",
                     "orco_16_het",
                     "orco_516_ko",
                     "orco_25_ko")

figcolors$labels = c("Wild-type",
                     expression(italic("orco" ^ "16/+")),
                     expression(italic("orco" ^ "5/16+")),
                     expression(italic(paste("orco"^"2/5"))))

orco_CO2_scitraks <- dat %>%
  ggplot(aes(x = stim, y = distance*100, fill=group, color = group)) +
  geom_line(aes(group = trialn),color='gray80',na.rm = TRUE) + 
  geom_point(color='gray80',shape=16,cex=2) +
  stat_summary(fun = mean, geom = "line", aes(group = 1),na.rm = TRUE,cex=1) +
  stat_summary(fun = mean, geom = "point", aes(group = 1),na.rm = TRUE, cex=2.5) +
  scale_fill_manual(values = figcolors$fill) +
  scale_color_manual(values = figcolors$fill) +
  theme_classic() +
  facet_wrap(~group,
             nrow = 1,
             scale="free_x") +
  ylab("Total distance traveled (cm)") +
  theme(panel.spacing.x = unit(1, "lines"),
        strip.background = element_blank(),
        legend.justification ="top",
        legend.title = element_blank(),
        legend.text.align = 0)+
  theme(axis.line.x = element_blank(),
        axis.text = element_text(size=8,color="black"),
        axis.ticks = element_line(color="black"),
        axis.title = element_text(size=9),
        legend.text = element_text(size=9)) +
  coord_capped_cart(ylim=c(0,800),left='both')