#plot heat-seek behavior for Figures 1,2,4,9
library(tidyverse)
library(rstatix) #load for statistical testing
library(ggbeeswarm2) #for scatter plots
library(lemon) #load for axis style
library(multcompView) #load for letter statistics
library(plotrix) #for stderr

#set path to this folder
path <- "set/path/to/tmorita_2025_data/"
source(paste0(path, "scripts/figure_group_colors.R"))

temp <- c(26,28.5,31,33.5,36,38.5,40,45,50,55,60)

##########
#functions ---------------------------------------------------------------
##########
#plot heatseek quantification data as line plot
plot_heatseeking_quant_line <- function(dat,color_param,ylim){
  dat$group <- factor(dat$group, 
                         levels = color_param$breaks)
  summary <- dat %>% 
    group_by(group,pelt_temp_step) %>% 
    summarise(heatseek_mean = mean(heatseek),
              heatseek_sd = sd(heatseek),
              heatseek_se = std.error(heatseek))
  
  plot_output <- ggplot(summary, 
                        aes(x = pelt_temp_step,
                            y = heatseek_mean,
                            group = group,
                            color = group,
                            fill=group)) +
    geom_line() +
    theme_classic() +
    geom_pointrange(aes(ymin = heatseek_mean - heatseek_se,
                        ymax = heatseek_mean + heatseek_se)) +
    labs(x = expression(paste("Peltier Temperature (", degree, "C)")), y = "% on peltier") +
    scale_fill_manual(values = color_param$fill,
                      labels = color_param$labels) +
    scale_color_manual(values = color_param$col,
                       labels = color_param$labels) +
    theme(legend.justification ="top",
        legend.title = element_blank(),
        legend.text.align = 0,
        strip.text.x = element_blank()) +
    theme(axis.line = element_line(color="black"),
          axis.text = element_text(size=8,color="black"),
          axis.ticks = element_line(color="black"),
          axis.title = element_text(size=9),
          legend.text = element_text(size=9)) +
    coord_capped_cart(ylim=ylim, left='both')
  return(plot_output)
}
#plot heatseek quantification data as scatter plot
plot_heatseek_quant_scatter <- function(dat,color_param,ylimit){
  groups <- color_param$breaks %in% dat$group
  dat$group <- factor(dat$group,levels=color_param$breaks)
  
  plt <- dat %>%
    ggplot(aes(x = as.character(pelt_temp_step), y = heatseek, fill=group, color = group)) +
    geom_boxplot(position=position_dodge(width = 1.1),
                 alpha=0.4,
                 outlier.shape = NA,
                 coef = 0)+
    geom_beeswarm(dodge.width=1.1,shape=1,cex=2) +
    scale_fill_manual(values = color_param$fill,
                      labels = color_param$labels) +
    scale_color_manual(values = color_param$col,
                       labels = color_param$labels) +
    theme_classic() +
    facet_wrap(~pelt_temp_step,
               nrow = 1,
               scale="free_x") +
    labs(x=expression(paste("Peltier (", degree,"C)")),y = "% on peltier") +
    scale_x_discrete(breaks=color_param$breaks,
                     labels=color_param$labels) +
    theme(panel.spacing.x = unit(1, "lines"),
          strip.background = element_blank(),
          legend.justification ="top",
          legend.title = element_blank(),
          legend.text.align = 0)+
    theme(axis.line = element_line(color="black"),
          axis.text = element_text(size=8,color="black"),
          axis.ticks = element_line(color="black"),
          axis.title = element_text(size=9),
          legend.text = element_text(size=9)) +
    coord_capped_cart(ylim=ylimit,left='both')
  return(plt)
}

#stats, one-way ANOVA
stats_heatseek <- function(dat,temp){
  dat <- dat %>% filter(pelt_temp_step == temp)
  #ANOVA, one-way
  anova.rr <- aov(heatseek ~ group, data=dat)
  summary(anova.rr)
  #Tukey posthoc test
  tukey.rr <- TukeyHSD(anova.rr)
  #Compact letter display
  cld.rr <-multcompLetters4(anova.rr, tukey.rr)
  #summary table
  exp_summary <- group_by(dat, group) %>%
    summarise(mean=mean(heatseek), sd=sd(heatseek), n=n()) %>%
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
stats_heatseek_all_temp <- function(dat,temp){
  out <-list()
  for (i in 1:length(temp)){
    s <- stats_heatseek(dat,temp[i])
    out[[length(out)+1]] <- s
    }
  return(out)
}

##################################
#plot orco heatseek quantification ----------------------------------------
##################################
f <- paste0(path,"/figure2/orco_quant_all.tsv")
dat <- read_tsv(file = f)
#plot heat-seek line plot for Orco experiment
plot_orco_heatseek_quant_line <- plot_heatseeking_quant_line(dat,figcolors$orco,c(0,45))
#do stats
stats_orco_heatseek_quant <- stats_heatseek_all_temp(dat,temp)

###################################
#plot Ir25a heatseek quantification ---------------------------------------
###################################
f <- paste0(path,"/figure1/ir25a_quant_all.tsv")
dat <- read_tsv(file = f)
#plot heat-seek line plot for Ir25a experiment
plot_ir25a_heatseek_quant_line <- plot_heatseeking_quant_line(dat,figcolors$ir25a,c(0,45))
#do stats
stats_ir25a_heatseek_quant <- stats_heatseek_all_temp(dat,temp)

##################################
#plot Ir76b heatseek quantification----------------------------------------
##################################
f <- paste0(path,"/figure1/ir76b_quant_all.tsv")
dat <- read_tsv(file = f)
#plot heat-seek line plot for Ir76b experiment
plot_ir76b_heatseek_quant_line <- plot_heatseeking_quant_line(dat,figcolors$ir76b,c(0,45))
#do stats
stats_ir76b_heatseek_quant <- stats_heatseek_all_temp(dat,temp)

##################################
#plot Ir8a heatseek quantification ----------------------------------------
##################################
f <- paste0(path,"/figure1/ir8a_quant_all.tsv")
dat <- read_tsv(file = f)
#plot heat-seek line plot for Ir8a experiment
plot_ir8a_heatseek_quant_line <- plot_heatseeking_quant_line(dat,figcolors$ir8a,c(0,45))
#do stats
stats_ir8a_heatseek_quant <- stats_heatseek_all_temp(dat,temp)

##################################
#tip cut heatseek quantification ------------------------------------------
##################################
f <- paste0(path,"/figure4/heatseek_antennatipcut.tsv")
dat <- read_tsv(file = f)
#plot heat-seek line plot for antenna tipcut experiment
plot_tipcut_heatseek_quant_line <- plot_heatseeking_quant_line(dat,figcolors$antcut,c(0,20))
#do stats
stats_tipcut_heatseek_quant <- stats_heatseek_all_temp(dat,temp)

##################################
#tarsi cut heatseek quantification ------------------------------------------
##################################
f <- paste0(path,"/figure4/heatseek_tarsicut.tsv")
dat <- read_tsv(file = f)
tarsicut_heatseek_plot <- plot_heatseek_quant_scatter(dat,figcolors$tarsicut,c(0,30))
tarsicut_heatseek_stat <- stats_heatseek_all_temp(dat,c(26,36,55))

##################################
#Ir140 Orco double mutant heatseek quantification ---------------------------
##################################
f <- paste0(path,"/figure9/ir140_orco_doubleknockout.tsv")
dat <- read_tsv(file = f)
ir140_orco_heatseek_plot <- plot_heatseek_quant_scatter(dat,figcolors$ir140_orco,c(0,35))
ir140_orco_heatseek_stat <- stats_heatseek_all_temp(dat,c(26,36,55))

##################################
#Ir140 single mutant heatseek quantification -------------------------------
##################################
f <- paste0(path,"/figure9/ir140_knockout.tsv")
dat <- read_tsv(file = f)
ir140_heatseek_plot <- plot_heatseek_quant_scatter(dat,figcolors$ir140,c(0,35))
ir140_heatseek_stat <- stats_heatseek_all_temp(dat,c(26,36,55))

####################
# heat+DEET heatseek ------------------------------------------------------
####################
f <- paste0(path,"/figure2/orco_heatwithdeet.tsv")
dat <- read_tsv(file = f)

dat$group <- factor(dat$group,
                                          levels=c('orl_wt',
                                                   'orco_16_het',
                                                   'orco_516_ko',
                                                   'orco_25_ko'))
dat$treatment <- factor(dat$treatment,
                                              levels=c('solvent',
                                                       'deet'))

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

figcolors$labels = c("Wildtype",
                     expression(italic("orco" ^ "16/+")),
                     expression(italic("orco" ^ "5/16+")),
                     expression(italic(paste("orco"^"2/5"))))

plot_heatwithDEET_heatseek_quant <- dat %>%
  ggplot(aes(x = treatment, y = heatseek, fill=group, color = group)) +
  geom_boxplot(alpha=0.4,
               outlier.shape = NA,
               coef = 0) +
  geom_beeswarm(shape=1,
                cex=2) +
  scale_fill_manual(values = figcolors$fill) +
  scale_color_manual(values = figcolors$fill) +
  theme_classic() +
  facet_wrap(~group,
             nrow = 1,
             scale="free_x") +
  ylab("% on peltier") +
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
  coord_capped_cart(ylim=c(0,100),left='both')


