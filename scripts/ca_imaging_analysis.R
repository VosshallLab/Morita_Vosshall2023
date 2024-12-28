#A script to plot figures for Figure 6

library(tidyverse)
library(viridis)
library(scales)
library(lemon)
library(rstatix)
library(multcompView)

###############
###functions###
###############
#plot heatmap
plot_heatmap <- function(d){
  d %>%
    ggplot(aes(x=time,y=neuronID,fill = deltaRR, color=deltaRR)) + 
    geom_tile()+
    theme_void() +
    scale_fill_viridis(option = "A",limits=c(-0.1, 1.5),oob=squish)+
    scale_color_viridis(option = "A",limits=c(-0.1, 1.5),oob=squish)+
    labs(y="ROI",) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(angle = 90),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          line = element_blank(),
          legend.title = element_blank(),
          panel.background = element_blank()) +
    theme(axis.title = element_text(size=9),
          legend.text = element_text(size=8),
          legend.title = element_text(size=8))
}

#plot summary
plot_summary <- function(d){
  d %>%
    ggplot(aes(x=factor(tissue,levels = c("antenna", "palp", "prob", "legs")),y=peak)) +
    geom_violin(color = NA,
                fill = "black",
                alpha=0.3,
                scale = "width") +
    geom_boxplot(color = "black",
                 fill = NA,
                 width=0.2,
                 alpha=0.5,
                 outlier.shape = NA)+
    theme_classic() +
    theme(axis.title.x=element_blank(),
          #axis.text.x=element_blank(),
          axis.line.x.bottom=element_blank(),
          axis.ticks.x=element_blank(),
          #axis.line.y = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none") +
    theme(axis.line = element_line(color="black"),
          axis.text = element_text(size=8,color="black"),
          axis.ticks = element_line(color="black"),
          axis.title = element_text(size=9),
          legend.text = element_text(size=9)) +
    theme(legend.position = "none") +
    coord_capped_cart(ylim=c(-0.2,2.5),left='both')
}

#stats
do_stats <- function(mat,dep_var,ind_var){
  res.kruskal <- mat %>% kruskal_test(as.formula(str_c(dep_var,"~",ind_var)))
  pwc <- mat %>% dunn_test(as.formula(str_c(dep_var,"~",ind_var)), p.adjust.method = "bonferroni",detailed = TRUE)
  letter_matrix <- c()
  c <- as.vector(mapply('paste',pwc$group1,pwc$group2,sep="-"))
  test <- pwc$p.adj
  names(test) <- c
  letters <- multcompLetters(test)
  letter_matrix <- rbind(letter_matrix,letters$Letters)
  return(list(kuskal = res.kruskal, pwc = pwc, letter = letter_matrix))
}
###############

#import data
path <- "set/path/to/tmorita_2025_data/"
dat <- read_tsv(paste0(path,"figure6/ca_imaging_full_data.tsv"))
dat_summary <- read_tsv(paste0(path,"figure6/ca_imaging_summary_dat.tsv"))

#plot antenna heat response heatmap 
antenna_heat_heatmap <-
  dat %>% 
  filter(tissue == 'antenna' & stim == 'heat') %>%
  plot_heatmap()

#plot antenna kcl response heatmap
antenna_kcl_heatmap <-
  dat %>% 
  filter(tissue == 'antenna' & stim == 'kcl') %>%
  plot_heatmap()

#plot maxillary palp heat response heatmap 
maxillarypalp_heat_heatmap <-
  dat %>% 
  filter(tissue == 'palp' & stim == 'heat') %>%
  plot_heatmap()

#plot maxillary palp kcl response heatmap
maxillarypalp_kcl_heatmap <-
  dat %>% 
  filter(tissue == 'palp' & stim == 'kcl') %>%
  plot_heatmap()

#plot proboscis heat response heatmap 
proboscis_heat_heatmap <-
  dat %>% 
  filter(tissue == 'prob' & stim == 'heat') %>%
  plot_heatmap()

#plot proboscis kcl response heatmap
proboscis_kcl_heatmap <-
  dat %>% 
  filter(tissue == 'prob' & stim == 'kcl') %>%
  plot_heatmap()

#plot leg heat response heatmap 
leg_heat_heatmap <-
  dat %>% 
  filter(tissue == 'legs' & stim == 'heat') %>%
  plot_heatmap()

#plot leg kcl response heatmap
leg_kcl_heatmap <-
  dat %>% 
  filter(tissue == 'legs' & stim == 'kcl') %>%
  plot_heatmap()

#plot peak heat data
peak_heat_summary_plot <-
  dat_summary %>% 
  filter(stim == 'heat') %>%
  plot_summary()

#plot peak kcl data
peak_kcl_summary_plot <-
  dat_summary %>% 
  filter(stim == 'kcl') %>%
  plot_summary()

#summary statistic for heat response
heat_peak_stats <- 
  dat_summary %>%
  filter(stim == 'heat') %>%
  do_stats("peak","tissue")

#summary statistic for kcl response
kcl_peak_stats <- 
  dat_summary %>%
  filter(stim == 'kcl') %>%
  do_stats("peak","tissue")
