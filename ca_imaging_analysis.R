#calcium imaging figure
library(tidyverse)
library(matrixStats)
library(rstatix)
library(multcompView)
library(lemon)
library(data.table)

#set path to this folder
path <- "set/path/to/tmorita_2023_data/"
source(paste0(path, "/figure_group_colors.R"))

##########
#functions ----------------------------------------------------------------
##########
# stats
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
#auc function
auc_function <- function(y){
  n <- length(y)
  0.5*(y[1]+y[n]+2*sum(y[-c(1,n)]))
}
#plot average trace
plot_response <- function(full_dat,t){
  trace_plot <- 
    full_dat %>%
    filter(tissue==t) %>%
    ggplot(aes(x=time,y=deltaRR,col=genotype,fill=genotype)) +
    stat_summary(geom = "line",fun = mean,na.rm = TRUE) +
    stat_summary(geom="ribbon",fun.data = mean_se, alpha=0.4,color=NA) +
    theme_classic()+
    scale_colour_manual(breaks = breaks,
                        values = fill,
                        labels = labels,
                        name = NULL) + 
    scale_fill_manual(breaks = breaks,
                      values = alpha(fill,0.4),
                      labels = labels,
                      name = NULL) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
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
    coord_capped_cart(ylim=c(-0.1,.8),left='both')
}  
#plot temperature
plot_temp <- function(full_dat,t){
  trace_temperature <- 
    full_dat %>%
    filter(tissue==t) %>%
    unite(traialID,c(batch,sampleID,sampleType)) %>%
    group_by(traialID,time) %>%
    summarize(temperature = mean(temperature)) %>%
    group_by(time) %>%
    ggplot(aes(x=time,y=temperature)) +
    stat_summary(geom = "line",fun = mean,na.rm = TRUE) +
    stat_summary(geom="ribbon",fun.data = mean_se, alpha=0.4,color=NA) +
    theme_classic() +
    theme(#axis.title.x=element_blank(),
      #axis.text.x=element_blank(),
      axis.line.x.bottom=element_blank(),
      #axis.ticks.x=element_blank(),
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
    coord_capped_cart(ylim=c(20,35),left='both')
}
#plot quantification
plot_quantification <- function(summary_dat,x,y,t,ylimit){
  summary_dat %>%
    filter(tissue==t) %>%
    ggplot(aes_string(x=x,y=y,col=x,fill=x)) +
    geom_violin(color = NA,
                alpha=0.3) +
    geom_boxplot(width=0.2,
                 alpha=0.5,
                 outlier.shape = NA)+
    theme_classic()+
    scale_colour_manual(breaks = breaks,
                        values = fill,
                        labels = labels,
                        name = NULL) + 
    scale_fill_manual(breaks = breaks,
                      values = alpha(fill,0.4),
                      labels = labels,
                      name = NULL) +
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
    coord_capped_cart(ylim=ylimit,left='both')
}

#figure colors
breaks = c("wt",
           "ht",
           "ko")
fill <- c("black",
          "gray50",
          "firebrick3")
labels = c("Wildtype",
           expression(italic("Orco" ^ "16/+")),
           expression(italic("Orco" ^ "5/16")))


f <- paste0(path,"figure6/CaImaging_full_data.tsv")
full_dat <- read_tsv(file = f)

summary_dat <- dat %>% 
  group_by(genotype,batch,sampleID,tissue,sampleType,cell) %>% 
  summarise(peakresp = max(deltaRR[50:150]),
            auc = auc_function(deltaRR[50:150]))
summary_dat <- data.frame(summary_dat)
summary_dat$genotype <- factor(summary_dat$genotype,
                            levels = c("wt", "ht", "ko"))
#plots
fl_trace <- plot_response(full_dat,"fl")
ml_trace <- plot_response(full_dat,"ml")
fl_temp <- plot_temp(full_dat,"ml")
ml_temp <- plot_temp(full_dat,"ml")
fl_peak <- plot_quantification(summary_dat,"genotype","peakresp","fl",c(-.3,3))
ml_peak <- plot_quantification(summary_dat,"genotype","peakresp","ml",c(-.3,3))
fl_auc <- plot_quantification(summary_dat,"genotype","auc","fl",c(-30,250))
ml_auc <- plot_quantification(summary_dat,"genotype","auc","ml",c(-30,250))

#stats
summary_dat %>% filter(tissue == 'fl') %>% do_stats("peakresp","genotype")
summary_dat %>% filter(tissue == 'fl') %>% do_stats("auc","genotype")
summary_dat %>% filter(tissue == 'ml') %>% do_stats("peakresp","genotype")
summary_dat %>% filter(tissue == 'ml') %>% do_stats("auc","genotype")
