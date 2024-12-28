#plot neurotranscriptome female data for Figure 8

library(tidyverse)
library(tximport)
library(DESeq2)
library(biomaRt)
library(GenomicFeatures)
library(reshape2)
library(data.table)
library(lemon)
library(ggbeeswarm2)
library(multcompView) #load for letter statistics

#set path to this folder
path <- "set/path/to/tmorita_2025_data/"
source(paste0(path, "scripts/figure_group_colors.R"))

##########
#functions --------------------------------------------
##########
#plot scatter plot
plot_gene_expression_scatter <- function(ncts,sample_table,gene_info,ylimits){
  dat <- merge(sample_table, t(ncts[gene_info$id,,drop=FALSE]), by='row.names')
  dat_long <- dat %>% pivot_longer(cols = !c("Row.names","sex", "tissue", "behavior", "replicate"),
                                   names_to = "gene",
                                   values_to = "counts")
  dat_long$tissue <- factor(dat_long$tissue,levels=tissue_id)
  dat_long$gene <- factor(dat_long$gene,levels=gene_info$id)
  
  
  plt <- dat_long %>% 
    # filter(behavior %in% c('SF','O')) %>%
    na.omit() %>%
    ggplot(aes(x=tissue, y=counts, fill=tissue, color=tissue)) +
    geom_boxplot(alpha=0.4,
                 outlier.shape = NA,
                 coef = 0)+
    geom_beeswarm(shape=1,cex=2)+
    facet_wrap(~gene,
               nrow = 1,
               scale="free_x") +
    theme_classic() +
    scale_color_manual(values=c("black",
                                "gray30",
                                "gray60")) +
    scale_fill_manual(values=c("black",
                               "gray30",
                               "gray60")) +
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
#stats
stat_tissue_expression <- function(dat){
  #ANOVA, one-way
  anova.rr <- aov(counts ~ tissue, data=dat)
  summary(anova.rr)
  #Tukey posthoc test
  tukey.rr <- TukeyHSD(anova.rr)
  #Compact letter display
  cld.rr <-multcompLetters4(anova.rr, tukey.rr)
  #summary table
  exp_summary <- group_by(dat, tissue) %>%
    summarise(mean=mean(counts), sd=sd(counts), n=n()) %>%
    arrange(desc(mean)) %>%
    merge(data.frame(letters=cld.rr$tissue$Letters), 
          by.x = "tissue", 
          by.y = "row.names")
  stat_out <- list(anova.rr=anova.rr,
                   anova_summary = summary(anova.rr),
                   tukey.rr=tukey.rr,
                   exp_summary=exp_summary)
  return(stat_out)
}

##########
#load data --------------------------------------------
##########
#list of tissues to look at
tissue_id <- c("FL",
               "ML",
               "HL")

tissue_name <- c("Foreleg",
                 "Midleg",
                 "Hindleg")

dir <- paste0(path,"figure8/neurotranscriptome_female_quants/")


#list of behavioral status to look at
behavior_id <- c("SF",
               "O")

behavior_name <- c("Sugar fed",
                 "Gravid")


#generate sample table
sample_names <- list.dirs(path=dir, full.names = FALSE, recursive = FALSE)
sampleID <- sapply(strsplit(sample_names,"_quant"), `[`, 1)
sample_sex <- sapply(strsplit(sample_names,"_"), `[`, 1)
sample_tissue <- sapply(strsplit(sample_names,"_"), `[`, 2)
sample_behaviorStatus <- sapply(strsplit(sample_names,"_"), `[`, 3)
sample_replicateN <- sapply(strsplit(sample_names,"_"), `[`, 4)
sample_table <- data.frame(names = sampleID,
                           sex = sample_sex,
                           tissue = sample_tissue,
                           behavior = sample_behaviorStatus,
                           replicate = sample_replicateN)
sample_table <- sample_table %>% column_to_rownames('names')

#select tissue sample to be analyzed
sample_table <- 
  sample_table %>% 
  filter(tissue %in% tissue_id) %>%
  filter(behavior %in% behavior_id)

#select files to be analyzed
files <- file.path(dir, paste0(rownames(sample_table),"_quant"), "quant.sf")
names(files) <- paste0(rownames(sample_table))
files <- files[rownames(sample_table)]

f <- paste0(path,"/figure8/AaegLVP_VB58-Jove19.gtf")
txdb <- makeTxDbFromGFF(f,format="gtf")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

#import files and normalize across tissues to be analyzed
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
dds <- DESeqDataSetFromTximport(txi, sample_table, ~tissue)
dds <- estimateSizeFactors(dds)
ncts <- counts(dds, normalized=TRUE)

##########################
#genes to analyze and plot -------------------------------------------------
##########################
gene_info_ir <- 
  data.frame(name = c("Ir166",
                      "Ir7d",
                      "Ir129",
                      "Ir139",
                      "Ir138",
                      "Ir75i",
                      "Ir140"),
             id = c("Ir166",
                    "Ir7d",
                    "Ir129",
                    "Ir139",
                    "Ir138",
                    "Ir75i",
                    "Ir140"))
plt_ir <- plot_gene_expression_scatter(ncts,
                                       sample_table,
                                       gene_info_ir,
                                       ylimits=c(0,500))


gene_info_ir25a <- 
  data.frame(name = c("Ir25a"),
             id = c("Ir25a"))
plt_ir25a <- plot_gene_expression_scatter(ncts,
                                          sample_table,
                                          gene_info_ir25a,
                                          ylimits=c(0,5000))

####################################################
#pairwise differential expression and assign letters --------------------------
####################################################
dds <- DESeq(dds)
myResFLvsML <- results(dds, contrast = c("tissue", "FL", "ML"))
myResFLvsML$newPadj <- p.adjust(myResFLvsML$pvalue,method="fdr")
myResMLvsHL <- results(dds, contrast = c("tissue", "ML", "HL"))
myResMLvsHL$newPadj <- p.adjust(myResMLvsHL$pvalue,method="fdr")
myResFLvsHL <- results(dds, contrast = c("tissue", "FL", "HL"))
myResFLvsHL$newPadj <- p.adjust(myResFLvsHL$pvalue,method="fdr")

gene_info_ir_ir25a <- rbind(gene_info_ir,gene_info_ir25a)

pwc_matrix <- 
  cbind("FL-ML"=myResFLvsML[gene_info_ir_ir25a$id,'newPadj',drop=F],
        "ML-HL"=myResMLvsHL[gene_info_ir_ir25a$id,'newPadj',drop=F],
        "FL-HL"=myResFLvsHL[gene_info_ir_ir25a$id,'newPadj',drop=F])
colnames(pwc_matrix) <- sapply(strsplit(colnames(pwc_matrix),"[.]"), `[`, 1)

pwc_matrix_logical <- pwc_matrix %>% as.matrix() < 0.01
pwc_matrix_logical_letters <- apply(pwc_matrix_logical, 1, multcompLetters)
