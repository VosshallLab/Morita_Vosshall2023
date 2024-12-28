#DE analysis using most updated transcript annotation files for Figure 8

library(tidyverse)
library(tximport)
library(DESeq2)
library(GenomicFeatures)
library(reshape2)
library(data.table)
library(lemon)
library(ggbeeswarm2) #for scatter plots
library(ggnewscale)

#set path to this folder
path <- "set/path/to/tmorita_2025_data/"
source(paste0(path, "scripts/figure_group_colors.R"))

##########
#functions ---------------------------------------------------------
##########
#perform differential expression using DESeq2
#and adjust p value
diffExp <- function(txi,sampleTable,graph,cutoff){
  #perform differential expression using DEseq2
  dds <- DESeqDataSetFromTximport(txi, sampleTable, ~genotype)
  #filter out genes < 10
  dds <- estimateSizeFactors(dds)
  idx <- rowSums( counts(dds, normalized=TRUE) >= 10 ) >= 10
  dds <- dds[idx,]
  #
  dds <- DESeq(dds)
  normCounts <- counts(dds, normalized = TRUE)
  myRes <- results(dds, contrast = c("genotype", "orco516", "ORL"))
  myRes <- myRes[order(myRes$pvalue), ]
  
  if(graph == TRUE){
    #DE summary
    plotDispEsts(dds)
    summary(myRes)
    summary(myRes_lfc)
    plotMA(myRes)
    plotMA(myRes_lfc)
  }
  
  #adjust pvalue
  myResAsDF <- as.data.frame(myRes)
  myResAsDF$newPadj <- p.adjust(myResAsDF$pvalue,method="fdr")
  myResAsDF <- myResAsDF[!is.na(myResAsDF$padj), ]
  myResAsDF <- myResAsDF[order(myResAsDF$pvalue), ]
  group <- rep("NS",dim(myResAsDF)[1])
  orco_enriched <- which(myResAsDF$newPadj < cutoff & myResAsDF$log2FoldChange > 0)
  group[orco_enriched] <- "orco"
  wt_enriched <- which(myResAsDF$newPadj < cutoff & myResAsDF$log2FoldChange < 0)
  group[wt_enriched] <- "wt"
  myResAsDF <- cbind(myResAsDF,group)
  return(list(myResAsDF=myResAsDF,
              normCounts=normCounts))
}

#plot fold change
plot_foldchange_ordered <- function(dat,ylimit){
  
  tmp <- dat %>% 
    filter(group != 'ND') %>%
    arrange(by=log2FoldChange)
  tmp <- tmp %>% mutate(xpos = 1:n())
  
  tmp |> 
    ggplot()+
    geom_hline(yintercept=0, linetype="dashed",color='gray60') +
    geom_point(aes(xpos, log2FoldChange,
                   size = log10(baseMean)),
               color = 'black',
               fill = 'white',
               filter(tmp, group == 'NS'),
               pch = 21,
               #size = 8,
               alpha = 1) +
    scale_size(limits=c(1,5),
               range = c(0, 5))+
    new_scale_colour() +
    new_scale_fill() +
    geom_point(aes(xpos, log2FoldChange,
                   fill = log10(newPadj),
                   size = log10(baseMean)),
               color = 'black',
               filter(tmp, group == 'wt'),
               pch = 21,
               #size = 8,
               alpha = 1,
               stroke = 1.5) +
    scale_fill_gradient(low="black", high="white",limits=c(-30,0)) +
    new_scale_colour() +
    new_scale_fill() +
    geom_point(aes(xpos, log2FoldChange,
                   fill = log10(newPadj),
                   size = log10(baseMean)),
               color = 'firebrick3',
               filter(tmp, group == 'orco'),
               pch = 21,
               alpha = 1,
               stroke = 1.5) +
    scale_fill_gradient(low="firebrick3", high="white",limits=c(-30,0)) +
    theme_classic() +
    theme(axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    theme(axis.ticks = element_line(color = "black"),
          axis.text=element_text(size=8,color="black"),
          axis.title=element_text(size=9,color="black"),
          legend.text=element_text(size=9)) +
    coord_capped_cart(ylim = ylimit,left='both')
}

plot_gene_expression <- function(txi,gene_list,gene_name,ylimit){
  chemoR <- txi[gene_list,]
  chemoR <- rownames_to_column(data.frame(chemoR))
  chemoR_table <- melt(data.table(chemoR),id.vars = 'rowname')
  tmp <- samples %>% column_to_rownames('sampleN')
  chemoR_table$genotype <- tmp[chemoR_table$variable,'genotype']
  chemoR_table$rowname <- factor(chemoR_table$rowname,
                                 levels = c(gene_list))
  
  plt <- chemoR_table %>% ggplot(aes(x=rowname, y=value, fill=factor(genotype,levels = c('ORL','orco516')))) +
    geom_boxplot(position=position_dodge(1)) +
    theme_classic() +
    scale_fill_manual(labels=c("ORL" = "wild type",
                               "orco516" = bquote(italic("orco")^italic("5/16"))),
                      values=c("grey",
                               "firebrick3")) +
    scale_x_discrete(breaks=gene_list,
                     labels=gene_name) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.title = element_blank(),
          legend.justification = "top",
          axis.text.x = element_text(angle = 45,vjust = 1, hjust=1,face = "italic"),
          axis.line.x = element_blank(),
          axis.title.x = element_blank()) +
    ylab('Transcripts per million') +
    coord_capped_cart(ylim = ylimit,left='both')
  return(plt)
}

plot_gene_expression_normcounts <- function(full_dat,ylimit){
  
  # chemoR <- full_dat[gene_list,] %>% 
  #   cbind(gene_name) %>%
  #   arrange(baseMean)
  
  chemoR <- full_dat%>%
    arrange(baseMean)
  
  gene_list <- chemoR$Row.names
  gene_name <- chemoR$name
  
  chemoR <- chemoR %>%
    dplyr::select(starts_with('sample'))
  
  chemoR <- rownames_to_column(data.frame(chemoR))
  
  chemoR_table <- melt(data.table(chemoR),id.vars = 'rowname')
  tmp <- samples %>% column_to_rownames('sampleN')
  chemoR_table$genotype <- tmp[chemoR_table$variable,'genotype']
  chemoR_table$rowname <- factor(chemoR_table$rowname,
                                 levels = c(gene_list))
  chemoR_table$genotype <- factor(chemoR_table$genotype,
                                  levels = c('ORL','orco516'))
  
  plt <- chemoR_table %>% ggplot(aes(x=rowname, y=value, fill=genotype,color=genotype)) +
    geom_boxplot(position=position_dodge(width = 1.1),
                 alpha=0.4,
                 outlier.shape = NA,
                 coef = 0)+
    geom_beeswarm(dodge.width=1.1,shape=1,cex=1)+
    facet_wrap(~rowname,
               nrow = 1,
               scale="free_x") +
    scale_fill_manual(labels=c("ORL" = "wild type",
                               "orco516" = bquote(italic("orco")^italic("5/16"))),
                      values=c("black",
                               "firebrick3")) +
    scale_color_manual(labels=c("ORL" = "wild type",
                                "orco516" = bquote(italic("orco")^italic("5/16"))),
                       values=c("black",
                                "firebrick3")) +
    scale_x_discrete(breaks=gene_list,
                     labels=gene_name) +
    theme_classic()+
    theme(panel.spacing.x = unit(1, "lines"),
          strip.background = element_blank(),
          legend.justification ="top",
          legend.title = element_blank(),
          legend.text.align = 0,
          strip.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45,hjust=0.95))+
    theme(axis.ticks = element_line(color = "black"),
          axis.text=element_text(size=8,color="black"),
          axis.title=element_text(size=9,color="black"),
          legend.text=element_text(size=9)) +
    ylab('Normalized read counts') +
    coord_capped_cart(ylim = ylimit,left='both')
  return(plt)
}


#############
#Load dataset -------------------------------------------------------
#############
#define sample IDs
samples<- data.frame(sampleN = paste0("sample",1:10),
                     sample = c("WT1","WT2","WT3","WT4","WT5","KO1","KO2","KO3","KO4","KO5"),
                     tissue = rep("tarsi",10),
                     genotype = c(rep("ORL",5),rep("orco516",5)))

#merge transcript IDs to gene IDs
f <- paste0(path,"/figure8/AaegLVP_VB58-Jove19.gtf")
txdb <- makeTxDbFromGFF(f,format="gtf")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

#load quantified mapped reads
dir <- paste0(path,"/figure8/orco_vs_wt_quants")
files <- file.path(dir, paste0(samples$sample,"_quant"), "quant.sf")
names(files) <- paste0("sample", 1:10)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

################
#Analyze dataset -------------------------------------------------------
################
#make metadata file for all samples
sampleTable <- samples %>% column_to_rownames(var='sampleN')
sampleTable$genotype <- factor(sampleTable$genotype, levels=c('ORL','orco516'))
DEresults <- diffExp(txi,sampleTable,graph = FALSE,cutoff=0.01)

full_dat <- merge(DEresults$normCounts,DEresults$myResAsDF,by = 0, all=TRUE) %>% 
  column_to_rownames('Row.names') %>% 
  mutate(group = ifelse(is.na(group), 'ND', group))

##########
#analyze by chemoreceptors ---------------------------------------------
##########
#load chemoreceptor gene list
f <- paste0(path,"/figure8/chemoreceptor_gene_list.tsv")
chemoreceptor_genelist <- read_tsv(f) %>% column_to_rownames('id')
chemoreceptor_dat <- full_dat[rownames(chemoreceptor_genelist),]
chemoreceptor_dat <- merge(chemoreceptor_dat, chemoreceptor_genelist, by=0)
rownames(chemoreceptor_dat) <- chemoreceptor_dat$Row.names

irs <- chemoreceptor_dat %>% filter(gene_category == 'ir')
grs <- chemoreceptor_dat %>% filter(gene_category == 'gr')
ors <- chemoreceptor_dat %>% filter(gene_category == 'or')
trps <- chemoreceptor_dat %>% filter(gene_category == 'trp')
ppks <- chemoreceptor_dat %>% filter(gene_category == 'ppk')

ir_fc_plot <- plot_foldchange_ordered(irs,c(-2,2))
or_fc_plot <- plot_foldchange_ordered(ors,c(-2,2))
gr_fc_plot <- plot_foldchange_ordered(grs,c(-2,2))
ppk_fc_plot <- plot_foldchange_ordered(ppks,c(-2,2))
trp_fc_plot <- plot_foldchange_ordered(trps,c(-2,2))

#load neuromodulator gene list
f <- paste0(path,"/figure8/neuromodulator_list.tsv")
neuromodulator_genelist <- read_tsv(f)
#neuropeptide receptors
npr_list <- neuromodulator_genelist %>% filter(type == 'neuropeptide receptor')
npr_dat <- full_dat[unique(npr_list$geneID),] %>% drop_na() %>% rownames_to_column("Row.names")
npr_dat2 <- npr_dat
npr_fc_plot <- plot_foldchange_ordered(npr_dat,c(-3,3))
npr_fc_plot <- npr_fc_plot + scale_y_continuous(breaks = seq(-3,3,by=1))

#neuropeptide
np_list <- neuromodulator_genelist %>% filter(type == 'neuropeptide')
np_dat <- full_dat[unique(np_list$geneID),] %>% drop_na() %>% rownames_to_column("Row.names")
np_dat2 <- np_dat
np_fc_plot <- plot_foldchange_ordered(np_dat,c(-3,3))
np_fc_plot <- np_fc_plot + scale_y_continuous(breaks = seq(-3,3,by=1))

#neurotransimitter biosynthesis genes
nt_list <- neuromodulator_genelist %>% 
  filter(!type %in% c('neuropeptide', 'neuropeptide receptor')) %>%
  filter(!grepl('receptor',type))
nt_dat <- full_dat[unique(nt_list$geneID),] %>% drop_na() %>% rownames_to_column("Row.names")
nt_dat2 <- nt_dat
nt_fc_plot <- plot_foldchange_ordered(nt_dat,c(-3,3))
nt_fc_plot <- nt_fc_plot + scale_y_continuous(breaks = seq(-3,3,by=1))

#neurotransimitter receptor genes
ntr_list <- neuromodulator_genelist %>% 
  filter(!type %in% c('neuropeptide', 'neuropeptide receptor')) %>%
  filter(grepl('receptor',type))
ntr_dat <- full_dat[unique(ntr_list$geneID),] %>% drop_na() %>% rownames_to_column("Row.names")
ntr_dat2 <- ntr_dat
ntr_fc_plot <- plot_foldchange_ordered(ntr_dat,c(-3,3))
ntr_fc_plot <- ntr_fc_plot + scale_y_continuous(breaks = seq(-3,3,by=1))


##########
#orco plot ---------------------------------------------------
##########
orco_expression <-
  full_dat['Orco',] %>% 
  dplyr::select(starts_with('sample')) %>% 
  t() %>%
  merge(sampleTable,by=0)

orco_scatter_plot <-
  orco_expression %>%
  ggplot(aes(x=genotype,y=Orco,color=genotype,fill=genotype))+
  geom_boxplot(alpha=0.4,
               outlier.shape = NA,
               coef = 0)+
  geom_beeswarm(shape=1,cex=2) +
  theme_classic() +
  scale_color_manual(values=c("black",
                              "firebrick3")) +
  scale_fill_manual(values=c("black",
                             "firebrick3")) +
  theme(axis.title.x = element_blank(),
        legend.position = 'none')+
  theme(axis.ticks = element_line(color = "black"),
        axis.text=element_text(size=8,color="black"),
        axis.title=element_text(size=9,color="black"),
        legend.text=element_text(size=9)) +
  ylab('Normalized read counts') +
  coord_capped_cart(ylim = c(0,250),left='both')


#####################################
#differential expression gene classes -------------------------
#####################################
f <- paste0(path,"/figure8/AaegLVP_VB58-Jove19_geneTypes.tsv")
gene_types <- 
  read_tsv(f) %>% 
  column_to_rownames('geneID')

rownames(full_dat) %in% rownames(gene_types) %>% length()

#get list of differentially expresed genes
orcogenes <- full_dat %>% filter(group == 'orco') %>% rownames() #614
wtgenes <- full_dat %>% filter(group == 'wt') %>% rownames() #688

#protein coding
protgeneslist <-  gene_types %>% filter(geneClass == "protein_coding_gene")
nprotgenes <- sum(rownames(full_dat) %in% rownames(protgeneslist))
nprotorcogenes <- sum(rownames(protgeneslist) %in% orcogenes) #515
nprotwtgenes <- sum(rownames(protgeneslist) %in% wtgenes) #613

#non coding genes
ncgeneslist <-  gene_types %>% filter(geneClass == "ncRNA_gene")
nncgenes <- sum(rownames(full_dat) %in% rownames(ncgeneslist))
nncorcogenes <- sum(rownames(ncgeneslist) %in% orcogenes) #90
nncotwtgenes <- sum(rownames(ncgeneslist) %in% wtgenes) #68

#pseudogenes
pseudogeneslist <-  gene_types %>% filter(geneClass == "pseudogene")
npseudogenes <- sum(rownames(full_dat) %in% rownames(pseudogeneslist))
npseudoorcogenes <- sum(rownames(pseudogeneslist) %in% orcogenes) #9
npseudowtgenes <- sum(rownames(pseudogeneslist) %in% wtgenes) #7



