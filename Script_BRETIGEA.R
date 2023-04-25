##########################################################################
##ESTIMATE BRAIN CELL TYPE PROPORTIONS IN HIPPOCAMPUS DATA WITH BRETIGEA
##########################################################################

##Load packages
library(BRETIGEA)
library(limma)
library(illuminaRatv1.db)
library(reshape2)
library(ggplot2)
library(plyr)
library(dplyr)
library(grid)
library(Hmisc)

source("Script_centered_matrix.R")


##Read in expression matrix with probe names
#E.medians.relative.t.maxvar.m <- readRDS("E_medians_relative_t_maxvar_m.rds")

##Convert probe names to gene symbols through EntrezID to get up to date gene symbols
E.genes <- merge(conversion.df, E.medians.relative.t.maxvar.m, by = "ProbeID")

##Get cell marker genes from BRETIGEA. Select top markers
cell_types <- unique(markers_df_brain$cell)
names(cell_types) <- cell_types

markers_touse <- lapply(cell_types, function(x){
  sel.df <- markers_df_brain[markers_df_brain$cell %in% x,]
  touse <- sel.df[1:150,]$markers
  touse2 <- capitalize(tolower(touse))
})


##Import relative expression data and select top 50 markers of interest. Did not use the converted gene ids, but the ones from the original mouse genes
E.genes2 <- lapply(cell_types, function(x){y <- E.genes[E.genes$GeneID_updated %in% markers_touse[[x]],]; return(y)})

##Order according to importance rank
E.genes2.o <- lapply(cell_types, function(x){E.genes2[[x]][order(match(E.genes2[[x]]$GeneID_updated, markers_touse[[x]])),]})

##Select top 50 genes
exp_mat_rel_markers.list <- lapply(E.genes2.o, function(x){y <- x[x$GeneID_updated %in% unique(x$GeneID_updated)[1:50],]; return(y)})

##Make into data frame
exp_mat_rel_markers.df <- ldply(exp_mat_rel_markers.list, data.frame)

##Add variables
exp_mat_rel_markers.df$Portion <- as.factor(unlist(sapply(strsplit(as.character(exp_mat_rel_markers.df$sample_type), "-"), `[`, 2, simplify=FALSE)))
exp_mat_rel_markers.df$Age <- as.factor(unlist(sapply(strsplit(as.character(exp_mat_rel_markers.df$sample_type), "-"), `[`, 3, simplify=FALSE)))

##Run ANOVA analysis
marker.anova <- function(gene.list){
  df.selected <- exp_mat_rel_markers.df[toupper(exp_mat_rel_markers.df$GeneID_updated) %in% toupper(gene.list),]
  df.lm <- lm(Expression ~ Portion + Age, df.selected)
  df.anova <- anova(df.lm)
  return(df.anova)
}

marker.summary <- function(gene.list){
  df.selected <- exp_mat_rel_markers.df[toupper(exp_mat_rel_markers.df$GeneID_updated) %in% toupper(gene.list),]
  df.lm <- lm(Expression ~ Portion + Age, df.selected)
  return(summary(df.lm))
}

anova.results <- lapply(markers_touse, marker.anova)
summary.results <- lapply(markers_touse, marker.summary)

##Move rownames to new column
anova.results2 <- lapply(anova.results, function(x){as.data.frame(x)})
anova.results3 <- lapply(anova.results2, function(x){x$variable <- rownames(x); return(x)})
anova.results.df <- ldply(anova.results3, data.frame)

##Perform Bonferroni correction according to 6 tests (number of celltypes), adjusted p-value should not exceed 1
anova.results.df$Bonferroni <- unlist(lapply(1:length(anova.results.df$Pr..F.), function(x){min(anova.results.df$Pr..F.[x]*6, 1)}))

























