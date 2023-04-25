#####################################################
##SCRIPT TO MAKE CENTERED MATRIX FOR HIPPOCAMPUS DATA
##Dorsal, intermediate, and ventral hippocampus gene expression
##Ages P0, P9, P18, P60
##Illumina microarray
#####################################################

#############################################
### Packages.
#############################################
library(limma)
library(ascii)
library(illuminaRatv1.db)
library(genefilter)
library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)


#####################################
##Load normalized expression matrix
#####################################

source("Script_limma.R")
#data.N3 <- readRDS("data.N3.rds")

E <- exprs(data.N3)

##Calculate median expression for each condition (time point and hippocampal part)
E.medians <- aggregate(t(E), list(Condition=sub("^..", "", colnames(E))), median)
rownames(E.medians) <- E.medians$Condition
E.medians <- E.medians[,-1]

##Calculate median value for each probe
colMedians <- apply(E.medians, 2, median)

##Center data around median value, rescale to 0. Transpose and eliminate probes without variance.
E.medians.relative <- sweep(E.medians, 2, colMedians)
E.medians.relative.t <- t(E.medians.relative)

E.relative.tfin.m <- melt(E.medians.relative.t)

##Samples ordered
samples <- c("-DORSAL-P0", "-DORSAL-P9", "-DORSAL-P18", "-DORSAL-P60", "-INTER-P0", "-INTER-P9", "-INTER-P18", 
             "-INTER-P60", "-VENT-P0", "-VENT-P9", "-VENT-P18", "-VENT-P60")
E.relative.tfin.m$Var2 <- factor(E.relative.tfin.m$Var2, ordered = TRUE, levels = samples)

  
##Choose probe with highest variance for genes with more than one probe
E.variance <- data.frame(apply(E, 1, var, na.rm = TRUE))
colnames(E.variance) <- "variance"
E.variance$ProbeID <- rownames(E.variance)
E.variance2 <- merge(E.variance, conversion.df)
max.var.df <- E.variance2 %>% group_by(GeneID_updated) %>% top_n(1, variance)
maxvar.probes <- max.var.df$ProbeID

##Subset centered matrix with selected probes. Go from 12047 to 9152 probes.
E.medians.relative.t.maxvar <- E.medians.relative.t[rownames(E.medians.relative.t) %in% maxvar.probes,]
E.medians.relative.t.maxvar.m <- melt(E.medians.relative.t.maxvar)
E.medians.relative.t.maxvar.m$Var2 <- factor(E.medians.relative.t.maxvar.m$Var2, ordered = TRUE, levels = samples)
colnames(E.medians.relative.t.maxvar.m) <- c("ProbeID", "sample_type", "Expression")

saveRDS(E.medians.relative.t.maxvar.m, file = "E_medians_relative_t_maxvar_m.rds")

