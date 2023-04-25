#############################################
## MITCH ANALYSIS RAT REACTOME
#############################################

library(limma)
library(ggplot2)
library(plyr)
library(reshape2)
library(stringr)
library(mitch)
library(dplyr)
library(ReactomeContentService4R)
library(Hmisc)
library(illuminaRatv1.db)
library(ComplexHeatmap)
library(circlize)
library(factoextra)

source("Script_centered_matrix.R")

###############################
##Get Reactome gene list in gmt format for rat. 
###############################

rno_pathways <- getSchemaClass(class = "Pathway", species = "RNO", all = TRUE)
rno_reactome <- gmt_import("http://download.baderlab.org/EM_Genesets/February_05_2021/Rat/symbol/Pathways/Rat_Reactome_February_05_2021_symbol.gmt")

##Make all gene names upper case
rno_reactome2 <- lapply(rno_reactome, toupper)

##Remove unnecessary text from names of pathways
names(rno_reactome2) <- sapply(names(rno_reactome2), function(x) strsplit(x, "\\%")[[1]][1])

##########################################################
##Load DE lists and convert into right format for mitch
##########################################################

##Read in fit2
#fit2 <- readRDS("fit2.rds")

##Make list of toptables for each contrast
conversion.df$Gene <- conversion.df$GeneID_updated
DE.list <- lapply(contrasts, function(x){y <- topTable(fit2, coef=x, sort.by='none', n=Inf); return(y)})
DE.list2 <- lapply(DE.list, merge, y = conversion.df, by = "ProbeID", all = FALSE)

##Make rownames gene symbol. Remove probe and Entrez ID
DE.list3 <- lapply(DE.list2, function(x){y <- ddply(x, .(Gene), function(z) return(z[which(abs(z$logFC)==max(abs(z$logFC))),]))})
DE.list4 <- lapply(DE.list3, function(x){x <- x[!is.na(x$Gene),]; rownames(x) <- toupper(x$Gene); x})
DE.list5 <- lapply(DE.list4, function(x) {x$ProbeID <- NULL; x})
DE.list6 <- lapply(DE.list5, function(x) {x$TargetID <- NULL; x})
DE.list7 <- lapply(DE.list6, function(x) {x$Gene <- NULL; x$PROBE_ID.x <- NULL; x$PROBE_ID.y <- NULL; x$GeneID_updated <- NULL; x$GeneID <- NULL; x})

###################################
##Mitch analysis for all contrasts
###################################

mitch.list <- mitch_import(DE.list7, "limma")
mitch_allres_new <- mitch_calc(mitch.list, rno_reactome2, priority = "effect", resrows=10)
#saveRDS(mitch_allres_new, file = "mitch_allres_new.rds")
#mitch_allres_new <- readRDS("mitch_allres_new.rds")

##Get toplevel terms for pathways (terms that are high in the Reactome hierarchy)
rno_pathways <- getSchemaClass(class = "Pathway", species = "RNO", all = TRUE)
rno_pathways$Pathway <- toupper(rno_pathways$displayName)
pathway_names <- names(mitch_allres_new$input_genesets)
names(pathway_names) <- pathway_names
toplevel_func <- function(pathway){
  print(pathway)
  pathway_id <- rno_pathways[rno_pathways$Pathway %in% pathway,]$stId
  toplevel_res <- getPathways(pathway_id, top.level = TRUE)
  toplevel_term <- toplevel_res$name
  print(toplevel_term)
  return(toplevel_term)
}
pathway_names2 <- pathway_names[pathway_names %in% rno_pathways$Pathway]
pathway_names3 <- pathway_names2[!pathway_names2 %in% rno_pathways[duplicated(rno_pathways$Pathway),]$Pathway]

toplevels <- unlist(lapply(pathway_names3, toplevel_func))
toplevels.df <- data.frame("set" = names(toplevels), "toplevel" = toplevels, row.names = NULL)
toplevels.df.dups <- toplevels.df
toplevels.df.dups$set[!toplevels.df.dups$set %in% pathway_names3] <- sub(".$", "", toplevels.df.dups$set[!toplevels.df.dups$set %in% pathway_names3])

##Simplify classification
toplevels.df <- toplevels.df.dups %>% 
  mutate(Main = str_replace_all(toplevel, c("Cell Cycle" = "Cell Cycle/Apoptosis",
                                            "Blood coagulation" = "Other", 
                                            "Drug Absorption, Distribution, Metabolism and Excretion \\(ADME\\)" = "Other", 
                                            "Drug ADME" = "Other", 
                                            "Drug Pharmacokinetics \\(PK\\)" = "Other", 
                                            "Metabolism and Transport of Drugs \\(ADME\\)" = "Other", 
                                            "Metabolism of proteins" = "Metabolism", 
                                            "Reproduction" = "Other", 
                                            "Metabolism of RNA" = "Metabolism", 
                                            "Muscle contraction" = "Other",
                                            "Digestion and absorption" = "Other",
                                            "Gene expression$" = "Metabolism",
                                            "Gene expression \\(Transcription\\)" = "Metabolism",
                                            "Hemostasis" = "Other",
                                            "Signaling in Immune system" = "Immune System",
                                            "Immune System signaling" = "Immune System",
                                            "Signal Transduction" = "Signaling",
                                            "Signaling Pathways" = "Signaling",
                                            "Chromatin organization" = "Cell Cycle/Apoptosis",
                                            "DNA Repair" = "Cell Cycle/Apoptosis",
                                            "DNA Replication" = "Cell Cycle/Apoptosis",
                                            "Vesicle-mediated transport" = "Transport",
                                            "Protein localization" = "Transport",
                                            "Signaling$" = "Signaling",
                                            "Transport of small molecules" = "Transport",
                                            "Circadian Clock" = "Other",
                                            "Hemostasis" = "Other",
                                            "Programmed Cell Death" = "Cell Cycle/Apoptosis",
                                            "Cell-Cell communication" = "Signaling",
                                            "Cellular responses to external stimuli" = "Signaling",
                                            "Cellular responses to stimuli" = "Signaling",
                                            "Sensory Perception" = "Neuronal System")))

##Merge results with toplevel term
mitch_allres_new2 <- merge(mitch_allres_new$enrichment_result, unique(toplevels.df[,c("set", "Main")]))

####################################
##Plot heatmap
####################################

##Select only columns with effect size and reactome terms after selecting only significant terms. Exclude interaction contrasts. 
mitch_allres_new2_sig <- mitch_allres_new2[mitch_allres_new2$p.adjustMANOVA < 0.05,]
mitch_allres_new3 <- mitch_allres_new2_sig[,c(1,4:15, 25:39)]

##Make Reactome terms rownames and remove column
rownames(mitch_allres_new3) <- mitch_allres_new3$set
mitch_allres_new3$set <- NULL

##Determine optimal k
fviz_nbclust(mitch_allres_new3, kmeans, method = "gap_stat", k.max = 100)

##Plot all terms and all contrasts in one big heatmap - tried different values for row_km
colorder <- c("s.P0vsP9", "s.VP0vsP9", "s.IP0vsP9", "s.DP0vsP9", "s.P9vsP18", "s.VP9vsP18", "s.IP9vsP18", "s.DP9vsP18", 
              "s.P18vsP60", "s.VP18vsP60", "s.IP18vsP60", "s.DP18vsP60", "s.VD", "s.VD0", "s.VD9", "s.VD18","s.VD60", 
              "s.VI", "s.VI0", "s.VI9", "s.VI18", "s.VI60", "s.ID", "s.ID0", "s.ID9",  
              "s.ID18", "s.ID60")

legendorder <- c("Developmental Biology", "Extracellular matrix organization", 
                 "Organelle biogenesis and maintenance", "Cell Cycle/Apoptosis",
                 "Neuronal System", "Signaling", "Autophagy", "Metabolism", "Transport", "Other")
legendorder <- c("Developmental Biology", "Extracellular matrix organization", 
                 "Organelle biogenesis and maintenance", "Cell Cycle/Apoptosis",
                 "Neuronal System", "Signaling", "Immune System", "Metabolism", "Transport", "Other")
legendcols2 <- c("#01665e", "#35978f", "#80cdc1","#c7eae5","#b2182b","#d6604d","#f4a582", 
                 "#bf812d", "#dfc27d", "grey")
legendcols4 <- setNames(legendcols2, legendorder)

##Plot heatmap and cluster groups
pdf(file = "Reactome_heatmap3_sig15sel.pdf", height = 25, width = 13)
set.seed(41000)
row_ha = rowAnnotation(Main = mitch_allres_new2[mitch_allres_new2$set %in% rownames(mitch_allres_new3),]$Main,
                       col = list(Main = legendcols4), 
                       annotation_legend_param = list(at = names(legendcols4), labels = names(legendcols4))) 
ht = Heatmap(mitch_allres_new3, 
        clustering_distance_rows = "euclidean", 
        cluster_columns = FALSE,
        column_order = colorder,
        name = "Effect size", #title of legend
        column_title = "Contrasts", row_title = "Reactome terms",
        row_names_gp = gpar(fontsize = 3), # Text size for row names
        right_annotation = row_ha,
        row_km = 15, row_km_repeats = 600, clustering_method_rows = "complete"
)
#saveRDS(ht, file = "ht_mitch_heatmap3_sig.rds")
#ht <- readRDS("ht_mitch_heatmap3_sig.rds")
ht
dev.off()


#ht <- readRDS("ht_mitch_heatmap3_sig.rds")

###########################
##Get terms for each cluster
###########################
set.seed(41000)
htd = draw(ht)
c_terms <- lapply(1:length(row_order(htd)), function(x){y <- t(t(row.names(mitch_allres_new3[row_order(htd)[[x]],]))); return(y)})
names(c_terms) <- 1:length(row_order(htd))
c_terms.df <- ldply(c_terms, data.frame)
colnames(c_terms.df) <- c("Cluster", "Term")

#####################
##Get gene IDs for Reactome pathways of interest
#####################

##Get Reactome codes for pathways of interest
Reactomepathways <- getSchemaClass(class = "Pathway", species = "human", all = TRUE)

##Get Reactome codes for pathways in clusters
cluster_codes.df <- Reactomepathways[toupper(Reactomepathways$displayName) %in% c_terms.df$Term,c("displayName", "stId")]
cluster_codes.df$displayName <- toupper(cluster_codes.df$displayName)

##Add Reactome code information to dataframe with which pathways are in which clusters
colnames(cluster_codes.df) <- c("Term", "Code")
c_terms.df2 <- merge(c_terms.df, cluster_codes.df)

##Get Reactome codes for relevant pathways
cluster_codes <- cluster_codes.df$Code
names(cluster_codes) <- cluster_codes

##Get gene IDs for each pathway
gene_ids <- lapply(cluster_codes, event2Ids)
gene_ids2 <- lapply(gene_ids, function(x) {y <- x[["geneSymbol"]]; return(y)})

##Get all genes that are differentially expressed in one contrast or another
DE.list7.sig <- lapply(DE.list7, function(x){y <- x[x$adj.P.Val < 0.05,]; return(y)})
sig_genes <- unique(unlist(lapply(DE.list7.sig, function(x){y <- rownames(x); return(y)})))

########################
##Get significant genes for each Mitch term
########################

##Make list of all significant genes in each reactome code
termsig_genes <- lapply(gene_ids2, function(x) {y <- x[x %in% sig_genes]; return(y)})

##Make a list of all significant genes in each cluster (merging reactome codes for each cluster)
clustergenes_func <- function(cluster){
  codes <- c_terms.df2[c_terms.df2$Cluster %in% cluster,]$Code
  termsig_genesincluster <- termsig_genes[codes]
  genes <- unique(unlist(termsig_genesincluster))
  return(genes)
}

clustersig_genes.list <- lapply(1:length(c_terms), clustergenes_func)

#saveRDS(clustersig_genes.list, file = "clustersig_genes_list2_sig.rds")
#clustersig_genes.list <- readRDS("clustersig_genes_list2_sig.rds")

##################################
##Subset normalized expression matrix with the maximum variance probes 
##################################

##Expression matrix with only probes of maximum variance. Convert 
E.max.var.raw <- data.frame(E[rownames(E) %in% maxvar.probes,])
E.max.var.raw$ProbeID <- rownames(E.max.var.raw)
E.max.var.raw2 <- merge(E.max.var.raw, conversion.df, by = "ProbeID")
E.max.var.raw3 <- E.max.var.raw2[!is.na(E.max.var.raw2$GeneID_updated),]
E.max.var.raw4 <- E.max.var.raw3[,c(2:59)]
rownames(E.max.var.raw4) <- toupper(E.max.var.raw3$GeneID_updated)
E.max.var.t <- t(E.max.var.raw4)


##Get expression matrix per module
module.expmat.list <- lapply(clustersig_genes.list, function(x){y <- E.max.var.t[,colnames(E.max.var.t) %in% x]; return(y)})











