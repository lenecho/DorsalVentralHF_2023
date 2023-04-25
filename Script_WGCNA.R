###########################
##Run WGCNA per Mitch module and on whole data set
###########################

library(WGCNA)

source("Script_Mitch.R")

##Add expression matrix of all differentially expressed genes to list of exp mat per module
E.max.var.t.sig <- E.max.var.t[,colnames(E.max.var.t) %in% sig_genes]
module.expmat.list[[length(module.expmat.list) + 1]] <- E.max.var.t.sig

##Check for genes and samples with too many missing values
gsg.list <- lapply(module.expmat.list, goodSamplesGenes)
lapply(gsg.list, function(x) x$allOK)

##Make clustering trees of samples for each module expression matrix
sampleTrees.list <- lapply(module.expmat.list, function(x) {y = hclust(dist(x), method = "average"); return(y)})

samples = rownames(module.expmat.list[[1]])

##Get traits data
traits <- pData(data.N3)
traits$Labels <- paste(traits$HF.part, traits$Time.point, sep = "_")
traits2 <- traits$Labels
names(traits2) <- traits$Sample.switched

##Make numeric representation of traits for each sample (HF part and age). Assign colors per trait
traitRows = match(samples, traits$Sample.switched);
traits.factors <- traits %>% mutate_if(is.character,as.factor)
traits.factors <- traits.factors[,-1]
traits.factors$HF.part <- factor(traits.factors$HF.part, levels = c("DORSAL", "INTER", "VENT"))
traits.factors$Time.point <- factor(traits.factors$Time.point, levels = c("P0", "P9", "P18", "P60"))
traits.factors$Labels <- NULL
traits.numeric <- traits.factors %>% mutate_if(is.factor,as.numeric)
traitColors = numbers2colors(traits.numeric, signed = FALSE)

# Clean up any ongoing parallel computing
unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
unregister_dopar()

# Choose a set of soft thresholding powers
powers = c(1:20)  

# choose power based on SFT criterion for each tree
sft.list <- lapply(module.expmat.list, pickSoftThreshold, powerVector = powers)

# SFT index as a function of different powers
pdf("Power_selection_WGCNA2_sig.pdf")
par(mfrow = c(1, 2))

lapply(sft.list, function(x){
  plot(x$fitIndices[, 1], -sign(x$fitIndices[, 3]) * x$fitIndices[, 2], 
       xlab = "Soft Threshold (power)", ylab = "SFT, signed R^2", type = "n", main = paste("Scale independence"))
  text(x$fitIndices[, 1], -sign(x$fitIndices[, 3]) * x$fitIndices[, 2], 
       labels = powers, col = "red")
  
})
dev.off()

selected_powers <- c(8,9,6,7,6,7,7,7,6,9,9,6,7)
selected_powers <- c(4,7,13,6,7,7,7,5,14,5,7,13,7,5,6)


##Run network construction and module detection
net.list <- lapply(1:length(module.expmat.list), function(x){
  y <- blockwiseModules(module.expmat.list[[x]], power = selected_powers[x],
                        TOMType = "signed", minModuleSize = 20, deepSplit = 4,
                        maxBlockSize = 8300, networkType = "signed",
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs = TRUE,
                        saveTOMFileBase = "TOM1",
                        verbose = 3); return(y)
})

#saveRDS(net.list, file = "net_list.rds")


##Get list of eigengenes per module
nSamples = nrow(module.expmat.list[[1]])

eigengene_func <- function(x){
  print(x)
  moduleLabels = net.list[[x]]$colors
  moduleColors = labels2colors(moduleLabels)
  MEs0 = moduleEigengenes(module.expmat.list[[x]], moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  moduleTraitCor = cor(MEs, traits.numeric, use = "p")
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
  MEs$HF.part <- sapply(strsplit(rownames(MEs), "\\."), "[[", 2)
  MEs$Age <- sapply(strsplit(rownames(MEs), "\\."), "[[", 3)
  MEs$Age <- factor(MEs$Age, levels = c("P0", "P9", "P18", "P60"))
  MEs.m <- melt(MEs)
  return(MEs.m)
}

eigengene.list <- lapply(1:length(net.list), eigengene_func)

names(eigengene.list) <- c(paste(rep("Mitch", times = length(eigengene.list)-1), 1:(length(eigengene.list)-1), sep = ""), "Main")
eigengene.df <- ldply(eigengene.list, data.frame)

modulecolor.list <- lapply(net.list, function(x) x$colors)
modulegene.list <- modulecolor.list

##Make data frame listing genes in each group for each Mitch module
modulegene.df.list <- lapply(modulegene.list, function(x){y <- data.frame("Submodule" = x, "Genes" = names(x)); return(y)})
names(modulegene.df.list) <- paste("Mitch_module", 1:length(modulegene.df.list), sep = "_")
modulegene.df <- ldply(modulegene.df.list, data.frame)
colnames(modulegene.df) <- sub(".id", "Module", colnames(modulegene.df))
modulegene.df$Module <- sub(paste0("Mitch_module_", length(modulegene.df.list)), "Modules_allDEgenes", modulegene.df$Module)

##Add EntrezIDs to data frame
entrez_symbol$Genes <- toupper(entrez_symbol$GeneID_updated)
modulegene.df2 <- merge(modulegene.df, entrez_symbol)

##Add Illumina probes to data frame instead of EntrezID
tt.3 <- E.max.var.raw4
tt.3$Genes <- rownames(tt.3)
modulegene.df3 <- merge(modulegene.df, tt.3)
modulegene.df4 <- modulegene.df3

##Function to match color to module number
color_number_func <- function(Mitchmodule){
  module_numbers <- unname(modulegene.list[[Mitchmodule]])
  names(module_numbers) <- modulecolor.list[[Mitchmodule]]
  module_numbers_unique <- module_numbers[!duplicated(module_numbers)]
  return(module_numbers_unique)
}

modulecolor.list2 <- lapply(1:length(modulecolor.list), color_number_func)
modulecolor.list2 <- lapply(1:length(modulecolor.list), function(x) {names(modulecolor.list[[x]]) <- unname(modulegene.list[[x]])})

##Get pearson correlation between the eigengenes in each Mitch module and the overall eigengenes for all differentially expressed genes in the dataset
moduleLabels_main = net.list[[length(net.list)]]$colors
moduleColors_main = labels2colors(moduleLabels_main)
MEs0_main = moduleEigengenes(module.expmat.list[[length(net.list)]], moduleColors_main)$eigengenes
MEs_main = orderMEs(MEs0_main)

eigengene_func <- function(x){
  moduleLabels = net.list[[x]]$colors
  moduleColors = labels2colors(moduleLabels)
  MEs0 = moduleEigengenes(module.expmat.list[[x]], moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  moduleTraitCor = cor(MEs_main, MEs, use = "all.obs")
  colnames(moduleTraitCor) <- paste(colnames(moduleTraitCor), "cor", sep = "_")
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
  colnames(moduleTraitPvalue) <- paste(colnames(moduleTraitPvalue), "pval", sep = "_")
  module_corr_results <- cbind(moduleTraitCor, moduleTraitPvalue)
  return(module_corr_results)
}  

cor_result.list <- lapply(1:(length(net.list) - 1), eigengene_func)
cor_result.list2 <- lapply(1:length(cor_result.list), function(x){
  colnames(cor_result.list[[x]]) <- paste("Mitch", x, colnames(cor_result.list[[x]]), sep = "_"); return(cor_result.list)})


#####################
##List of genes in WGCNA modules for all DE genes
#####################

##Get list of genes for each expression signature
all_expsignatures_raw <- data.frame(net.list[[length(net.list)]]$colors)
all_expsignatures_raw$Gene <- rownames(all_expsignatures_raw)
all_expsignatures.list <- split(all_expsignatures_raw$Gene, all_expsignatures_raw[,1])

##Convert geneid to Entrezid
all_expsignatures_entrez.list <- lapply(all_expsignatures.list, function(x){y <- conversion.df[toupper(conversion.df$Gene) %in% x,]$EntrezID; return(y)})
