##################################################
##DIFFERENTIAL EXPRESSION ANALYSIS USING LIMMA
##################################################

##Load R packages
library(lumi)
library(limma)
library(illuminaRatv1.db)
library(genefilter)
library(geneplotter)
library(plyr)
library(ArrayExpress)

require(lumi)


##Import raw data from Array Express
test <- getAE("E-MTAB-2617", type = "full", path="./data")

dataFn <- test$rawArchive

## read sample probe profile - 22522 probes
data <- lumiR(dataFn, convertNuID=FALSE, inputAnnotation=FALSE, QC=FALSE, dec=",",
              columnNameGrepPattern=list(exprs="AVG_Signal", se.exprs="BEAD_STDERR", detection="Detection", beadNum="Avg_NBEADS"))

#################
##Preprocessing
#################
probes <- unlist(mget(featureNames(data), illuminaRatv1ARRAYADDRESS, ifnotfound=NA))
keep <- !is.na(probes)
data <- data[keep,]
featureNames(data) <- probes[keep]
fData(data)$ProbeID <- probes[keep]

##Read genomestudio samplesheet
desc <- read.table(test$sdrf, header=TRUE, sep="\t")
sids <- desc$Assay.Name
rownames(desc) <- sids
desc$Sample_Name <- rownames(desc)
desc$location.switch <- c("dorsal hippocampus"="D", "intermediate hippocampus"="I", "ventral hippocampus"="V")[desc$Factor.Value.organism.part.]
desc$location.switch[desc$Sample_Name == "G1-DORSAL-P9"] <- "V"
desc$location.switch[desc$Sample_Name == "G1-VENT-P9"] <- "D"
desc$time.switch <- desc$Factor.Value.age.

data <- data[,sids] # ensure data is in same order as phenodata


##Remove outliers: G4-INTER-P18 had low intensity, G3-INTER-P60 high intensity
outNames <- c("G4-INTER-P18", "G3-INTER-P60")
if (!is.null(outNames)){
  data.0 <- data
  desc.0 <- desc
  cat('Removing outliers\n')
  out <- which(sampleNames(data) %in% outNames)
  data2 <- data[,-out]
  desc <- desc[-out,]
}


##Filter genes based on intensity of signal and annotation. Results in 12047 probes.
geneFilter <- TRUE
if (geneFilter == TRUE){
  presentLim <- 2
  present <- detectionCall(data, Th=0.01, "probe")
  keepProbes <- present >= presentLim
  data.filt <- data2[keepProbes,]
  ## genelevel filter
  annotation(data.filt) <- "illuminaRatv1"
}

#Perform variance transformation
data.T <- lumiT(data.filt, "log2")

#Perform between array normalization
data.N2 <- lumiN(data.T, "quantile")

##Incorporate switches in data object (two samples were switched). 
data.N3 <- data.N2

desc$Animal.group <- sapply(as.character(desc$Sample_Name), function(x) strsplit(x, "-")[[1]][1])
desc$Sample.switched <- paste(desc$Animal.group, paste("-", 
                                                       paste(desc$location.switch, 
                                                             paste("-", paste("P", desc$time.switch, sep = ""), sep = ""), sep = ""), sep = ""), sep = "")
desc$Sample.switched <- sub("V", "VENT", sub("I", "INTER", sub("D", "DORSAL", desc$Sample.switched)))
desc.short <- desc[, c("Sample_Name", "Sample.switched")]
desc.short2 <- cbind(desc.short, colnames(exprs(data.N3)))
sampleNames(data.N3) <- desc.short2$Sample.switched
colnames(exprs(data.N3)) <- desc.short2$Sample.switched
desc.short2$HF.part <- sapply(desc.short2$Sample.switched, function(x) strsplit(x, "-")[[1]][2])
desc.short2$Time.point <- sapply(desc.short2$Sample.switched, function(x) strsplit(x, "-")[[1]][3])
desc.switched <- desc.short2[, c("Sample.switched", "HF.part", "Time.point")]
rownames(desc.switched) <- desc.switched$Sample.switched
phenoData(data.N3) <- AnnotatedDataFrame(desc.switched)

nn.switch <- sampleNames(data.N3)
group <- as.factor(sapply(nn.switch, function(x) strsplit(x, "-")[[1]][1]))
pt <- as.factor(sapply(nn.switch, function(x) strsplit(x, "-")[[1]][3]))
animal <- factor(paste(group, pt, sep="_"))

bm <- paste(pData(data.N3)[,"HF.part"], pData(data.N3)[,"Time.point"], sep="")
cl <- as.factor(bm)

################################
##Design matrix and limma model
################################

des <- model.matrix(~-1+cl)
colnames(des) <- levels(cl)
dc <- duplicateCorrelation(data.N3, design=des, block=animal)
cc <- dc$cor
fit <- lmFit(data.N3, design=des, method="ls", block=animal, cor=cc)

cmat <- makeContrasts(VI0 = VENTP0 - INTERP0,
                      VD0 = VENTP0 - DORSALP0,
                      ID0 = INTERP0 - DORSALP0,
                      VI9 = VENTP9 - INTERP9,
                      VD9 = VENTP9 - DORSALP9,
                      ID9 = INTERP9 - DORSALP9,
                      VI18 = VENTP18 - INTERP18,
                      VD18 = VENTP18 - DORSALP18,
                      ID18 = INTERP18 - DORSALP18,
                      VI60 = VENTP60 - INTERP60,
                      VD60 = VENTP60 - DORSALP60,
                      ID60 = INTERP60 - DORSALP60,
                      DiffP0vsP9 = (VENTP9 - VENTP0) - (DORSALP9 - DORSALP0),
                      DiffP9vsP18 = (VENTP18 - VENTP9) - (DORSALP18 - DORSALP9),
                      DiffP18vsP60 = (VENTP60 - VENTP18) - (DORSALP60 - DORSALP18),
                      VIDiffP0vsP9 = (VENTP9 - VENTP0) - (INTERP9 - INTERP0),
                      VIDiffP9vsP18 = (VENTP18 - VENTP9) - (INTERP18 - INTERP9),
                      VIDiffP18vsP60 = (VENTP60 - VENTP18) - (INTERP60 - INTERP18),
                      IDDiffP0vsP9 = (INTERP9 - INTERP0) - (DORSALP9 - DORSALP0),
                      IDDiffP9vsP18 = (INTERP18 - INTERP9) - (DORSALP18 - DORSALP9),
                      IDDiffP18vsP60 = (INTERP60 - INTERP18) - (DORSALP60 - DORSALP18),
                      P0vsP9 = (VENTP9 + INTERP9 + DORSALP9)/3 - (VENTP0 + INTERP0 + DORSALP0)/3,
                      P9vsP18 = (VENTP18 + INTERP18 + DORSALP18)/3 - (VENTP9 + INTERP9 + DORSALP9)/3,
                      P18vsP60 = (VENTP60 + INTERP60 + DORSALP60)/3 - (VENTP18 + INTERP18 + DORSALP18)/3,
                      VD = (VENTP0 + VENTP9 + VENTP18 + VENTP60)/4 - (DORSALP0 + DORSALP9 + DORSALP18 + DORSALP60)/4,
                      VI = (VENTP0 + VENTP9 + VENTP18 + VENTP60)/4 - (INTERP0 + INTERP9 + INTERP18 + INTERP60)/4,
                      ID = (INTERP0 + INTERP9 + INTERP18 + INTERP60)/4 - (DORSALP0 + DORSALP9 + DORSALP18 + DORSALP60)/4,
                      VP0vsP9 = VENTP9 - VENTP0, VP9vsP18 = VENTP18 - VENTP9, VP18vsP60 = VENTP60 - VENTP18,
                      IP0vsP9 = INTERP9 - INTERP0, IP9vsP18 = INTERP18 - INTERP9, IP18vsP60 = INTERP60 - INTERP18,
                      DP0vsP9 = DORSALP9 - DORSALP0, DP9vsP18 = DORSALP18 - DORSALP9, DP18vsP60 = DORSALP60 - DORSALP18,
                      levels=des)

fit2 <- contrasts.fit(fit, cmat)
fit2 <- eBayes(fit2)

##Get results with topTable for all contrasts. Create an annotation conversion table between probeID, EntrezID, old gene symbol (from chip), and new gene symbol (from org.Rn.eg.db using EntrezID) 
a.dc <- decideTests(fit2, adjust.method="fdr", p.value=0.05, lfc=1)
tt.coef1 <- topTable(fit2, coef=1, sort.by='none', n=Inf)
tt.coef1.sel <- tt.coef1[, c("ProbeID", "PROBE_ID", "AveExpr", "logFC", "adj.P.Val")]
probes <- as.character(tt.coef1.sel$PROBE_ID)
EntrezID <- unlist(mget(probes, illuminaRatv1ENTREZID, ifnotfound=NA))
GeneName <- unlist(mget(probes, illuminaRatv1GENENAME, ifnotfound=NA))
GeneID <- unlist(mget(probes, illuminaRatv1SYMBOL, ifnotfound=NA))
conversion_raw.df <- data.frame(EntrezID, GeneName, GeneID, ProbeID=tt.coef1.sel$ProbeID)
conversion_raw.df$PROBE_ID <- rownames(conversion_raw.df)
entrez_symbol <- as.data.frame(org.Rn.egSYMBOL)
colnames(entrez_symbol) <- c("EntrezID", "GeneID_updated")
conversion_raw.df2 <- merge(conversion_raw.df, entrez_symbol, all.x = TRUE)
conversion.df <- conversion_raw.df2[match(tt.coef1.sel$PROBE_ID, conversion_raw.df2$PROBE_ID),]

tt <- cbind(conversion.df, tt.coef1.sel)
colnames(tt) <- sub("adj.P.Val", "adj.P.Val V0 - I0", sub("logFC", "logFC V0 - I0", colnames(tt)))

contrasts <- colnames(fit2$coefficients)
names(contrasts) <- contrasts

coef.names <- contrasts[2:length(contrasts)]
for (n in coef.names){
  tt.n <- topTable(fit2, coef=n, sort.by='none', n=Inf)
  tt.n <- tt.n[, c("logFC", "adj.P.Val")]
  colnames(tt.n) <- paste(c("logFC", "adj.P.Val"), n)
  tt <- cbind(tt, tt.n)
}


saveRDS(data.N3, file = "data.N3.rds")
saveRDS(fit2, file = "fit2.rds")

