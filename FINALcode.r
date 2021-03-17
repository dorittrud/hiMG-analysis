
#Load the required Packages
library(DESeq)
library(limma)
library(edgeR)
library(reshape2)
library(biomaRt)
library(annotate)
library(gplots)
library(RColorBrewer)
library(flashClust)
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
library(Hmisc)
library(ggplot2)

####Set the working directory and load required packages for Analysis in R.
setwd("C:/Users/Kristopher/Dropbox/Microglia/FINAL_ANALYSIS/")

####Previously, raw gene count data was generated using subread and merged into a single data frame. We now read this dataframe into R.
load("C:\\Users\\Kristopher\\Dropbox\\Microglia\\FINAL_ANALYSIS\\RAWdata.RData")
dim(RAWdata)
row.names(RAWdata) = RAWdata$Geneid
RAWdata = RAWdata[order(RAWdata$Geneid, decreasing = FALSE),]
RAWdata = RAWdata[,-1]

#Read in Gene Annots
GeneAnnots = read.delim("AllGeneAnnots.txt", sep="\t", header=TRUE)
row.names(GeneAnnots) = GeneAnnots$Geneid
GeneAnnots = GeneAnnots[order(GeneAnnots$Geneid, decreasing = FALSE),]
table(row.names(RAWdata) == row.names(GeneAnnots))
backup.GeneAnnots = GeneAnnots
GeneAnnots = GeneAnnots[,2]

#Copy Raw data to a new data frame for manipulation
allDATA = RAWdata

#Read in annotations for each data subset
FINAL_ANNOTS.txt = read.delim("FINAL_ANNOTS.txt", sep="\t", header=TRUE)

#make sure that the annotations and data are properly ordered
table(names(allDATA) == FINAL_ANNOTS.txt$Rep)

#add row.names to each sample annotation file
row.names(FINAL_ANNOTS.txt) = FINAL_ANNOTS.txt$Rep

####Prepare DGElist objects for input into Voom
DGEdata.allDATA <- DGEList(counts =  allDATA, group=FINAL_ANNOTS.txt$group, genes=GeneAnnots)

#Filter the data for genes that are expressed in relevant samples and normalize the data sets
X = DGEdata.allDATA
isexprDATA <- rowSums(cpm(X[,1:182]) > 5) >=4;
table(isexprDATA)
X <- X[isexprDATA,keep.lib.sizes=FALSE]
X <- calcNormFactors(X, method = "TMM")
DGEdata.allDATA = X
rm(X)

#From the isexprdata table, we see that 21395 genes pass the filter
nGenes_DGEdata.allDATA = 21395

####Create design matrices for each of the analyses, this specifies the comparisons that we would like to make.
design.allDATA = FINAL_ANNOTS.txt[,-c(1:9)]

####Voom for allData
VOOM.allDATA <- voom(DGEdata.allDATA, design=design.allDATA,plot=TRUE)

####Calculate Differential expression using a global linear model
#specify the comparisons that we want to make for differentially expressed genes
contrasts.hiMG = makeContrasts(
		hiMG.VS.hPSC = hiMG_Trudler_et.al-(hPSC_Trudler_et.al+hPSC_Abud_et.al+hPSC_Gifford_et.al)/3,
		levels = design.allDATA)


####Fit the data to a integrated linear model and calculate differntially expressed genes
FIT.VOOM.allDATA = lmFit(VOOM.allDATA,design.allDATA)
FIT.VOOM.allDATA.v1 = contrasts.fit(FIT.VOOM.allDATA,contrasts.hiMG)
FIT.VOOM.allDATA.v1 = eBayes(FIT.VOOM.allDATA.v1)


#### We use the 'treat method' (McCarthy and Smyth 2009) to calculate p-values from empirical Bayes moderated t-statistics to require genes 
#### to have a log-FC that is significantly greater than 1 (equivalent to a 2-fold difference between cell types on the original scale).
tfit.v1 <- treat(FIT.VOOM.allDATA.v1, lfc=1)
dt.v1 <- decideTests(tfit.v1)
summary(dt.v1)


####extract differential expression statistics for each of the comparisons
DEX.hiMG.VS.hPSC = data.frame(topTreat(tfit.v1,coef="hiMG.VS.hPSC",n=nGenes_DGEdata.allDATA))

#add GeneIds
DEX.hiMG.VS.hPSC$Geneid = row.names(DEX.hiMG.VS.hPSC)

####Generate a Column that specifies which group has higher expression
DEX.hiMG.VS.hPSC$UP_IN = ifelse(DEX.hiMG.VS.hPSC$logFC <= 0,"hPSC","hiMG")
DEX.hiMG.VS.hPSC = DEX.hiMG.VS.hPSC[order(DEX.hiMG.VS.hPSC$UP_IN,DEX.hiMG.VS.hPSC$adj.P.Val, decreasing = FALSE),]
DEX.hiMG.VS.hPSC$SIG = ifelse(abs(DEX.hiMG.VS.hPSC$logFC) >= 1 & DEX.hiMG.VS.hPSC$adj.P.Val <= 0.01,1,0)

#Filter and Organize genes that are up regulated in hiMG vs hPSC 
SIGup.DEX.hiMG.VS.hPSC = DEX.hiMG.VS.hPSC[DEX.hiMG.VS.hPSC$SIG == 1 & DEX.hiMG.VS.hPSC$UP_IN != "hPSC",]
SIGup.DEX.hiMG.VS.hPSC$DirectionalRank = seq(1:nrow(SIGup.DEX.hiMG.VS.hPSC))
SIGup.DEX.hiMG.VS.hPSC = SIGup.DEX.hiMG.VS.hPSC[,c(7,1,2,6,8,10)]

#Filter and Organize genes that are down regulated in hiMG vs hPSC 
SIGdown.DEX.hiMG.VS.hPSC = DEX.hiMG.VS.hPSC[DEX.hiMG.VS.hPSC$SIG == 1 & DEX.hiMG.VS.hPSC$UP_IN == "hPSC",]
SIGdown.DEX.hiMG.VS.hPSC$DirectionalRank = seq(1:nrow(SIGdown.DEX.hiMG.VS.hPSC))
SIGdown.DEX.hiMG.VS.hPSC = SIGdown.DEX.hiMG.VS.hPSC[,c(7,1,2,6,8,10)]

#Recombine the filtered sets and write to file. 
allSIG.DEX.hiMG.VS.hPSC = rbind(SIGup.DEX.hiMG.VS.hPSC,SIGdown.DEX.hiMG.VS.hPSC)
colnames(allSIG.DEX.hiMG.VS.hPSC)[2] = "GeneSymbol"
allSIG.DEX.hiMG.VS.hPSC = merge(allSIG.DEX.hiMG.VS.hPSC, backup.GeneAnnots[,-2], by.x = "Geneid", by.y = "Geneid", sort = FALSE, all.x = TRUE)
write.table(allSIG.DEX.hiMG.VS.hPSC, file = "TableS2A_allSIG.DEX.hiMG.VS.hPSC.txt", sep = "\t", row.names = FALSE, quote = FALSE)

####Create lcpm data for plotting the upregulated genes.
PLOT.lcpm <- data.frame(cpm(DGEdata.allDATA, log=TRUE, prior.count = 2))
table(names(PLOT.lcpm) == FINAL_ANNOTS.txt$Rep)
PLOT.lcpm$Geneid = row.names(PLOT.lcpm)
PLOT.lcpm = merge(backup.GeneAnnots[,1:2], PLOT.lcpm, by.x = "Geneid", by.y = "Geneid", sort = FALSE)
plot.SIGup.DEX.hiMG.VS.hPSC = merge(SIGup.DEX.hiMG.VS.hPSC[,1:2], PLOT.lcpm[,-2], by.x = "Geneid", by.y = "Geneid", sort = FALSE)

#center the lcpm values on the mean for each gene and write to file
names(plot.SIGup.DEX.hiMG.VS.hPSC)[1:2] = c("UNIQID","NAME")
plot.SIGup.DEX.hiMG.VS.hPSC[,3:184] = plot.SIGup.DEX.hiMG.VS.hPSC[,3:184] - apply(plot.SIGup.DEX.hiMG.VS.hPSC[,3:184],1,mean)
write.table(plot.SIGup.DEX.hiMG.VS.hPSC, file = "plot.SIGup.DEX.hiMG.VS.hPSC.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#plot.SIGup.DEX.hiMG.VS.hPSC was clustered using average linkage euclidean distance on genes and arrays in GeneCluster 3.0, and visualized in Java Treeview


####We also want to do functional enrichments for the upregulated genes using GREAT.
GREAT_BACKGROUND.bed = merge(PLOT.lcpm[,1:2], backup.GeneAnnots[,-2], by.x = "Geneid", by.y = "Geneid", sort = FALSE)
GREAT_BACKGROUND.bed = GREAT_BACKGROUND.bed[GREAT_BACKGROUND.bed$Great.Gene == "YES",]
SIGup.DEX.hiMG.VS.hPSC.bed = merge(SIGup.DEX.hiMG.VS.hPSC[,1:2], GREAT_BACKGROUND.bed[,-2], by.x = "Geneid", by.y = "Geneid", sort = FALSE)
GREAT_BACKGROUND.bed = GREAT_BACKGROUND.bed[order(GREAT_BACKGROUND.bed$chr,GREAT_BACKGROUND.bed$Start, decreasing = FALSE),c(3:5,1)]
SIGup.DEX.hiMG.VS.hPSC.bed = SIGup.DEX.hiMG.VS.hPSC.bed[order(SIGup.DEX.hiMG.VS.hPSC.bed$chr,SIGup.DEX.hiMG.VS.hPSC.bed$Start, decreasing = FALSE),c(3:5,1)]
options(scipen=999)
write.table(GREAT_BACKGROUND.bed, file = "GREAT_BACKGROUND.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(SIGup.DEX.hiMG.VS.hPSC.bed, file = "SIGup.DEX.hiMG.VS.hPSC.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
options(scipen=1)


#Enrichements were determined in GREAT using the following settings: 
# Human: GRCh37 (UCSC hg19, Feb/2009)
# Single nearest gene within 1000.0 kb
# DO NOT Include curated regulatory domains


GREATshown_NAMES = c("Ontology","Term.Name","Hyper.Rank","Hyper.Raw.P-Value","Hyper.FDR.Q-Val","Hyper.Fold.Enrichment","Hyper.Foreground.Region.Hits","Hyper.Total.Regions","Hyper.Region.Set.Coverage","Hyper.Foreground.Gene.Hits","Total.Genes.Annotated")
GREATall_NAMES = c("Ontology","ID","Desc","Rank","HyperP","HyperBonfP","HyperFdrQ","RegionFoldEnrich","ExpFgRegions","FgRegionsHit","BgRegionsHit","FgSetTermCov","TermCov","NumFgGenesHit","NumBgGenesHit","TotalGenes","FgRegionNames","BgRegionNames","FgGeneNames","BgGeneNames")

greatExportAll_SIGup.DEX.hiMG.VS.hPSC.tsv = read.delim("greatExportAll_SIGup.DEX.hiMG.VS.hPSC.tsv", sep="\t", header=FALSE,skip=4,col.names=c(GREATall_NAMES))
shown.SIGup.DEX.hiMG.VS.hPSC.tsv = read.delim("temp.txt",sep="\t",header=FALSE,skip=0,col.names=c(GREATshown_NAMES))

shown.SIGup.DEX.hiMG.VS.hPSC.tsv = shown.SIGup.DEX.hiMG.VS.hPSC.tsv[,1:6]
shown.SIGup.DEX.hiMG.VS.hPSC.tsv$UNIQID = paste(shown.SIGup.DEX.hiMG.VS.hPSC.tsv$Ontology,shown.SIGup.DEX.hiMG.VS.hPSC.tsv$Hyper.Rank,sep=";")
greatExportAll_SIGup.DEX.hiMG.VS.hPSC.tsv$UNIQID = paste(greatExportAll_SIGup.DEX.hiMG.VS.hPSC.tsv$Ontology,greatExportAll_SIGup.DEX.hiMG.VS.hPSC.tsv$Rank,sep=";")

greatExportAll_SIGup.DEX.hiMG.VS.hPSC.tsv = greatExportAll_SIGup.DEX.hiMG.VS.hPSC.tsv[,c(21,2,17,19)]
GREATresults_SIGup.DEX.hiMG.VS.hPSC = merge(greatExportAll_SIGup.DEX.hiMG.VS.hPSC.tsv, shown.SIGup.DEX.hiMG.VS.hPSC.tsv, by.x = "UNIQID", by.y = "UNIQID", sort = FALSE)
GREATresults_SIGup.DEX.hiMG.VS.hPSC = GREATresults_SIGup.DEX.hiMG.VS.hPSC[,c(5,6,2,7,9,10,3,4)]

GREATresults_SIGup.DEX.hiMG.VS.hPSC$minuslog10FDRq = -log(GREATresults_SIGup.DEX.hiMG.VS.hPSC$Hyper.FDR.Q.Val,10)
GREATresults_SIGup.DEX.hiMG.VS.hPSC = GREATresults_SIGup.DEX.hiMG.VS.hPSC[,c(1:6,9,7,8)]
write.table(GREATresults_SIGup.DEX.hiMG.VS.hPSC, file = "TableS2B_GREATresults_SIGup.DEX.hiMG.VS.hPSC.txt", sep = "\t", row.names = FALSE, quote = FALSE)


#Now we want to take an unbiased look at the data, so we will generate a sample corrleation matrix using the top 10pct most variable genes in the detected genes set
PLOT.lcpm$VAR <- apply(PLOT.lcpm[,-c(1:2)], 1, var) 
PLOT.lcpm = PLOT.lcpm[order(PLOT.lcpm$VAR, decreasing = TRUE),]

X = round(0.1*nGenes_DGEdata.allDATA)

PLOT.lcpm.top10pct = PLOT.lcpm[1:X,-185]
PLOT.lcpm = PLOT.lcpm[,-185]
names(PLOT.lcpm.top10pct)[1:2] = c("UNIQID","NAME")

PearsonMatrix.Voom.rcorr <- rcorr(as.matrix(PLOT.lcpm.top10pct[,-c(1:2)]), type = "pearson")
PearsonMatrix.Voom.rcorr.r <- data.frame(PearsonMatrix.Voom.rcorr$r)
PearsonMatrix.Voom.rcorr.r$UNIQID = row.names(PearsonMatrix.Voom.rcorr.r)
PearsonMatrix.Voom.rcorr.r$NAME = row.names(PearsonMatrix.Voom.rcorr.r)
PearsonMatrix.Voom.rcorr.r = PearsonMatrix.Voom.rcorr.r[,c(183,184,1:182)]
write.table(PearsonMatrix.Voom.rcorr.r, file = "PearsonMatrix.Voom.rcorr.r.txt", sep = "\t", row.names = FALSE, quote = FALSE)

PLOT.lcpm.top10pct.cntrd = PLOT.lcpm.top10pct
PLOT.lcpm.top10pct.cntrd[,3:184] = PLOT.lcpm.top10pct.cntrd[,3:184] - apply(PLOT.lcpm.top10pct.cntrd[,3:184],1,mean)


PLOT.lcpm.cntrd = PLOT.lcpm 
PLOT.lcpm.cntrd[,3:184] = PLOT.lcpm.cntrd[,3:184] - apply(PLOT.lcpm.cntrd[,3:184],1,mean)
PLOT.lcpm.cntrd
names(PLOT.lcpm.cntrd)[1:2] = c("UNIQID","NAME")

write.table(PLOT.lcpm.cntrd, file = "PLOT.lcpm.cntrd.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#Now we want to determine what genes are in common across the diverse Microglia samples by cmoparing them to cluster 1 samples and non-microglia cluster 2 samples.

contrasts.Signature = makeContrasts(
	allMicroglia.VS.Cluster1 = 
		(pMGL_Muffat_et.al+
		hiMG_Trudler_et.al+
		iMGL_Abud_et.al+
		iMGL.coculture_Abud_et.al+
		iMGL.TGFB.withdrawl_Abud_et.al+
		pMGL.coculture_Muffat_et.al+
		pf.Microglia.coculture_Muffat_et.al+
		pf.Microglia_Muffat_et.al+
		pf.Microglia_Abud_et.al+
		pf.Microglia_Trudler_et.al+
		pa.Microglia_Abud_et.al+
		adult.Microglia.ExVivo_Gosselin_et.al+
		adult.Microglia.InVitro_Gosselin_et.al)/13-(
			adult.Adipose_Encode+
			adult.AdrenalGland_Encode+
			adult.Bladder_Encode+
			adult.Brain_Encode+
			adult.Brain_Gosselin_et.al+
			adult.Colon_Encode+
			adult.Endometrium_Encode+
			adult.Esophagus_Encode+
			adult.Heart_Encode+
			adult.Kidney_Encode+
			adult.Liver_Encode+
			adult.Lung_Encode+
			adult.Pancreas_Encode+
			adult.SalivaryGland_Encode+
			adult.SmIntestine_Encode+
			adult.Thyroid_Encode+
			fetal.Heart_Encode+
			fetal.Lung_Encode+
			fetal.Muscle_Encode+
			fetal.Placenta_Encode+
			fetal.Stomach_Encode+
			fetal.Thymus_Encode+
			hPSC_Abud_et.al+
			hPSC_Trudler_et.al+
			hPSC_Gifford_et.al)/25,
	allMicroglia.VS.Cluster2 = (
		pMGL_Muffat_et.al+
		hiMG_Trudler_et.al+
		iMGL_Abud_et.al+
		iMGL.coculture_Abud_et.al+
		iMGL.TGFB.withdrawl_Abud_et.al+
		pMGL.coculture_Muffat_et.al+
		pf.Microglia.coculture_Muffat_et.al+
		pf.Microglia_Muffat_et.al+
		pf.Microglia_Abud_et.al+
		pf.Microglia_Trudler_et.al+
		pa.Microglia_Abud_et.al+
		adult.Microglia.ExVivo_Gosselin_et.al+
		adult.Microglia.InVitro_Gosselin_et.al)/13-(
			adult.Spleen_Encode+
			adult.Monocytes.ExVivo_Gosselin_et.al+
			adult.LymphNode_Encode+
			pa.Monocytes_Abud_et.al+
			pa.Monocytes.Inflammatory_Abud_et.al+
			pa.MyeloidDendritic_Abud_et.al)/6,
levels = design.allDATA)


####Fit the data to a integrated linear model and calculate differntially expressed genes
FIT.VOOM.allDATA.v2 = contrasts.fit(FIT.VOOM.allDATA,contrasts.Signature)
FIT.VOOM.allDATA.v2 = eBayes(FIT.VOOM.allDATA.v2)


#### We use the 'treat method' (McCarthy and Smyth 2009) to calculate p-values from empirical Bayes moderated t-statistics to require genes to have a log-FC that is significantly greater than 1 (equivalent to a 2-fold difference between cell types on the original scale).
tfit.v2 <- treat(FIT.VOOM.allDATA.v2, lfc=1)
dt.v2 <- decideTests(tfit.v2)
summary(dt.v2)



####extract differential expression statistics for each of the comparisons
DEX.allMicroglia.VS.Cluster1 = data.frame(topTreat(tfit.v2,coef="allMicroglia.VS.Cluster1",n=nGenes_DGEdata.allDATA))
DEX.allMicroglia.VS.Cluster2 = data.frame(topTreat(tfit.v2,coef="allMicroglia.VS.Cluster2",n=nGenes_DGEdata.allDATA))

###Add Geneids
DEX.allMicroglia.VS.Cluster1$Geneid = row.names(DEX.allMicroglia.VS.Cluster1)
DEX.allMicroglia.VS.Cluster2$Geneid = row.names(DEX.allMicroglia.VS.Cluster2)


####Generate a Column that specifies which group has higher expression
DEX.allMicroglia.VS.Cluster1$UP_IN = ifelse(DEX.allMicroglia.VS.Cluster1$logFC <= 0,"Cluster1","allMicroglia")
DEX.allMicroglia.VS.Cluster2$UP_IN = ifelse(DEX.allMicroglia.VS.Cluster2$logFC <= 0,"Cluster2","allMicroglia")

DEX.allMicroglia.VS.Cluster1 = DEX.allMicroglia.VS.Cluster1[order(DEX.allMicroglia.VS.Cluster1$UP_IN,DEX.allMicroglia.VS.Cluster1$adj.P.Val, decreasing = FALSE),]
DEX.allMicroglia.VS.Cluster2 = DEX.allMicroglia.VS.Cluster2[order(DEX.allMicroglia.VS.Cluster2$UP_IN,DEX.allMicroglia.VS.Cluster2$adj.P.Val, decreasing = FALSE),]

DEX.allMicroglia.VS.Cluster1$SIG = ifelse(abs(DEX.allMicroglia.VS.Cluster1$logFC) >= 1 & DEX.allMicroglia.VS.Cluster1$adj.P.Val <= 0.01,1,0)
DEX.allMicroglia.VS.Cluster2$SIG = ifelse(abs(DEX.allMicroglia.VS.Cluster2$logFC) >= 1 & DEX.allMicroglia.VS.Cluster2$adj.P.Val <= 0.01,1,0)


SIGup.DEX.allMicroglia.VS.Cluster1 = DEX.allMicroglia.VS.Cluster1[DEX.allMicroglia.VS.Cluster1$SIG == 1 & DEX.allMicroglia.VS.Cluster1$UP_IN == "allMicroglia",]
SIGup.DEX.allMicroglia.VS.Cluster1$DirectionalRank = seq(1:nrow(SIGup.DEX.allMicroglia.VS.Cluster1))
SIGup.DEX.allMicroglia.VS.Cluster1 = SIGup.DEX.allMicroglia.VS.Cluster1[,c(7,1,2,6,8,10)]

SIGup.DEX.allMicroglia.VS.Cluster2 = DEX.allMicroglia.VS.Cluster2[DEX.allMicroglia.VS.Cluster2$SIG == 1 & DEX.allMicroglia.VS.Cluster2$UP_IN == "allMicroglia",]
SIGup.DEX.allMicroglia.VS.Cluster2$DirectionalRank = seq(1:nrow(SIGup.DEX.allMicroglia.VS.Cluster2))
SIGup.DEX.allMicroglia.VS.Cluster2 = SIGup.DEX.allMicroglia.VS.Cluster2[,c(7,1,2,6,8,10)]
	
	
names(SIGup.DEX.allMicroglia.VS.Cluster1)[3:6] = paste(
		names(SIGup.DEX.allMicroglia.VS.Cluster1)[3:6],
		"vs.Cluster1",sep=".")
names(SIGup.DEX.allMicroglia.VS.Cluster1)[1:2] = c("UNIQID","NAME")

names(SIGup.DEX.allMicroglia.VS.Cluster2)[3:6] = paste(
		names(SIGup.DEX.allMicroglia.VS.Cluster2)[3:6],
		"vs.Cluster2",sep=".")
names(SIGup.DEX.allMicroglia.VS.Cluster2)[1:2] = c("UNIQID","NAME")


FINAL_SIGNATURE = merge(SIGup.DEX.allMicroglia.VS.Cluster1, SIGup.DEX.allMicroglia.VS.Cluster2[,-2], by.x = "UNIQID", by.y = "UNIQID", sort = FALSE)

plot.FINAL_SIGNATURE = merge(FINAL_SIGNATURE[,1:2], PLOT.lcpm[,-2], by.x = "UNIQID", by.y = "Geneid", sort = FALSE)
names(plot.FINAL_SIGNATURE)[1:2] = c("UNIQID","NAME")

plot.FINAL_SIGNATURE[,3:184] = plot.FINAL_SIGNATURE[,3:184] - apply(plot.FINAL_SIGNATURE[,3:184],1,mean)
write.table(plot.FINAL_SIGNATURE, file = "plot.FINAL_SIGNATURE.txt", sep = "\t", row.names = FALSE, quote = FALSE)

plot.FINAL_SIGNATURE2 = plot.FINAL_SIGNATURE[,c(1,2,5,6,8,9,62,63,64,65,66,67,68,69,70,75,76,77,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,115,116,117,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173)]
write.table(plot.FINAL_SIGNATURE2, file = "plot.FINAL_SIGNATURE2.txt", sep = "\t", row.names = FALSE, quote = FALSE)

plot.FINAL_SIGNATURE2c = plot.FINAL_SIGNATURE[,c(1,2,34,35,42,43,103,104,105,106,107,108,109,110,111,112,113,114,174,175,176,177,178,179,180,181,182,183,184)]

temp1 = plot.FINAL_SIGNATURE2
X = temp1[,-c(1:2)]
X[X > 0] <- 1
X[X <= 0] <- 0

temp1$ABOVEmean = apply(X,1,sum)

# 84*.9
# [1] 75.6

temp1 = temp1[temp1$ABOVEmean >= 76,c(1,87)]


temp2 = plot.FINAL_SIGNATURE2c
X = temp2[,-c(1:2)]
X[X <= 0] <- -1
X[X > 0] <- 0
X[X == -1] <- 1
temp2$BELOWmean = apply(X,1,sum)


temp2 = temp2[temp2$BELOWmean >= 24,c(1,30)]


colnames(temp1)[2] = "pctABOVEmean_Microglia"
colnames(temp2)[2] = "pctBELOWmean_Cluster2"

ABoveBelow_90_Signature = merge(temp1, temp2, by.x = "UNIQID", by.y = "UNIQID", sort = TRUE)

ABoveBelow_90_Signature$pctABOVEmean_Microglia = round((ABoveBelow_90_Signature$pctABOVEmean_Microglia/84)*100)
ABoveBelow_90_Signature$pctBELOWmean_Cluster2 = round((ABoveBelow_90_Signature$pctBELOWmean_Cluster2/27)*100)

ABoveBelow_90_Signature = merge(ABoveBelow_90_Signature, FINAL_SIGNATURE, by.x = "UNIQID", by.y = "UNIQID", sort = FALSE)
ABoveBelow_90_Signature = ABoveBelow_90_Signature[ABoveBelow_90_Signature$pctBELOWmean_Cluster2 >= 90,]
ABoveBelow_90_Signature =ABoveBelow_90_Signature[,c(1,4,2:3,5:12)]
ABoveBelow_90_Signature$AVG.RANK = apply(ABoveBelow_90_Signature[,c(6,10)],1,mean)

ABoveBelow_90_Signature = ABoveBelow_90_Signature[order(ABoveBelow_90_Signature$AVG.RANK, decreasing = FALSE),]
ABoveBelow_90_Signature[,1:2]

write.table(SIGup.DEX.allMicroglia.VS.Cluster1, file = "TableS3a.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(SIGup.DEX.allMicroglia.VS.Cluster2, file = "TableS3b.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(ABoveBelow_90_Signature, file = "TableS3c.txt", sep = "\t", row.names = FALSE, quote = FALSE)


plotABoveBelow_90_Signature = merge(ABoveBelow_90_Signature[,c(1:2,13)], plot.FINAL_SIGNATURE[,-2], by.x = "UNIQID", by.y = "UNIQID", sort = FALSE)
plotABoveBelow_90_Signature = plotABoveBelow_90_Signature[order(plotABoveBelow_90_Signature$AVG.RANK, decreasing = FALSE),]
plotFINAL_ABoveBelow_90_Signature = plotABoveBelow_90_Signature[,c(1,2,79,80,81,82,119,120,4,5,123,124,125,126,121,122,72,73,74,75,8,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,11,12,13,14,15,16,17,18,127,128,129,130,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,37,38,39,40,41,42,45,46,35,36,43,44,113,114,115,175,184,185,176,177,178,179,180,181,182,183,104,105,106,107,108,109,110,111,112,6,7,83,84,85,86,87,88,89,90,91,63,64,65,66,67,68,92,93,94,95,96,97,69,70,71,98,99,100,101,102,103,78,116,117,118,76,77,9,10,155,164,165,166,167,168,169,170,171,172,173,156,174,157,158,159,160,161,162,163,131,140,141,142,143,144,145,146,147,148,149,132,150,151,152,153,154,133,134,135,136,137,138,139)]

write.table(plotFINAL_ABoveBelow_90_Signature, file = "plotFINAL_ABoveBelow_90_Signature.txt", sep = "\t", row.names = FALSE, quote = FALSE)



