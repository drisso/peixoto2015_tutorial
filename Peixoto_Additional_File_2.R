## ----setup, echo=FALSE, message=FALSE------------------------------------
knitr::opts_chunk$set(fig.align="center", cache=TRUE, message=FALSE, echo=TRUE, results="markup", fig.show="asis", size="small", warning=FALSE, tidy=TRUE, tidy.opts=list(width.cutoff=80))
options(width=65)
library(BiocStyle)
library(EDASeq)
library(RUVSeq)
library(RColorBrewer)
colors <- brewer.pal(9, "Set1")


## ----Barco_mrna----------------------------------------------------------
data <- read.table("Peixoto_Input_for_Additional_file_2/GSE60261.txt", header = TRUE, sep="\t")
dicerKOmRNA <- as.matrix(data[,2:ncol(data)])
x <- as.factor(rep(c("KO", "WT"), each = 3))
colLib <- colors[x]

boxplot(log(dicerKOmRNA+1), las=2)

filter <- apply(dicerKOmRNA, 1, function(x) length(x[which(x>0)])>3)
table(filter)
filtered <- dicerKOmRNA[filter,]

plotRLE(filtered, col=colLib, outline=FALSE, ylim=c(-.4, .4), las=2)
plotPCA(filtered, col=colLib)

groups <- matrix(data=c(1:3, 4:6), nrow=2, byrow=TRUE)
s <- RUVs(round(filtered), 1:nrow(filtered), k=1, groups)

plotRLE(s$normalizedCounts, col=colLib, outline=FALSE, ylim=c(-.4, .4), las=2)
plotPCA(s$normalizedCounts, col=colLib)

## ----Barco_small---------------------------------------------------------
data <- read.table("Peixoto_Input_for_Additional_file_2/GSE60262.txt", header = TRUE, sep="\t")
dicerKOsmallRNA <- as.matrix(data[,2:ncol(data)])
x <- as.factor(rep(c("KO", "WT"), each = 3))
colLib <- colors[x]

boxplot(log(dicerKOsmallRNA+1), las=2)

filter <- apply(dicerKOsmallRNA, 1, function(x) length(x[which(x>0)])>3)
table(filter)
filtered <- dicerKOsmallRNA[filter,]

plotRLE(filtered, col=colLib, outline=FALSE, ylim=c(-1, 1), las=2)
plotPCA(filtered, col=colLib)

groups <- matrix(data=c(1:3, 4:6), nrow=2, byrow=TRUE)
s <- RUVs(round(filtered), 1:nrow(filtered), k=2, groups, round=FALSE)

plotRLE(s$normalizedCounts, col=colLib, outline=FALSE, ylim=c(-1, 1), las=2)
plotPCA(s$normalizedCounts, col=colLib)

## ----h2az----------------------------------------------------------------
data <- read.table("Peixoto_Input_for_Additional_file_2/GSE58797.txt", header = TRUE, sep="\t")
H2AZRNA <- as.matrix(data[,2:10])
x <- as.factor(rep(c("Control", "shRNA", "FC"), each = 3))
colLib <- colors[x]

boxplot(log(H2AZRNA+1), las=2)

filter <- apply(H2AZRNA, 1, function(x) length(x[which(x>0)])>3)
table(filter)
filtered <- H2AZRNA[filter,]

H2AZFPKMuq <- betweenLaneNormalization(filtered, which="upper", round=FALSE)

boxplot(log(H2AZFPKMuq+1), las=2)

plotRLE(H2AZFPKMuq, col=colLib, outline=FALSE, ylim=c(-.2, .2), las=2)

plotPCA(H2AZFPKMuq, col=colLib)

groups <- matrix(data=c(1:3, 4:6, 7:9), nrow=3, byrow=TRUE)
s <- RUVs(round(H2AZFPKMuq), 1:nrow(H2AZFPKMuq), k=1, groups, round=FALSE)

plotRLE(s$normalizedCounts, col=colLib, outline=FALSE, ylim=c(-.2, .2), las=2)
plotPCA(s$normalizedCounts, col=colLib)

## ----Fischer1------------------------------------------------------------
data <- read.table("Peixoto_Input_for_Additional_file_2/GSE61915.txt", header = TRUE, sep="\t")
agingRNASE <- as.matrix(data[,2:ncol(data)])
x <- as.factor(c(rep("Young",5), rep("Old",6)))
colLib <- colors[x]

boxplot(log(agingRNASE+1), las=2)

filter <- apply(agingRNASE, 1, function(x) length(x[which(x>0)])>3)
table(filter)
filtered <- agingRNASE[filter,]

plotRLE(filtered, col=colLib, outline=FALSE, ylim=c(-.5, .5), las=2)
plotPCA(filtered, col=colLib)

groups <- matrix(data=c(1:5, -1, 6:11), nrow=2, byrow=TRUE)
s <- RUVs(round(filtered), 1:nrow(filtered), k=1, groups)

plotRLE(s$normalizedCounts, col=colLib, outline=FALSE, ylim=c(-.5, .5), las=2)
plotPCA(s$normalizedCounts, col=colLib)

## ----Fischer2------------------------------------------------------------
data <- read.table("Peixoto_Input_for_Additional_file_2/GSE53380.txt", header = TRUE, sep="\t")
KO_NORRNA <- as.matrix(data[,2:ncol(data)])
x <- as.factor(c(rep("Control", 6), rep("NOR", 5), rep("KO",6), rep("KO_NOR",7)))
colLib <- colors[x]

boxplot(log(KO_NORRNA+1), las=2)

filter <- apply(KO_NORRNA, 1, function(x) length(x[which(x>0)])>3)
table(filter)
filtered <- KO_NORRNA[filter,]

uq <- betweenLaneNormalization(filtered, which="upper")
boxplot(log(uq+1), las=2)

plotRLE(uq, col=colLib, outline=FALSE, ylim=c(-.6, .6), las=2)
plotPCA(uq, col=colLib)

groups <- matrix(data=c(1:6, -1, 7:11, -1, -1, 12:17, -1, 18:24), nrow=4, byrow=TRUE)
s <- RUVs(uq, 1:nrow(uq), k=2, groups)

plotRLE(s$normalizedCounts, col=colLib, outline=FALSE, ylim=c(-.6, .6), las=2)
plotPCA(s$normalizedCounts, col=colLib)

## ----tsai----------------------------------------------------------------
data <- read.table("Peixoto_Input_for_Additional_file_2/GSE65159.txt", header = TRUE, sep="\t")
AD_RNA <- as.matrix(data[,2:ncol(data)])
x <- as.factor(rep(c("Control6wk","AD6wk", "Control2wk","AD2wk"), each = 3))
colLib <- colors[x]

boxplot(log(AD_RNA+1), las=2)

filter <- apply(AD_RNA, 1, function(x) length(x[which(x>0)])>3)
table(filter)
filtered <- AD_RNA[filter,]

uq <- betweenLaneNormalization(filtered, which="upper")
boxplot(log(uq+1), las=2)

plotRLE(uq, col=colLib, outline=FALSE, ylim=c(-1, 1), las=2)
plotPCA(uq, col=colLib)

groups <- matrix(data=c(1:3, 4:6, 7:9, 10:12), nrow=4, byrow=TRUE)
s <- RUVs(uq, 1:nrow(uq), k=1, groups)

plotRLE(s$normalizedCounts, col=colLib, outline=FALSE, ylim=c(-1, 1), las=2)
plotPCA(s$normalizedCounts, col=colLib)

## ----reijmers------------------------------------------------------------
data <- read.table("Peixoto_Input_for_Additional_file_2/GSE58343.txt", header = TRUE, sep="\t")
Ribosome_RNA <- as.matrix(data[,2:ncol(data)])
x <- as.factor(c(rep("HC", 10), rep("FC", 12)))
colLib <- colors[x]

boxplot(log(Ribosome_RNA+1), las=2)

filter <- apply(Ribosome_RNA, 1, function(x) length(x[which(x>0)])>3)
table(filter)
filtered <- Ribosome_RNA[filter,]

boxplot(log(filtered+1), las=2)

plotRLE(filtered, col=colLib, outline=FALSE, ylim=c(-1, 2), las=2)
plotPCA(filtered, col=colLib)

groups <- matrix(data=c(1, 3, 6, 8,
                        4, 9, -1, -1,
                        2, 5, 7, 10,
                        11, 13, 17, 19,
                        14, 20, -1, -1,
                        16, 22, -1, -1,
                        12, 15, 18, 21), nrow=7, byrow=TRUE)
s <- RUVs(round(filtered), 1:nrow(filtered), k=1, groups, round=FALSE)

plotRLE(s$normalizedCounts, col=colLib, outline=FALSE, ylim=c(-1, 2), las=2)
plotPCA(s$normalizedCounts, col=colLib)

## ----sessionInfo---------------------------------------------------------
sessionInfo()

