## ----setup, echo=FALSE---------------------------------------------------
knitr::opts_chunk$set(fig.align="center", cache=TRUE, message=FALSE, echo=TRUE, results="markup", fig.show="asis", size="small", warning=FALSE, tidy=TRUE, tidy.opts=list(width.cutoff=80))
options(width=65)
library(BiocStyle)

## ----preliminaries-------------------------------------------------------
library(limma)
library(edgeR)
library(EDASeq)
library(RUVSeq)
library(ffpe)
library(RColorBrewer)

## ----readFC--------------------------------------------------------------
fc <- read.table("Peixoto_Additional_inputtext/Peixoto_CC_FC_RT.txt", row.names=1, header=TRUE)
negControls <- read.table("Peixoto_Additional_inputtext/Peixoto_NegativeControls.txt", sep='\t', header=TRUE, as.is=TRUE)
positive <- read.table("Peixoto_Additional_inputtext/Peixoto_positive_controls.txt", as.is=TRUE, sep='\t', header=TRUE)

x <- as.factor(rep(c("CC", "FC", "RT"), each=5))
names(x) <- colnames(fc)

filter <- apply(fc, 1, function(x) length(x[which(x>10)])>5)
filtered <- as.matrix(fc)[filter,]

negCon <- intersect(negControls[,2], rownames(filtered))
FCup <- intersect(positive[positive[,3]=="UP",1], rownames(filtered))
FCdown <- intersect(positive[positive[,3]=="DOWN",1], rownames(filtered))
RTup <- intersect(positive[positive[,4]=="UP",1], rownames(filtered))
RTdown <- intersect(positive[positive[,4]=="DOWN",1], rownames(filtered))

colors <- brewer.pal(9, "Set1")
colLib <- colors[x]

## ----uq------------------------------------------------------------------
uq <- betweenLaneNormalization(filtered, which="upper")

plotRLE(uq, col=colLib, outline=FALSE, las=3, ylim=c(-.2, .2), ylab="Relative Log Expression", cex.axis=1, cex.lab=1)
plotPCA(uq, col=colLib, cex=1, cex.axis=1, cex.lab=1, xlim=c(-.6, .9), ylim=c(-.7, .6))

## ----replicateMatrix-----------------------------------------------------
groups <- matrix(data=c(1:5, 6:10, 11:15), nrow=3, byrow=TRUE)
groups

## ----ruv-----------------------------------------------------------------
s <- RUVs(uq, negCon, k=5, groups)

plotRLE(s$normalizedCounts, col=colLib, outline=FALSE, las=3, ylim=c(-.2, .2), ylab="Relative Log Expression", cex.axis=1, cex.lab=1)
plotPCA(s$normalizedCounts, col=colLib, cex=1, cex.axis=1, cex.lab=1, xlim=c(-.6, .9), ylim=c(-.7, .6))

## ----uqDesign------------------------------------------------------------
design <- model.matrix(~x)
design

## ----edgerUQ-------------------------------------------------------------
y <- DGEList(counts=filtered, group=x)
y <- calcNormFactors(y, method="upperquartile")

y <- estimateGLMCommonDisp(y, design, verbose=TRUE)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)

## ----uq-fc---------------------------------------------------------------
lrt <- glmLRT(fit, coef=2)
topUQFC <- topTags(lrt, n=Inf)$table

hist(topUQFC$PValue, main="", xlab="p-value", breaks=100, ylim=c(0, 1400))

## ----uq-fc-volcano-------------------------------------------------------
plot(topUQFC[,1], -log10(topUQFC$PValue), pch=20, col="gray", cex=.5, ylab="-log10(p-value)", xlab="log2(FC/CC)", ylim=c(0, 85), xlim=c(-2, 4), cex.lab=1, cex.axis=1)
de <- rownames(topUQFC[topUQFC$FDR<=0.01,])
points(topUQFC[de,1], -log10(topUQFC[de, "PValue"]), pch=20, col=colors[2], cex=1)
points(topUQFC[FCup,1], -log10(topUQFC[FCup, "PValue"]), pch=1, col=colors[1], cex=1, lwd=2)
points(topUQFC[FCdown,1], -log10(topUQFC[FCdown, "PValue"]), pch=1, col=colors[1], cex=1, lwd=2)
points(topUQFC[negCon,1], -log10(topUQFC[negCon, "PValue"]), pch=1, col=colors[3], cex=1, lwd=2)

## ----uq-rt---------------------------------------------------------------
lrt <- glmLRT(fit, coef=3)
topUQRT <- topTags(lrt, n=Inf)$table

hist(topUQRT$PValue, main="", xlab="p-value", breaks=100, ylim=c(0, 900))

plot(topUQRT[,1], -log10(topUQRT$PValue), pch=20, col="gray", cex=.5, ylab="-log10(p-value)", xlab="log2(RT/CC)", ylim=c(0, 65), xlim=c(-2, 2.5), cex.lab=1, cex.axis=1)
de <- rownames(topUQRT[topUQRT$FDR<=0.01,])
points(topUQRT[de,1], -log10(topUQRT[de, "PValue"]), pch=20, col=colors[2], cex=1)
points(topUQRT[RTup,1], -log10(topUQRT[RTup, "PValue"]), pch=1, col=colors[1], cex=1, lwd=2)
points(topUQRT[RTdown,1], -log10(topUQRT[RTdown, "PValue"]), pch=1, col=colors[1], cex=1, lwd=2)
points(topUQRT[negCon,1], -log10(topUQRT[negCon, "PValue"]), pch=1, col=colors[3], cex=1, lwd=2)

## ----designRUV-----------------------------------------------------------
design <- model.matrix(~x + s$W)
design

## ----edgerRUV------------------------------------------------------------
y <- DGEList(counts = filtered, group = x)
y <- calcNormFactors(y, method = "upperquartile")

y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)

## ----ruv-fc--------------------------------------------------------------
lrt <- glmLRT(fit, coef=2)
topRsFC <- topTags(lrt, n=Inf)$table

hist(topRsFC$PValue, main="", xlab="p-value", breaks=100, ylim=c(0, 1400))

plot(topRsFC[,1], -log10(topRsFC$PValue), pch=20, col="gray", cex=.5, ylab="-log10(p-value)", xlab="log2(FC/CC)", ylim=c(0, 85), xlim=c(-2, 4), cex.lab=1, cex.axis=1)
de <- rownames(topRsFC[topRsFC$FDR<=0.01,])
points(topRsFC[de,1], -log10(topRsFC[de, "PValue"]), pch=20, col=colors[2], cex=1)
points(topRsFC[FCup,1], -log10(topRsFC[FCup, "PValue"]), pch=1, col=colors[1], cex=1, lwd=2)
points(topRsFC[FCdown,1], -log10(topRsFC[FCdown, "PValue"]), pch=1, col=colors[1], cex=1, lwd=2)
points(topRsFC[negCon,1], -log10(topRsFC[negCon, "PValue"]), pch=1, col=colors[3], cex=1, lwd=2)

## ----ruv-rt--------------------------------------------------------------
lrt <- glmLRT(fit, coef=3)
topRsRT <- topTags(lrt, n=Inf)$table

hist(topRsRT$PValue, main="", xlab="p-value", breaks=100, ylim=c(0, 900))

plot(topRsRT[,1], -log10(topRsRT$PValue), pch=20, col="gray", cex=.5, ylab="-log10(p-value)", xlab="log2(RT/CC)", ylim=c(0, 65), xlim=c(-2, 2.5), cex.lab=1, cex.axis=1)
de <- rownames(topRsRT[topRsRT$FDR<=0.01,])
points(topRsRT[de,1], -log10(topRsRT[de, "PValue"]), pch=20, col=colors[2], cex=1)
points(topRsRT[RTup,1], -log10(topRsRT[RTup, "PValue"]), pch=1, col=colors[1], cex=1, lwd=2)
points(topRsRT[RTdown,1], -log10(topRsRT[RTdown, "PValue"]), pch=1, col=colors[1], cex=1, lwd=2)
points(topRsRT[negCon,1], -log10(topRsRT[negCon, "PValue"]), pch=1, col=colors[3], cex=1, lwd=2)

## ----olm-----------------------------------------------------------------
olm <- read.table("Peixoto_Additional_inputtext/Peixoto_OLM_HC.txt", row.names=1, header=TRUE)
stopifnot(all(rownames(olm)==rownames(fc)))

x <- as.factor(rep(c("HC", "OLM"), each=6))
names(x) <- colnames(olm)
colLib <- colors[x]

filter <- apply(olm, 1, function(x) length(x[which(x>10)])>5)
filtered <- as.matrix(olm[filter,])

negCon <- intersect(negControls[,2], rownames(filtered))
OLMup <- intersect(positive[positive[,5]=="UP",1], rownames(filtered))
OLMdown <- intersect(positive[positive[,5]=="DOWN",1], rownames(filtered))

## ----uq-olm--------------------------------------------------------------
uqOLM <- betweenLaneNormalization(filtered, which="upper")

plotRLE(uqOLM, col=colLib, outline=FALSE, las=3, ylim=c(-.2, .2), ylab="Relative Log Expression", cex.axis=1, cex.lab=1)

plotPCA(uqOLM, col=colLib, cex=1, cex.axis=1, cex.lab=1, xlim=c(-.7, .7), ylim=c(-.7, .7))

## ----ruv-olm-------------------------------------------------------------
groups <- matrix(data=c(1:6, 7:12), nrow=2, byrow=TRUE)

sOLM <- RUVs(uqOLM, negCon, k=4, groups)

plotRLE(sOLM$normalizedCounts, col=colLib, outline=FALSE, las=3, ylim=c(-.2, .2), ylab="Relative Log Expression", cex.axis=1, cex.lab=1)

plotPCA(sOLM$normalizedCounts, col=colLib, cex=1, cex.axis=1, cex.lab=1, xlim=c(-.7, .7), ylim=c(-.7, .7))

## ----edger-uq-olm--------------------------------------------------------
design <- model.matrix(~x)
y <- DGEList(counts=filtered, group=x)
y <- calcNormFactors(y, method="upperquartile")

y <- estimateGLMCommonDisp(y, design, verbose=TRUE)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)

lrt <- glmLRT(fit, coef=2)
topUQOLM <- topTags(lrt, n=Inf)$table

hist(topUQOLM$PValue, main="", xlab="p-value", breaks=100, ylim=c(0, 600))

plot(topUQOLM[,1], -log10(topUQOLM$PValue), pch=20, col="gray", cex=.5, ylab="-log10(p-value)", xlab="log2(OLM/HC)", ylim=c(0, 80), xlim=c(-3, 3), cex.lab=1, cex.axis=1)
de <- rownames(topUQOLM[topUQOLM$FDR<=0.01,])
points(topUQOLM[de,1], -log10(topUQOLM[de, "PValue"]), pch=20, col=colors[2], cex=1)
points(topUQOLM[OLMup,1], -log10(topUQOLM[OLMup, "PValue"]), pch=1, col=colors[1], cex=1, lwd=2)
points(topUQOLM[OLMdown,1], -log10(topUQOLM[OLMdown, "PValue"]), pch=1, col=colors[1], cex=1, lwd=2)
points(topUQOLM[negCon,1], -log10(topUQOLM[negCon, "PValue"]), pch=1, col=colors[3], cex=1, lwd=2)

## ----edger-ruv-olm-------------------------------------------------------
design <- model.matrix(~x + sOLM$W)

y <- DGEList(counts = filtered, group = x)
y <- calcNormFactors(y, method = "upperquartile")

y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)

lrt <- glmLRT(fit, coef = 2)
topSOLM <- topTags(lrt, n=Inf)$table

hist(topSOLM$PValue, main="", xlab="p-value", breaks=100, ylim=c(0, 600))

plot(topSOLM[,1], -log10(topSOLM$PValue), pch=20, col="gray", cex=.5, ylab="-log10(p-value)", xlab="log2(OLM/HC)", ylim=c(0, 80), xlim=c(-3, 3), cex.lab=1, cex.axis=1)
de <- rownames(topSOLM[topSOLM$FDR<=0.01,])
points(topSOLM[de,1], -log10(topSOLM[de, "PValue"]), pch=20, col=colors[2], cex=1)
points(topSOLM[OLMup,1], -log10(topSOLM[OLMup, "PValue"]), pch=1, col=colors[1], cex=1, lwd=2)
points(topSOLM[OLMdown,1], -log10(topSOLM[OLMdown, "PValue"]), pch=1, col=colors[1], cex=1, lwd=2)
points(topSOLM[negCon,1], -log10(topSOLM[negCon, "PValue"]), pch=1, col=colors[3], cex=1, lwd=2)

## ----combined------------------------------------------------------------
counts <- cbind(olm, fc[,!grepl("RT", colnames(fc))])
batch <- as.factor(c(rep("lab1", 12), rep("lab2", 10)))
x <- as.factor(c(rep(c("CC", "OLM"), each=6), rep(c("CC", "FC"), each=5)))
names(x) <- names(batch) <- colnames(counts)
colors <- brewer.pal(9, "Set1")
colLib <- colors[x]
colBatch <- colors[4:5][batch]

filter <- apply(counts, 1, function(x) length(x[which(x>10)])>5)
filtered <- as.matrix(counts[filter,])
negCon <- intersect(negControls[,2], rownames(filtered))

## ----uq-combined---------------------------------------------------------
uqCombined <- betweenLaneNormalization(filtered, which="upper")

plotRLE(uqCombined, col=colLib, outline=FALSE, las=3, ylim=c(-.2, .2), ylab="Relative Log Expression", cex.axis=1, cex.lab=1)

plotPCA(uqCombined, col=colLib, cex=1, cex.axis=1, cex.lab=1, xlim=c(-.4, .6), ylim=c(-.4, .5))

## ----ruv-combined--------------------------------------------------------
groups <- matrix(data=c(1:6, 13:17, 7:12, rep(-1, 5), 18:22, rep(-1, 6)), nrow=3, byrow=TRUE)

sCombined <- RUVs(uqCombined, negCon, k=6, groups)

plotRLE(sCombined$normalizedCounts, col=colLib, outline=FALSE, las=3, ylim=c(-.2, .2), ylab="Relative Log Expression", cex.axis=1, cex.lab=1)

plotPCA(sCombined$normalizedCounts, col=colLib, cex=1, cex.axis=1, cex.lab=1, xlim=c(-.4, .6), ylim=c(-.4, .5))

## ----edger-uq-combined---------------------------------------------------
design <- model.matrix(~x -1)
y <- DGEList(counts=filtered, group=x)
y <- calcNormFactors(y, method="upperquartile")

y <- estimateGLMCommonDisp(y, design, verbose=TRUE)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)

## (OLM+FC)/2 vs. CC
lrt <- glmLRT(fit, contrast=c(-1, 1/2, 1/2))

## FC vs. CC

## lrt <- glmLRT(fit, contrast=c(-1, 1, 0))

## OLM vs. CC

## lrt <- glmLRT(fit, contrast=c(-1, 0, 1))

## FC vs. OLM

## lrt <- glmLRT(fit, contrast=c(0, 1, -1))

topUQCombined <- topTags(lrt, n=Inf)$table

hist(topUQCombined$PValue, main="", xlab="p-value", breaks=100, ylim=c(0, 1200))

plot(topUQCombined[,1], -log10(topUQCombined$PValue), pch=20, col="gray", cex=.5, ylab="-log10(p-value)", xlab="log2((OLM+FC)/2/CC)", ylim=c(0, 80), xlim=c(-3, 3), cex.lab=1, cex.axis=1)
de <- rownames(topUQCombined[topUQCombined$FDR<=0.01,])
points(topUQCombined[de,1], -log10(topUQCombined[de, "PValue"]), pch=20, col=colors[2], cex=1)
points(topUQCombined[FCup,1], -log10(topUQCombined[FCup, "PValue"]), pch=1, col=colors[1], cex=1, lwd=2)
points(topUQCombined[FCdown,1], -log10(topUQCombined[FCdown, "PValue"]), pch=1, col=colors[1], cex=1, lwd=2)
points(topUQCombined[OLMup,1], -log10(topUQCombined[OLMup, "PValue"]), pch=1, col=colors[1], cex=1, lwd=2)
points(topUQCombined[OLMdown,1], -log10(topUQCombined[OLMdown, "PValue"]), pch=1, col=colors[1], cex=1, lwd=2)
points(topUQCombined[negCon,1], -log10(topUQCombined[negCon, "PValue"]), pch=1, col=colors[3], cex=1, lwd=2)


## ----edger-ruv-combined--------------------------------------------------
design <- model.matrix(~x + sCombined$W -1)

y <- DGEList(counts = filtered, group = x)
y <- calcNormFactors(y, method = "upperquartile")

y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)

## (OLM+FC)/2 vs. CC
lrt <- glmLRT(fit, contrast=c(-1, 1/2, 1/2, rep(0, 6)))

## FC vs. CC

## lrt <- glmLRT(fit, contrast=c(-1, 1, 0, rep(0, 6)))

## OLM vs. CC

## lrt <- glmLRT(fit, contrast=c(-1, 0, 1, rep(0, 6)))

## FC vs. OLM

## lrt <- glmLRT(fit, contrast=c(0, 1, -1, rep(0, 6)))

topSCombined <- topTags(lrt, n=Inf)$table

hist(topSCombined$PValue, main="", xlab="p-value", breaks=100, ylim=c(0, 1200))

plot(topSCombined[,1], -log10(topSCombined$PValue), pch=20, col="gray", cex=.5, ylab="-log10(p-value)", xlab="log2((OLM+FC)/2/CC)", ylim=c(0, 80), xlim=c(-3, 3), cex.lab=1, cex.axis=1)
de <- rownames(topSCombined[topSCombined$FDR<=0.01,])
points(topSCombined[de,1], -log10(topSCombined[de, "PValue"]), pch=20, col=colors[2], cex=1)
points(topSCombined[FCup,1], -log10(topSCombined[FCup, "PValue"]), pch=1, col=colors[1], cex=1, lwd=2)
points(topSCombined[FCdown,1], -log10(topSCombined[FCdown, "PValue"]), pch=1, col=colors[1], cex=1, lwd=2)
points(topSCombined[OLMup,1], -log10(topSCombined[OLMup, "PValue"]), pch=1, col=colors[1], cex=1, lwd=2)
points(topSCombined[OLMdown,1], -log10(topSCombined[OLMdown, "PValue"]), pch=1, col=colors[1], cex=1, lwd=2)
points(topSCombined[negCon,1], -log10(topSCombined[negCon, "PValue"]), pch=1, col=colors[3], cex=1, lwd=2)


## ----microarray----------------------------------------------------------
data <- read.table("Peixoto_Additional_inputtext/Peixoto_FC_array_combined.txt", header=TRUE, as.is=TRUE, row.names=74)
data <- data[,-1]

platform <- rep("seq", ncol(data))
platform[grep("_", colnames(data))] <- "array"
platform <- as.factor(platform)
names(platform) <- colnames(data)

x <- substr(colnames(data), 1, 2)
x[x=="TT"] <- "RT"
x <- as.factor(x)
names(x) <- colnames(data)

## subset
include <- c(paste("CC30", c(3, 5, 6, 7), sep="_"), paste("FC30", c(3, 5, 6, 7), sep="_"), paste("CC", c(3, 5, 6, 7), sep=""), paste("FC", c(3, 5, 6, 7), sep=""))
data <- data[, include]
x <- droplevels(x[include])
platform <- droplevels(platform[include])

### filter
array <- as.matrix(data[, platform=="array"])
seq <- as.matrix(data[, platform=="seq"])
xA <- x[platform=="array"]
xS <- x[platform=="seq"]
filter <- apply(seq, 1, function(x) length(x[x>5])>3)

array <- array[filter,]
seq <- seq[filter,]

## ----limma---------------------------------------------------------------
design <- model.matrix(~xA)
fit <- lmFit(array, design)
fit <- eBayes(fit)
top <- topTable(fit, coef=2, n=Inf)
deLimma <- rownames(top)[top$adj.P.Val<=0.1]

## ----arrayUQ-------------------------------------------------------------
design <- model.matrix(~xS)
y <- DGEList(counts=seq, group=xS)
y <- calcNormFactors(y, method="upperquartile")

y <- estimateGLMCommonDisp(y, design, verbose=TRUE)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
topUQ <- topTags(lrt, n=nrow(seq))$table
deUQ <- rownames(topUQ[topUQ$FDR<0.1,])

## ----arrayRUV------------------------------------------------------------
negCon <- intersect(negControls[,2], rownames(seq))
norm <- betweenLaneNormalization(seq, which="upper")
groups <- matrix(data=c(1:4, 5:8), nrow=2, byrow=TRUE)
rS <- RUVs(norm, negCon, k=1, groups)

design <- model.matrix(~xS + rS$W)
y <- DGEList(counts=seq, group=xS)
y <- calcNormFactors(y, method="upperquartile")

y <- estimateGLMCommonDisp(y, design, verbose=TRUE)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
topRS <- topTags(lrt, n=nrow(seq))$table
deRS <- rownames(topRS[topRS$FDR<0.1,])

## ----comparison----------------------------------------------------------
limmaP <- top$P.Value
names(limmaP) <- rownames(top)

uqP <- topUQ$PValue
names(uqP) <- rownames(topUQ)

ruvSP <- topRS$PValue
names(ruvSP) <- rownames(topRS)

uq_limma = CATplot(uqP, limmaP, maxrank=1000, make.plot=F)
ruvS_limma = CATplot(ruvSP, limmaP, maxrank=1000, make.plot=F)


ul <- uq_limma[1:500,]
rl <- ruvS_limma[1:500,]
plot(ul[-(1:20),], ylim=c(0.15,0.55), col=colors[1], lwd=2, type="l", cex.axis=1, cex.lab=1)
lines(rl[-(1:20),], col=colors[2], lwd=2)
legend("topright", legend=c("UQ vs. limma", "RUVs vs limma"), col=colors, lwd=2, cex=1)

## ----sessionInfo---------------------------------------------------------
sessionInfo()

