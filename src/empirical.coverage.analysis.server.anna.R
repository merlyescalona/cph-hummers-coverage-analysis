library(ggplot2)
library(knitr)
suppressWarnings(suppressMessages(library(gplots)))
library(RColorBrewer)
colors = colorRampPalette(brewer.pal(9, "Set1"))(48)
# Helpful functions: http://stackoverflow.com/a/6463946
roundUp <- function(x) 10^ceiling(log10(x))
roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10)) {
    if(length(x) != 1) stop("'x' must be of length 1")
    10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}
# Data loading
# WD2="/media/merly/Baymax/research/cph-visit/coverage-analysis/"

WD2="/home/merly/"
WD=paste0(WD2,"anna/")
IMG=paste0(WD,"img3/")
LOG=paste0(WD,"files/log.txt")
# offtarget=paste0(WD2,"files2/offtarget.bed")
offtargetown=paste0(WD2,"files/offtarget.anna.bed")
targetsBedFile=paste0(WD2,"files/targets.anna.3.bed")
GFF=paste0(WD2,"files/Calypte_anna.gene.CDS.2750.3.gff")

# ------------------------------------------------------------------------------
samplesFilename=paste0(WD,"files/samples.txt")
sampleNames=read.table(samplesFilename, stringsAsFactors = F)
sampleNames=unlist(sampleNames)
################################################################################
# General information in the log.
write("Dataset: anna",LOG)
# ------------------------------------------------------------------------------
write(paste("Targets file:",targetsBedFile),LOG, append=T)
targets=read.table(targetsBedFile, stringsAsFactors=F)
targets$V4=targets$V3-targets$V2
colnames(targets)=c("Scaffold","Start","End","Size")
rownames(targets)=paste0(targets$Scaffold,rep("-", nrow(targets)),targets$Start,rep("-", nrow(targets)),targets$End)
numTargets=nrow(targets)
totalBases=sum(targets$Size)
scaffoldsSplit=split(targets, targets$Scaffold)
totalScaffolds=length(scaffoldsSplit)
totalTargetPerScaffold=sapply(scaffoldsSplit,nrow)
write(paste("\tNumber of scaffolds:",prettyNum(totalScaffolds, big.mark=",")),LOG, append=T)
write(paste("\tNumber of targets:",prettyNum(nrow(targets), big.mark=",")),LOG, append=T)
write(paste("\tGenome covered by targets (bp):",prettyNum(totalBases, big.mark=",")),LOG, append=T)
write("\tSummary targets size:",LOG, append=T)
summaryTargetsSizeTable=t(summary(targets$Size))
summaryTargetsSizeTable[1,]=prettyNum(summaryTargetsSizeTable[1,],big.mark=",")
write(kable(summaryTargetsSizeTable,format="markdown"),LOG, append=T)

################################################################################
gffdata=read.table(GFF,
  colClasses=c("character","character","character","numeric",
    "numeric","character","character","character","character"
  ))
nGenesSwift=length(unique(gffdata$V9))
genes=data.frame(genes=unique(gffdata$V9))
write(paste("\tNumber of genes:",prettyNum(nGenesSwift, big.mark=",")),LOG, append=T)
# ------------------------------------------------------------------------------
# num targets per gene
sgSplit=split(gffdata, gffdata$V9)
totalTargetPerGenesA=sapply(sgSplit,nrow)
genes$sizes=sapply(sgSplit,function(x){
  sum(x$V5-x$V4)
})
png(paste0(IMG,"1.data.description.1.numTargetsPerGene.png"),1366,768)
plot(totalTargetPerGenesA,pch=16, col="darkred", ,axes=F, xlab="Genes", ylab="Number of targeted regions", main="Anna")
axis(2); axis(1)
dev.off()
# Data description
## Total  number of target regions per scaffold
png(paste0(IMG,"1.data.description.1.numTargetsPerScaffold.png"),1366,768)
plot(totalTargetPerScaffold,pch=16, col="darkred", axes=F, xlab="Scaffolds", ylab="Number of targeted regions")
axis(2, las=2); axis(1, at=seq(0,roundUpNice(totalScaffolds),50), las=2)
dev.off()
# Size distribution of the target regions
png(paste0(IMG,"1.data.description.2.size.distribution.png"), 1366,768)
h=hist(targets$Size,
       breaks =1000,
       axes=F,
       col="darkgreen",border="white",
       main="Targeted regions size (zoom x-axis)",
       xlab="Size (bp)",
       xlim=c(0,1000), prob=T
)
d=density(targets$Size)
points(d,lwd=2, type="l")
axis(1); axis(2)
dev.off()

################################################################################
# Reading data
files2 <- list.files(
  path = paste0(WD,"bedtools2/hist/"),
  pattern="*.gz$")
files2=paste0(WD,"bedtools2/hist/",files2)
cov <- list(); cov_cumul <- list()

maxcov=list()
for (i in 1:length(files2)) {
  cov[[i]] <- read.table(gzfile(files2[i]), colClasses=c("character","numeric","numeric","numeric","numeric"))
  cov_cumul[[i]] <- 1-cumsum(cov[[i]][,5])
  maxcov[[i]]=nrow(cov[[i]])
}
maxX=max(unlist(maxcov))
# Breadth vs. coverage - targeted
png(paste0(IMG,"2.breadth.vs.coverage.1.all.png"),1400,900)
plot(cov[[1]][,2], cov_cumul[[1]],
     col='darkred', type='l', lwd=2,
     xlab="Depth of coverage",
     ylab="Breadth of coverages (%)",
     xlim=c(0,roundUpNice(maxX)),
     ylim=c(0,1.0), las=2
)
for (j in 1:length(sampleNames)){
  points(cov_cumul[[j]],
         type='l',
         lwd=3,
         col=colors[j])
}
abline(v = c(0,10,20, 50,80,100,200,300), col = c("red",rep("gray60",8)))
abline(h = c(0.25,0.50,0.75, 0.8,0.9), col = "gray60")
axis(1, at=c(10,20,50,80), labels=c(10,20,50,80), las=2)
axis(2, at=c(0.25,0.50,0.75,0.9), labels=c(0.25,0.50,0.75,0.9), las=2)
legend("topright", legend=sampleNames, col=colors, lty=1, lwd=4,cex=0.8)
dev.off()

# Breadth vs. coverage - targeted
png(paste0(IMG,"2.breadth.vs.coverage.2.zoom.500.png"),1400,1000)
plot(cov[[1]][,2], cov_cumul[[1]],
     col='darkred', type='l', lwd=2,
     xlab="Depth of coverage",
     ylab="Breadth of coverages (%)",
     xlim=c(0,100),
     ylim=c(0,1.0),las=2
)
for (j in 1:length(sampleNames)){
  points(cov_cumul[[j]],
         type='l',
         lwd=3,
         col=colors[j])
}
abline(v = c(0,10,20, 50,80,100,200,300), col = c("red",rep("gray60",8)))
abline(h = c(0.25,0.50,0.75, 0.8,0.9), col = "gray60")
axis(1, at=c(10,20,50,80), labels=c(10,20,50,80), las=2)
axis(2, at=c(0.25,0.50,0.75,0.9), labels=c(0.25,0.50,0.75,0.9), las=2)
legend("topright", legend=sampleNames, col=colors, lty=1, lwd=4,cex=0.8)
dev.off()



filesWO <- list.files(
  path = paste0(WD,"bedtools2/nohist/"),
  pattern="*.gz$")
filesWO=paste0(WD,"bedtools2/nohist/",filesWO)
mergedMatrix=matrix(0, nrow=nrow(targets),ncol = length(sampleNames))
mergedMatrixGenes=matrix(0, nrow=nrow(genes),ncol = length(sampleNames))


for (i in 1:length(filesWO)) {
# for (i in c(1,48)) {
  print(i)
  data <- read.table(
    gzfile(filesWO[i]),
    colClasses=c(
      "character", #scaffold
      "character", # genewise
      "character", # CDS
      "numeric", # start
      "numeric", # end
      "character", # -
      "character", # -
      "numeric", # -
      "character", ## gene
      "numeric", # cov1
      "numeric", # cov2
      "numeric", # size
      "numeric" # breadth

    ))[c(1,4,5,9,10:13)]
  data$targetname=paste0(data[,1],rep("-", nrow(data)),data[,2],rep("-", nrow(data)),data[,3])
  colnames(data)=c("scaffold","start","end","gene","cov1","numbases","size","breadth", "targetname")
  s=split(data, data$targetname)
  sumCoveragePerTarget=sapply(s,function(x){sum(x$cov1)})
  targetsToFilterOut=intersect(names(s),rownames(targets))
  mergedMatrix[,i]=unlist(sumCoveragePerTarget[targetsToFilterOut])
  ss=split(data,data$gene)
  sumCoveragePerGene=sapply(ss,function(x){sum(x$cov1)})
  targetsToFilterOut=intersect(names(ss),genes$genes)
  mergedMatrixGenes[,i]=unlist(sumCoveragePerGene[targetsToFilterOut])
  rm(data);rm(s);rm(sumCoveragePerTarget);rm(sumCoveragePerGene)
}
rownames(mergedMatrix)=rownames(targets)
colnames(mergedMatrix)=sampleNames
rownames(mergedMatrixGenes)=rownames(genes$genes)
colnames(mergedMatrixGenes)=sampleNames
coveragePerTarget=rowMeans(mergedMatrix/targets$Size) # These gives me nTargets elems / their sizes
coveragePerGene=rowMeans(mergedMatrixGenes/genes$sizes) # These gives me nTargets elems / their sizes
coveragePerSample=colSums(mergedMatrix)/totalBases # These gives me nSamples elems/totalBases

write.table(mergedMatrix,file=paste0(WD,"files/coverage.matrix.per.target.txt"), col.names=T,row.names=T)
write.table(mergedMatrixGenes,file=paste0(WD,"files/coverage.matrix.per.gene.txt"), col.names=T,row.names=T)
write.table(coveragePerGene,file=paste0(WD,"files/coverage.per.gene.txt"))
write.table(coveragePerTarget,file=paste0(WD,"files/coverage.per.target.txt"))
write.table(coveragePerSample,file=paste0(WD,"files/coverage.per.sample.txt"))

# mergedMatrix=read.table(paste0(WD,"files2/coverage.matrix.per.target.txt")
# mergedMatrixGenes=read.table(paste0(WD,"files2/coverage.matrix.per.gene.txt")
# coveragePerTarget=unlist(read.table(paste0(WD,"files2/coverage.per.target.txt", stringsAsFactors=F, colClasses=c("character","numeric"), header=T))
# coveragePerSample=unlist(read.table(paste0(WD,"files2/coverage.per.sample.txt", stringsAsFactors=F, colClasses=c("character","numeric"), header=T))
# coveragePerGene=unlist(read.table(paste0(WD,"files2/coverage.per.gene.txt", stringsAsFactors=F, colClasses=c("character","numeric"), header=T))

png(paste0(IMG,"3.coverage.1.genes.all.png"),800,400)
par(mar=c(5,4,3,2))
layout(matrix(c(1,1,2),1,3,byrow = T))
plot(coveragePerGene,
     pch=20,
     xlab="Genes",
     ylab="Depth of coverage", axes=F)
axis(1);axis(2)
boxplot(coveragePerGene, axes=F, width=3); axis(4)
dev.off()


capturedGene=mergedMatrixGenes/genes$sizes
percetageGenesNotCaputedFully=(sum(capturedGene<1)/(dim(capturedGene)[1]*dim(capturedGene)[2]))*100
capturedTarget=mergedMatrix/targets$Size
write(paste("Percentage of gene not fully fully captured:",prettyNum(percetageGenesNotCaputedFully, big.mark=",")),LOG, append=T)
percetageTargetNotCaputedFully=(sum(capturedTarget<1)/(dim(capturedTarget)[1]*dim(capturedTarget)[2]))*100
write(paste("Percentage of targets not fully fully captured:",prettyNum(percetageTargetNotCaputedFully, big.mark=",")),LOG, append=T)

###############################################################################
par(mar=c(5,4,3,2))
layout(matrix(c(1,1,2),1,3,byrow = T))
plot(coveragePerTarget,
     pch=20,
     xlab="Targets",
     ylab="Depth of coverage", axes=F)
axis(1);axis(2)
boxplot(coveragePerTarget, axes=F, width=3); axis(4)
###############################################################################
qnt <- quantile(coveragePerTarget, probs=c(.25, .75))
H <- 1.5 * IQR(coveragePerTarget)
y <- rep("black",length(coveragePerTarget))
y[coveragePerTarget < (qnt[1] - H)] <- "red"
y[coveragePerTarget > (qnt[2] + H)] <- "red"
# coverage target regions
png(paste0(IMG,"3.coverage.1.targets.all.png"),800,400)
par(mar=c(5,4,3,2))
layout(matrix(c(1,1,2),1,3,byrow = T))
plot(coveragePerTarget,
     pch=20, col=y,
     xlab="Targeted Regions",
     ylab="Depth of coverage", axes=F)
axis(1);axis(2)
boxplot(coveragePerTarget, axes=F, width=3); axis(4)
dev.off()
###############################################################################
png(paste0(IMG,"3.coverage.2.targets.filtered.outliers.png"), 1366,768)
ydata=coveragePerTarget
ydata[coveragePerTarget < (qnt[1] - H)] <- NA
ydata[coveragePerTarget > (qnt[2] + H)] <- NA
plot(ydata,
     pch=20,
     xlab="Targeted Regions",
     ylab="Depth of coverage", axes=F)
axis(1);axis(2)
abline(h=qnt[1], col="red", lty="dashed", lwd=2)
abline(h=qnt[2], col="red", lty="dashed", lwd=2)
abline(h=mean(ydata, na.rm = T), col="red", lwd=2)
text(labels = "mean",x = 24500,y = mean(ydata, na.rm = T)+5)
text(labels = "Q.25",x = 24500,y = qnt[1]+5)
text(labels = "Q.75",x = 24500,y = qnt[2]+5)
dev.off()
###############################################################################
png(paste0(IMG,"3.coverage.3.gene.depth.vs.size.dual.png"), 900,600)
layout(matrix(c(1,2),1,2,byrow = T))
y <- rep("black",length(sort(coveragePerGene)))
y[genes$sizes < 1 ] <- "red"
y[coveragePerGene < 1 ] <- "red"
plot(coveragePerGene ,genes$sizes, pch=20, col=y,
     xlab="Depth of coverage", ylab="Gene size",axes=F)
axis(1)
axis(2)
plot(coveragePerGene ,genes$sizes, pch=20, col=y,
     xlim=c(0,300),ylim=c(0,5000),
     xlab="Depth of coverage", ylab="Gene size", main="Zoom (0-300,0-5000)",
     axes=F)
axis(2);axis(1,at=seq(0, max(coveragePerGene ),100),las=2)
abline(v=seq(0, max(coveragePerGene),100), col="gray69")
dev.off()

png(paste0(IMG,"3.coverage.3.gene.depth.vs.size.nozoom.png"), 900,600)
y <- rep("black",length(sort(coveragePerGene)))
y[genes$sizes < 1 ] <- "red"
y[coveragePerGene < 1 ] <- "red"
plot(coveragePerGene ,genes$sizes, pch=20, col=y,
     xlab="Depth of coverage", ylab="Gene size",axes=F)
axis(1);axis(2)
dev.off()


cor.test(coveragePerGene,genes$sizes)
#
# 	Pearson's product-moment correlation
#
# data:  coveragePerGene and genes$sizes
# t = -1.8991, df = 2748, p-value = 0.05766
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.073481886  0.001176326
# sample estimates:
#         cor
# -0.03620329
#
#
# > lm(coveragePerGene~genes$sizes)
#
# Call:
# lm(formula = coveragePerGene ~ genes$sizes)
#
# Coefficients:
# (Intercept)  genes$sizes
#  103.925329    -0.002982
# summary(lmd)

# Call:
# lm(formula = coveragePerGene ~ genes$sizes)
#
# Residuals:
#    Min     1Q Median     3Q    Max
#  -95.9  -36.1   -6.5   27.6 4086.8
#
# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)
# (Intercept) 103.925329   2.850115  36.464   <2e-16 ***
# genes$sizes  -0.002982   0.001570  -1.899   0.0577 .
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 91.94 on 2748 degrees of freedom
# Multiple R-squared:  0.001311,	Adjusted R-squared:  0.0009473
# F-statistic: 3.606 on 1 and 2748 DF,  p-value: 0.05766

names(coveragePerGene)=genes$genes
lmd=lm(coveragePerGene~genes$sizes)
png(paste0(IMG,"3.coverage.3.gene.depth.vs.size.lm.correlation.png"), 1000,1000)
par(mfrow=c(2,2))
plot(lmd)
dev.off()

# The graphs on the first columns look at variance homogeneity among other things,
# normally you should see no pattern in the dots but just a random clouds of points.
# In this example this is clearly not the case since we see that the spreads
# of dots increase with higher values of cyl, our homogeneity assumptions is
# violated we can go back at the beginning and build new models this one cannot
# be interpreted… Sorry m1 you looked so great…
# For the record the graph on the top right check the normality assumptions, if
# your data are normally distributed the point should fall (more or less) in a
# straight line, in this case the data are normal.
# The final graph show how each y influence the model, each points is removed at
# a time and the new model is compared to the one with the point, if the point is
# very influential then it will have a high leverage value. Points with too high
# leverage value should be removed from the dataset to remove their outlying effect
# on the model.
###############################################################################
png(paste0(IMG,"3.coverage.3.target.depth.vs.size.png"), 900,600)
layout(matrix(c(1,2),1,2,byrow = T))
y <- rep("black",length(sort(coveragePerTarget)))
y[targets$Size < 1 ] <- "red"
y[coveragePerTarget < 1 ] <- "red"
plot(coveragePerTarget,targets$Size, pch=20, col=y,
     xlab="Depth of coverage", ylab="Target size",axes=F)
axis(1)
axis(2)
plot(coveragePerTarget,targets$Size, pch=20, col=y,
     xlim=c(0,1000),ylim=c(0,5000),
     xlab="Depth of coverage", ylab="Target size", main="Zoom (0-1000,0-5000)",
     axes=F)
axis(2);axis(1,at=seq(0, max(coveragePerTarget),100),las=2)
abline(v=seq(0, max(coveragePerTarget),100), col="gray69")
dev.off()

lmd=lm(coveragePerTarget~targets$Size)
png(paste0(IMG,"3.coverage.3.target.depth.vs.size.lm.correlation.png"), 1000,1000)
par(mfrow=c(2,2))
plot(lmd)
dev.off()
###############################################################################
## Frequency of coverage per target (0  < coverage < 1000)
png(paste0(IMG,"3.coverage.4.freq.coverage.target.png"),800,600)
layout(matrix(c(1),1,1,byrow = T))
hist(coveragePerTarget,breaks=5000,
     border="white",col="darkgreen",xlim=c(0,1000),
     xlab="Depth of coverage", main="Target region coverage frequency - zoom x-axis 0-1,000")
dev.off()
###############################################################################
## Frequency of coverage per gene (0  < coverage < 1000)
png(paste0(IMG,"3.coverage.4.freq.coverage.gene.png"),800,600)
layout(matrix(c(1),1,1,byrow = T))
hist(coveragePerGene,breaks=100,
    border="white",col="darkgreen",
    xlab="Depth of coverage", main="Genes coverage frequency")
dev.off()
###############################################################################
## Coverage per sample
names(coveragePerSample)=sampleNames
qnt <- quantile(coveragePerSample, probs=c(.25, .75))
H <- 1.5 * IQR(coveragePerSample)
y <- rep("black",length(coveragePerSample))
y[coveragePerSample < (qnt[1] - H)] <- "red"
y[coveragePerSample > (qnt[2] + H)] <- "red"
# samples close to ingroup ref: h1, h5,22,30
closeToRefIngroup=c("H1","H5","H22","H30")
closeToRefIngroupIndices=c(which(sampleNames==closeToRefIngroup[1]),
                           which(sampleNames==closeToRefIngroup[2]),
                           which(sampleNames==closeToRefIngroup[3]),
                           which(sampleNames==closeToRefIngroup[4]))

pchList=rep(20, length(sampleNames))
pchList[closeToRefIngroupIndices]=17
#-------------------------------------------------------------------------------
png(paste0(IMG,"3.coverage.4.sample.all.png"), 800,600)
par(mar=c(5,8,1,1))
plot(1:length(sampleNames),coveragePerSample,
  axes=F,
  xlab="Samples",ylab="Depth of coverage", pch=pchList, col=y)
abline(h=qnt[1], col="darkblue", lty="dashed", lwd=2)
abline(h=qnt[2], col="darkblue", lty="dashed", lwd=2)
abline(h=mean(coveragePerSample), col="darkblue", lwd=2)
axis(2,at=seq(0,300, 25),las=2,cex.axis=0.6)
axis(1, at=1:length(sampleNames), las=2,labels=sampleNames,
  cex.axis=0.8
)
text(labels = "mean",x = 48,y = mean(coveragePerSample)+5, cex=0.5)
text(labels = "Q.25",x = 48,y = qnt[1]+5, cex=0.5)
text(labels = "Q.75",x = 48,y = qnt[2]+5, cex=0.5)
legend("topright", legend = c("Outlier",  "Outlier close to ingroup ref", "Not-outlier","Not-outlier close to ingroup ref"),
       col=c(rep("red",2),rep("black",2)),
       pch=c(20,17,20,17))
dev.off()
###############################################################################
png(paste0(IMG,"3.coverage.5.overview.targets.png"),1000,1000)
heatmap(as.matrix(mergedMatrix),Rowv = NA,Colv=NA, labRow = NA, labCol = sampleNames)
dev.off()
###############################################################################
png(paste0(IMG,"3.coverage.5.overview.genes.png"),1000,1000)
heatmap(as.matrix(mergedMatrixGenes),Rowv = NA,Colv=NA, labRow = NA, labCol = sampleNames)
dev.off()

png(paste0(IMG,"3.coverage.5.overview.targets.small.dataset.png"),1000,1000)
heatmap(as.matrix(mergedMatrix[1:500,]),Rowv = NA,Colv=NA, labRow = NA, labCol = sampleNames)
dev.off()
###############################################################################
png(paste0(IMG,"3.coverage.5.overview.genes.small.dataset.png"),1000,1000)
heatmap(as.matrix(mergedMatrixGenes[1:500,]),Rowv = NA,Colv=NA, labRow = NA, labCol = sampleNames)
dev.off()
###############################################################################
captured=mergedMatrix
captured[captured>0]=1
notcaptured=1-captured
nFullyCaptured=(colSums(captured)/nrow(targets)) # At leat 1 base was recovered from that target
################################################################################
png(paste0(IMG,"4.capture.1.targets.fully.captured.png"), 800,600)
plot(nFullyCaptured,
  pch=20, col="darkred",
  ylim=c(0,1),
  axes=F, type = "b",lwd="3",
  ylab="Percentage of targeted regions fully captured", xlab="Samples" )
axis(2, at=seq(0,1.0,0.1), las=2)
axis(1, at=1:length(nFullyCaptured),labels=sampleNames,las=2) # At leat 1 base was recovered from that target
dev.off()
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# OFF TARGETS
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
offtarget=read.table(offtargetown,stringsAsFactors = F, colClasses = c("character","numeric","numeric"))
offtarget$V4=offtarget$V3-offtarget$V2
colnames(offtarget)=c("Scaffold","Start","End","Size")
totalOffTargetSize=sum(offtarget$Size)
prettyNum(totalOffTargetSize, big.mark = ",")
rownames(offtarget)=paste0(offtarget$Scaffold,rep("-", nrow(offtarget)),offtarget$Start,rep("-", nrow(offtarget)),offtarget$End)
files2 <- list.files(
path = paste0(WD,"bedtools2/hist/"),
pattern="*.gz$")
files2=paste0(WD,"bedtools2/hist2/",files2)
cov <- list(); cov_cumul <- list(); maxcov=list()
################################################################################
for (i in 1:length(files2)) {
  cov[[i]] <- read.table(files2[i], colClasses=c(
    "character","numeric","numeric","numeric","numeric"))
  cov_cumul[[i]] <- 1-cumsum(cov[[i]][,5])
  maxcov[[i]]=nrow(cov[[i]])
}
maxX=max(unlist(maxcov))
# Breadth vs. coverage - targeted
png(paste0(IMG,"2.breadth.vs.coverage.1.all.offtarget1.png"),1400,900)
plot(cov[[1]][,2], cov_cumul[[1]],
   col='darkred', type='l', lwd=2,
   xlab="Depth of coverage",
   ylab="Breadth of coverages (%)",
   xlim=c(0,roundUpNice(maxX)),
   ylim=c(0,1.0), las=2
)
for (j in 1:length(sampleNames)){
points(cov_cumul[[j]],
       type='l',
       lwd=3,
       col=colors[j])
}
abline(v = c(0,10,20, 50,80,100,200,300), col = c("red",rep("gray60",8)))
abline(h = c(0.25,0.50,0.75, 0.8,0.9), col = "gray60")
axis(1, at=c(10,20,50,80), labels=c(10,20,50,80), las=2)
axis(2, at=c(0.25,0.50,0.75,0.9), labels=c(0.25,0.50,0.75,0.9), las=2)
legend("topright", legend=sampleNames, col=colors, lty=1, lwd=4,cex=0.8)
dev.off()
################################################################################
# Breadth vs. coverage - zoom
png(paste0(IMG,"2.breadth.vs.coverage.2.zoom.500.offtarget1.png"),1400,1000)
plot(cov[[1]][,2], cov_cumul[[1]],
   col='darkred', type='l', lwd=2,
   xlab="Depth of coverage",
   ylab="Breadth of coverages (%)",
   xlim=c(0,20),
   ylim=c(0,1.0),las=2
)
for (j in 1:length(sampleNames)){
points(cov_cumul[[j]],
       type='l',
       lwd=3,
       col=colors[j])
}
abline(v = c(0:10, 15), col = "gray60")
abline(h = c(0.05,0.1,0.25,0.5, 0.90), col = "gray60")
axis(1, at=c(0:4,6:9), labels=c(0:4,6:9), las=2)
axis(2, at=c(0.05,0.1,0.25,0.5, 0.90), labels=c(0.05,0.1,0.25,0.5, 0.90), las=2)
legend("topright", legend=sampleNames, col=colors, lty=1, lwd=4)
dev.off()
################################################################################
filesWO <- list.files(
  path = paste0(WD,"bedtools2/nohist2/"),
  pattern="*.gz$")
filesWO=paste0(WD,"bedtools2/nohist2/",filesWO)
mergedMatrix=matrix(0, nrow=nrow(offtarget),ncol = length(sampleNames))
################################################################################
for (i in 1:length(filesWO)) {
  # for (i in c(1,48)) {
  print(i)
  data <- read.table(
    gzfile(filesWO[i]),
    colClasses=c(
      "character", #scaffold
      "numeric", # start
      "numeric", # end
      "numeric", # -
      "numeric", # cov1
      "numeric", # cov2
      "numeric" # breadth
    ))

  data$targetname=paste0(data[,1],rep("-", nrow(data)),data[,2],rep("-", nrow(data)),data[,3])
  colnames(data)=c("scaffold","start","end","pos","cov1","size","breadth", "targetname")
  s=split(data, data$targetname)
  sumCoveragePerTarget=sapply(s,function(x){sum(x$cov1)})
  targetsToFilterOut=intersect(names(sumCoveragePerTarget),rownames(offtarget))
  mergedMatrix[,i]=unlist(sumCoveragePerTarget[targetsToFilterOut])
  rm(data);rm(s);rm(sumCoveragePerTarget)
}
################################################################################
totalBases=sum(offtarget$Size)
rownames(mergedMatrix)=rownames(offtarget)
colnames(mergedMatrix)=sampleNames
coveragePerTarget=rowMeans(mergedMatrix/offtarget$Size) # These gives me nTargets elems / their sizes
coveragePerSample=colSums(mergedMatrix)/totalBases # These gives me nSamples elems/totalBases
################################################################################
write.table(mergedMatrix,file=paste0(WD,"files2/offtarget1.coverage.matrix.per.target.txt"), col.names=T, row.names=T)
write.table(coveragePerTarget,file=paste0(WD,"files2/offtarget1.coverage.per.target.txt"), col.names=T, row.names=T)
write.table(coveragePerSample,file=paste0(WD,"files2/offtarget1.coverage.per.sample.txt"), col.names=T, row.names=T)
################################################################################
mergedMatrix=read.table(paste0(WD,"files2/offtarget1.coverage.matrix.per.target.txt"))
coveragePerTarget=read.table(paste0(WD,"files2/offtarget1.coverage.per.target.txt"), stringsAsFactors=F, colClasses=c("character","numeric"), header=T)
coveragePerSample=read.table(paste0(WD,"files2/offtarget1.coverage.per.sample.txt"), stringsAsFactors=F, colClasses=c("character","numeric"), header=T)
coveragePerTarget=unlist(coveragePerTarget)
coveragePerSample=unlist(coveragePerSample)
################################################################################
capturedTarget=mergedMatrix/offtarget$Size
percetageTargetNotCaputedFully=(sum(capturedTarget<1)/(dim(capturedTarget)[1]*dim(capturedTarget)[2]))*100
write(paste("OFFTarget1-Percentage of targets not fully fully captured:",prettyNum(percetageTargetNotCaputedFully, big.mark=",")),LOG, append=T)
################################################################################
qnt <- quantile(coveragePerTarget, probs=c(.25, .75))
H <- 1.5 * IQR(coveragePerTarget)
y <- rep("black",length(coveragePerTarget))
y[coveragePerTarget < (qnt[1] - H)] <- "red"
y[coveragePerTarget > (qnt[2] + H)] <- "red"
# coverage target regions
png(paste0(IMG,"3.coverage.1.targets.all.offtarget1.png"),800,400)
par(mar=c(5,4,3,2))
layout(matrix(c(1,1,2),1,3,byrow = T))
plot(coveragePerTarget,
   pch=20, col=y,
   xlab="Targeted Regions",
   ylab="Depth of coverage", axes=F)
axis(1);axis(2)
boxplot(coveragePerTarget, axes=F, width=3); axis(4)
dev.off()
###############################################################################
png(paste0(IMG,"3.coverage.2.targets.filtered.outliers.offtarget1.png"), 1366,768)
ydata=coveragePerTarget
ydata[coveragePerTarget < (qnt[1] - H)] <- NA
ydata[coveragePerTarget > (qnt[2] + H)] <- NA
plot(ydata,
   pch=20,
   xlab="Targeted Regions",
   ylab="Depth of coverage", axes=F)
axis(1);axis(2)
abline(h=qnt[1], col="red", lty="dashed", lwd=2)
abline(h=qnt[2], col="red", lty="dashed", lwd=2)
abline(h=mean(ydata, na.rm = T), col="red", lwd=2)
text(labels = "mean",x = 24500,y = mean(ydata, na.rm = T)+5)
text(labels = "Q.25",x = 24500,y = qnt[1]+5)
text(labels = "Q.75",x = 24500,y = qnt[2]+5)
dev.off()
###############################################################################
png(paste0(IMG,"3.coverage.3.target.depth.vs.size.offtarget1.png"), 900,600)
layout(matrix(c(1,2),1,2,byrow = T))
y <- rep("black",length(sort(coveragePerTarget)))
y[offtarget$Size < 1 ] <- "red"
y[coveragePerTarget < 1 ] <- "red"
plot(coveragePerTarget,offtarget$Size, pch=20, col=y,
   xlab="Depth of coverage", ylab="Target size",axes=F)
axis(1)
axis(2)
plot(coveragePerTarget,offtarget$Size, pch=20, col=y,
   xlim=c(0,1000),ylim=c(0,5000),
   xlab="Depth of coverage", ylab="Target size", main="Zoom (0-1000,0-5000)",
   axes=F)
axis(2);axis(1,at=seq(0, max(coveragePerTarget),100),las=2)
abline(v=seq(0, max(coveragePerTarget),100), col="gray69")
dev.off()
################################################################################
## Frequency of coverage per target (0  < coverage < 1000)
png(paste0(IMG,"3.coverage.4.freq.coverage.target.offtarget1.png"),800,600)
layout(matrix(c(1),1,1,byrow = T))
hist(coveragePerTarget,breaks=50000,
   border="white",col="darkgreen",xlim=c(0,300),
   xlab="Depth of coverage", main="Target region coverage frequency - zoom x-axis 0-300")
dev.off()
###############################################################################
## Coverage per sample
names(coveragePerSample)=sampleNames
qnt <- quantile(coveragePerSample, probs=c(.25, .75))
H <- 1.5 * IQR(coveragePerSample)
y <- rep("black",length(coveragePerSample))
y[coveragePerSample < (qnt[1] - H)] <- "red"
y[coveragePerSample > (qnt[2] + H)] <- "red"
# samples close to ingroup ref: h1, h5,22,30
closeToRefIngroup=c("H1","H5","H22","H30")
closeToRefIngroupIndices=c(which(sampleNames==closeToRefIngroup[1]),
                         which(sampleNames==closeToRefIngroup[2]),
                         which(sampleNames==closeToRefIngroup[3]),
                         which(sampleNames==closeToRefIngroup[4]))

pchList=rep(20, length(sampleNames))
pchList[closeToRefIngroupIndices]=17
################################################################################
png(paste0(IMG,"3.coverage.4.sample.all.offtarget1.png"), 800,600)
par(mar=c(5,8,1,1))
plot(1:length(sampleNames),coveragePerSample,
  axes=F,
  xlab="Samples",ylab="Depth of coverage", pch=pchList, col=y)
abline(h=qnt[1], col="darkblue", lty="dashed", lwd=2)
abline(h=qnt[2], col="darkblue", lty="dashed", lwd=2)
abline(h=mean(coveragePerSample), col="darkblue", lwd=2)
axis(2,at=seq(0,300, 25),las=2,cex.axis=0.6)
axis(1, at=1:length(sampleNames), las=2,labels=sampleNames,
cex.axis=0.8
)
text(labels = "mean",x = 48,y = mean(coveragePerSample)+5, cex=0.5)
text(labels = "Q.25",x = 48,y = qnt[1]+5, cex=0.5)
text(labels = "Q.75",x = 48,y = qnt[2]+5, cex=0.5)
legend("topright", legend = c("Outlier",  "Outlier close to ingroup ref", "Not-outlier","Not-outlier close to ingroup ref"),
     col=c(rep("red",2),rep("black",2)),
     pch=c(20,17,20,17))
dev.off()
################################################################################
png(paste0(IMG,"3.coverage.5.overview.targets.offtarget1.png"),1000,1000)
heatmap(as.matrix(mergedMatrix),Rowv = NA,Colv=NA, labRow = NA, labCol = sampleNames)
dev.off()



################################################################################
# ANGSD
#
# WD2="/media/merly/Baymax/research/cph-visit/coverage-analysis/"
# WD=paste0(WD2,"anna/")
# sampleLabelsSwift=unlist(read.table(paste0(WD,"files/samples.txt"), stringsAsFactors = F))
# angsdAllSampleFileS=paste0(WD,"files/anna.originals.angsd..depthSample")
# angsdAllGlobalFileS=paste0(WD,"files/anna.originals.angsd..depthGlobal")
# angsdTargetedSampleFileS=paste0(WD,"files/anna.ontarget.angsd..depthSample")
# angsdTargetedGlobalFileS=paste0(WD,"files/anna.ontarget.angsd..depthGlobal")
# angsdUntargetedSampleFileS=paste0(WD,"files/anna.offtarget2.angsd..depthSample")
# angsdUntargetedGlobalFileS=paste0(WD,"files/anna.offtarget2.angsd..depthGlobal")
#
# allSampleS=data.frame((t(read.table(angsdAllSampleFileS))))
# allGlobalS=(scan(angsdAllGlobalFileS))
# rownames(allSampleS)=paste0(1:nrow(allSampleS))
# colnames(allSampleS)=sampleLabelsSwift
# targetedSampleS=data.frame((t(read.table(angsdTargetedSampleFileS))))
# targetedGlobalS=(scan(angsdTargetedGlobalFileS))
# rownames(targetedSampleS)=paste0(1:nrow(targetedSampleS))
# colnames(targetedSampleS)=sampleLabelsSwift
# untargetedSampleS=data.frame((t(read.table(angsdUntargetedSampleFileS))))
# untargetedGlobalS=(scan(angsdUntargetedGlobalFileS))
# rownames(untargetedSampleS)=paste0(1:nrow(untargetedSampleS))
# colnames(untargetedSampleS)=sampleLabelsSwift
# ################################################################################
# ## Depth of coverage per sample (all)
# png(paste0(IMG,"5.angsd.depth.all.png"), 1366,1366)
# heatmap(as.matrix(allSampleS), Rowv = NA,Colv=NA, labCol = sampleLabelsSwift)
# dev.off()
# png(paste0(IMG,"5.angsd.depth.targeted.png"), 1366,1366)
# heatmap(as.matrix(targetedSampleS), Rowv = NA,Colv=NA, labCol = sampleLabelsSwift)
# dev.off()
# png(paste0(IMG,"5.angsd.depth.untargeted.png"), 1366,1366)
# heatmap(as.matrix(untargetedSampleS),Rowv = NA,Colv=NA, labCol = sampleLabelsSwift)
# dev.off()
