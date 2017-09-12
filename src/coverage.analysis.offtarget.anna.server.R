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
WD2="/media/merly/Baymax/research/cph-visit/coverage-analysis/"

# WD2="/home/merly/"
WD=paste0(WD2,"anna/")
LOG=paste0(WD,"files/log.txt")
offtargetown=paste0(WD2,"files/offtarget.anna.bed")
GFF=paste0(WD2,"files/Calypte_anna.gene.CDS.2750.3.gff")

# ------------------------------------------------------------------------------
samplesFilename=paste0(WD,"files/samples.txt")
sampleNames=read.table(samplesFilename, stringsAsFactors = F)
sampleNames=unlist(sampleNames)
################################################################################
# ------------------------------------------------------------------------------
targets=read.table(targetsBedFile, stringsAsFactors=F)
targets$V4=targets$V3-targets$V2
colnames(targets)=c("Scaffold","Start","End","Size")
rownames(targets)=paste0(targets$Scaffold,rep("-", nrow(targets)),targets$Start,rep("-", nrow(targets)),targets$End)
numTargets=nrow(targets)
totalBases=sum(targets$Size)
scaffoldsSplit=split(targets, targets$Scaffold)
totalScaffolds=length(scaffoldsSplit)
totalTargetPerScaffold=sapply(scaffoldsSplit,nrow)
summaryTargetsSizeTable=t(summary(targets$Size))
summaryTargetsSizeTable[1,]=prettyNum(summaryTargetsSizeTable[1,],big.mark=",")

################################################################################
gffdata=read.table(GFF,
  colClasses=c("character","character","character","numeric",
    "numeric","character","character","character","character"
  ))
nGenesSwift=length(unique(gffdata$V9))
genes=data.frame(genes=unique(gffdata$V9))
# ------------------------------------------------------------------------------
# num targets per gene
sgSplit=split(gffdata, gffdata$V9)
totalTargetPerGenesA=sapply(sgSplit,nrow)
genes$sizes=sapply(sgSplit,function(x){
  sum(x$V5-x$V4)
})
offtarget=read.table(offtargetown,stringsAsFactors = F, colClasses = c("character","numeric","numeric"))
offtarget$V4=offtarget$V3-offtarget$V2
colnames(offtarget)=c("Scaffold","Start","End","Size")
totalOffTargetSize=sum(offtarget$Size)
prettyNum(totalOffTargetSize, big.mark = ",")
rownames(offtarget)=paste0(offtarget$Scaffold,rep("-", nrow(offtarget)),offtarget$Start,rep("-", nrow(offtarget)),offtarget$End)
files2 <- list.files(
path = paste0(WD,"bedtools2/hist/"),
pattern="*.gz$")
files2=paste0(WD,"bedtools2/hist/",files2)
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
# png(paste0(IMG,"2.breadth.vs.coverage.1.all.offtarget1.png"),1400,900)
png(paste0("2.breadth.vs.coverage.1.all.offtarget1.png"),1400,900)
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
png(paste0("2.breadth.vs.coverage.2.zoom.500.offtarget1.png"),1400,1000)
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
  path = paste0(WD,"bedtools2/nohist/"),
  pattern="*.gz$")
filesWO=paste0(WD,"bedtools2/nohist/",filesWO)
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
  colnames(data)=c("scaffold","start","end","cov1","numbases","size","breadth", "targetname")
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
################################################################################
write.table(mergedMatrix,file=paste0("offtarget1.coverage.matrix.per.target.txt"), col.names=T, row.names=T)
write.table(coveragePerTarget,file=paste0("offtarget1.coverage.per.target.txt"), col.names=T, row.names=T)
write.table(coveragePerSample,file=paste0("offtarget1.coverage.per.sample.txt"), col.names=T, row.names=T)
coveragePerTarget=unlist(coveragePerTarget)
coveragePerSample=unlist(coveragePerSample)
################################################################################
capturedTarget=mergedMatrix/offtarget$Size
percetageTargetNotCaputedFully=(sum(capturedTarget<1)/(dim(capturedTarget)[1]*dim(capturedTarget)[2]))*100
################################################################################
qnt <- quantile(coveragePerTarget, probs=c(.25, .75))
H <- 1.5 * IQR(coveragePerTarget)
y <- rep("black",length(coveragePerTarget))
y[coveragePerTarget < (qnt[1] - H)] <- "red"
y[coveragePerTarget > (qnt[2] + H)] <- "red"
# coverage target regions
png(paste0("3.coverage.1.targets.all.offtarget1.png"),800,400)
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
png(paste0("3.coverage.2.targets.filtered.outliers.offtarget1.png"), 1366,768)
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
png(paste0("3.coverage.3.target.depth.vs.size.offtarget1.png"), 900,600)
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
png(paste0("3.coverage.4.freq.coverage.target.offtarget1.png"),800,600)
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
png(paste0("3.coverage.4.sample.all.offtarget1.png"), 800,600)
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
png(paste0("3.coverage.5.overview.targets.offtarget1.png"),1000,1000)
heatmap(as.matrix(mergedMatrix),Rowv = NA,Colv=NA, labRow = NA, labCol = sampleNames)
dev.off()
