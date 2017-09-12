library(ggplot2)
suppressWarnings(suppressMessages(library(gplots)))
library(RColorBrewer)
suppressWarnings(suppressMessages(library(knitr)))
colors = colorRampPalette(brewer.pal(9, "Set1"))(46)
# Data loading
WD="/media/merly/Baymax/research/cph-visit/map2Ann/"
IMG="/media/merly/Baymax/research/cph-visit/map2Ann/img/"
targetsBedFile="/media/merly/Baymax/research/cph-visit/files/targets.bed"
targets=read.table(targetsBedFile, stringsAsFactors=F)
targets$V4=targets$V3-targets$V2
colnames(targets)=c("Scaffold","Start","End","Size")
rownames(targets)=paste0(targets$Scaffold,rep("-", nrow(targets)),targets$Start,rep("-", nrow(targets)),targets$End)
numTargets=nrow(targets)
totalBases=sum(targets$Size)
scaffoldsSplit=split(targets, targets$Scaffold)
totalScaffolds=length(scaffoldsSplit)
totalTargetPerScaffold=sapply(scaffoldsSplit,nrow)
# Filtering target regions with size 0
targets=targets[targets$Size>0]

# Data description
## Total  number of target regions per scaffold
png(paste0(IMG,"1.data.description.1.numTargetsPerScaffold.png"), width=15,height=6)
plot(totalTargetPerScaffold,pch=16, col="darkred", axes=F, xlab="Scaffolds", ylab="Number of targeted regions")
axis(2); axis(1)
dev.off()
# Size distribution of the target regions
png(paste0(IMG,"1.data.description.2.size.distribution.png"), width=15,height=6)
h=hist(targets$Size,
       breaks =500,
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


# kable(t(summary(targets$Size)), format="markdown")
samplesFilename=paste0(WD,"files/samples.txt")
sampleNames=read.table(samplesFilename, stringsAsFactors = F)
sampleNames=unlist(sampleNames)


# Reading data
files <- list.files(
  path = paste0(WD,"bedtools/hist/"),
  pattern="*.gz$")
files=paste0(WD,"bedtools/hist/",files)
cov <- list(); cov_cumul <- list()

for (i in 1:length(files)) {
  cov[[i]] <- read.table(files[i], colClasses=c("character","numeric","numeric","numeric","numeric"))
  cov_cumul[[i]] <- 1-cumsum(cov[[i]][,5])
}

# Breadth vs. coverage - targeted
png(paste0(IMG,"2.breadth.vs.coverage.1.all.png"), width=20,height=10)
plot(cov[[1]][,2], cov_cumul[[1]],
     col='darkred', type='l', lwd=2,
     xlab="Depth of coverage",
     ylab="Breadth of coverages (%)",
     ylim=c(0,1.0)
)
for (j in 1:46){
  points(cov_cumul[[j]],
         type='l',
         lwd=3,
         col=colors[j])
}
abline(v = c(20, 50,80,100,200,300), col = "gray60")
abline(h = c(0.25,0.50,0.9), col = "gray60")
axis(1, at=c(20,50,80), labels=c(20,50,80), las=2)
axis(2, at=c(0.25,0.5, 0.90), labels=c(0.25,0.5,0.90), las=2)
legend("topright", legend=sampleNames, col=colors, lty=1, lwd=4)
dev.off()

# Breadth vs. coverage - targeted
png(paste0(IMG,"2.breadth.vs.coverage.2.zoom.500.png"), width=20,height=10)
plot(cov[[1]][,2], cov_cumul[[1]],
     col='darkred', type='l', lwd=2,
     xlab="Depth of coverage",
     ylab="Breadth of coverages (%)",
     xlim=c(0,500),
     ylim=c(0,1.0)
)
for (j in 1:46){
  points(cov_cumul[[j]],
         type='l',
         lwd=3,
         col=colors[j])
}
abline(v = c(20, 50,80,100,200,300), col = "gray60")
abline(h = c(0.25,0.50,0.9), col = "gray60")
axis(1, at=c(20,50,80), labels=c(20,50,80), las=2)
axis(2, at=c(0.25,0.5, 0.90), labels=c(0.25,0.5,0.90), las=2)
legend("topright", legend=sampleNames, col=colors, lty=1, lwd=4)
dev.off()



filesWO <- list.files(
  path = paste0(WD,"bedtools/nohist/"),
  pattern="*.gz$")
filesWO=paste0(WD,"bedtools/nohist/",filesWO)

covwo <- list();
for (i in 1:length(filesWO)) {
  print(i)
  covwo[[i]] <- read.table(
    filesWO[i],
    colClasses=c(
      "character",
      "character",
      "character",
      "numeric",
      "numeric",
      "factor",
      "factor",
      "numeric",
      "factor",
      "numeric",
      "numeric",
      "numeric",
      "numeric"
    ))[,c(1,4,5,10:13)]
  covwo[[i]]$targetname=paste0(covwo[[i]][,1],rep("-", nrow(covwo[[i]])),covwo[[i]][,2],rep("-", nrow(covwo[[i]])),covwo[[i]][,3])
  colnames(covwo[[i]])=c("scaffold","start","end","cov1","cov2","size","breadth", "targetname")
}
mergedMatrix=matrix(0, nrow=nrow(targets),ncol = length(sampleNames))
for (i in 1:length(files)) {
  print(i)
  s=split(covwo[[i]], covwo[[i]]$targetname)
  sumCoveragePerTarget=sapply(s,function(x){sum(x$cov1)})
  mergedMatrix[,i]=unlist(sumCoveragePerTarget)
}
rownames(mergedMatrix)=rownames(targets)
colnames(mergedMatrix)=sampleNames
coveragePerTarget=rowMeans(mergedMatrix)/targets$Size # These gives me nTargets elems / their sizes
coveragePerSample=colMeans(mergedMatrix)/totalBases # These gives me nSamples elems/totalBases

qnt <- quantile(coveragePerTarget, probs=c(.25, .75))
H <- 1.5 * IQR(coveragePerTarget)
y <- rep("black",length(coveragePerTarget))
y[coveragePerTarget < (qnt[1] - H)] <- "red"
y[coveragePerTarget > (qnt[2] + H)] <- "red"

# coverage target regions
png(paste0(IMG,"3.coverage.1.targets.all.png"), width=15, height = 10)
layout(matrix(c(1,1,2,3,3,4),2,3,byrow = T))
plot(coveragePerTarget,
     pch=20, col=y,
     main="Sorted by ID",
     xlab="Targeted Regions",
     ylab="Depth of coverage", axes=F)
axis(1);axis(2)
boxplot(coveragePerTarget, axes=F, width=3); axis(4)
ydata=sort(coveragePerTarget)
y[sort(coveragePerTarget) < (qnt[1] - H)] <- "red"
y[sort(coveragePerTarget) > (qnt[2] + H)] <- "red"
plot(ydata,
     pch=20, col=y,
     main="Sorted by coverage",
     xlab="Targeted Regions",
     ylab="Depth of coverage", axes=F)
axis(1);axis(2)
boxplot(coveragePerTarget, axes=F, width=3); axis(4)
dev.off()

png(paste0(IMG,"3.coverage.2.targets.filtered.outliers.png"), width=15)
ydata=coveragePerTarget
y[coveragePerTarget < (qnt[1] - H)] <- "red"
y[coveragePerTarget > (qnt[2] + H)] <- "red"
ydata[coveragePerTarget < (qnt[1] - H)] <- NA
ydata[coveragePerTarget > (qnt[2] + H)] <- NA
plot(ydata,
     pch=20, col=y,
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

png(paste0(IMG,"3.coverage.3.depth.vs.size.png"), width=15)
layout(matrix(c(1,2),1,2,byrow = T))
y <- rep("black",length(sort(coveragePerTarget)))
y[targets$Size < 1 ] <- "red"
y[coveragePerTarget < 1 ] <- "red"
y[coveragePerTarget > 1000] <- "red"
plot(coveragePerTarget,targets$Size, pch=20, col=y,
     xlab="Depth of coverage", ylab="Target size")
abline(v=c(100,200,300,400,500, 1000,2000))
plot(coveragePerTarget,targets$Size, pch=20, col=y,
     xlim=c(0,5000),ylim=c(0,5000),
     xlab="Depth of coverage", ylab="Target size", main="Zoom (0-5000,0-5000)")
abline(v=c(100,200,300,400,500, 1000,2000))
dev.off()

## Frequency of coverage per target (0  < coverage < 1000)
png(paste0(IMG,"3.coverage.3.depth.vs.size.png"), width=15)
hist(coveragePerTarget[coveragePerTarget>0 && coveragePerTarget<1000],breaks=10000,
     border="white",col="darkgreen",
     xlim=c(0,1000),
     xlab="Depth of coverage", main="Target region coverage frequency (0  < coverage < 1000)")
dev.off()

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


png(paste0(IMG,"3.coverage.4.sample.all"), width=15)
plot(1:46,coveragePerSample, axes=F, xlab="Samples",ylab="Depth of coverage", pch=pchList, col=y)
abline(h=qnt[1], col="darkblue", lty="dashed", lwd=2)
abline(h=qnt[2], col="darkblue", lty="dashed", lwd=2)
abline(h=mean(coveragePerSample), col="darkblue", lwd=2)
axis(2); axis(1, at=1:46, las=2,labels=sampleNames)
text(labels = "mean",x = 46,y = mean(coveragePerSample)+5)
text(labels = "Q.25",x = 46,y = qnt[1]+5)
text(labels = "Q.75",x = 46,y = qnt[2]+5)
legend("topright", legend = c("Outlier",  "Outlier close to ingroup ref", "Not-outlier","Not-outlier close to ingroup ref"),
       col=c(rep("red",2),rep("black",2)),
       pch=c(20,17,20,17))
dev.off()

png(paste0(IMG,"3.coverage.5.overview.png"), heigth=20,width=20)
heatmap(as.matrix(mergedMatrix),Rowv = NA,Colv=NA, labRow = NA, labCol = sampleNames)
dev.off()


captured=mergedMatrix
captured[captured>0]=1
notcaptured=1-captured
nFullyCaptured=(colSums(captured)/nrow(targets)) # At leat 1 base was recovered from that target

png(paste0(IMG,"4.capture.1.targets.fully.captured"), width=15)
plot(nFullyCaptured,
  pch=20, col="darkred",
  ylim=c(0.95,1),
  axes=F, type = "b",lwd="3",
  ylab="Percentage of targeted regions fully captured", xlab="Samples" )
axis(2, at=seq(0.95,1.0,0.01), las=2)
axis(1, at=1:length(nFullyCaptured),labels=sampleNames,las=2) # At leat 1 base was recovered from that target
dev.off()
