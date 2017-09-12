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











nTargetedRegions=read.table("/home/merly/git/cph-visit/coverage-analysis/empirical/files/numTargets.per.sample.txt")
# nTargetedRegions
plot(1:nrow(nTargetedRegions), ((nTargetedRegions[,2]/numTargets))*100,axes=F, ylab="Percentage of targeted regions fully captured", xlab="Samples", pch=16, col="darkred", ylim=c(95,100))
axis(2, at=seq(95,100), las=2)
axis(1, at=seq(1,nrow(nTargetedRegions)),labels=nTargetedRegions[,1],las=2)


zeroCover=lapply(cov,function(x){(x[1,3]/x[1,4])*100})
zeroCoverageDF=data.frame(samples=labs, zerocoverage=zeroCover)
dlabs=c()
for (item in labs){dlabs=c(dlabs,item)}
plot(1:length(labs),zeroCover, axes=F, ylab="Percentage of not covered sites", xlab="Samples", pch=16, col="darkred")
axis(2)
axis(1, at=seq(1,length(dlabs)),labels=dlabs,las=2)

suppressMessages(suppressWarnings(dev.print(pdf, "bedtools.zerocoverage.pdf")))




######################################################################################################################################3# Reading data
filesZero <- list.files(
  path = "empirical/targetZero/",
  pattern="txt$")
ftz=paste0("empirical/targetZero/",filesZero)
targetZero <- list()
setsTargetZero <- list()
scaffolds <- list()
sizes<-list()
for (i in 1:length(ftz)) {
  targetZero[[i]] <- read.table(ftz[i],stringsAsFactors = T)
  setsTargetZero[[i]]<-paste0(targetZero[[i]][,1],"-",targetZero[[i]][,2],"-",targetZero[[i]][,3])
  scaffolds[[i]]<-paste0(targetZero[[i]][,1])
  sizes[[i]]<-(targetZero[[i]][,3]-targetZero[[i]][,2])+1
  names(sizes[[i]])<-setsTargetZero[[i]]
}

ss=lapply(setsTargetZero,length)
s=0
for(item in 1:length(ss)){
  s=s+ss[[item]]
}
allScaffolds=Reduce(union,scaffolds)
intersectedRegions=Reduce(intersect, setsTargetZero)
allRegions=Reduce(union, setsTargetZero)
allData=data.frame(matrix(0,nrow=length(labs), ncol=length(allRegions)))
rownames(allData)<-labs
colnames(allData)<-allRegions

histSizes=list()

for (i in 1:length(ftz)) {
  allData[i,intersect(colnames(allData),setsTargetZero[[i]])]=1
  histSizes[[i]]<-sizes[[i]][intersectedRegions]
}

heatmap(as.matrix(allData),Rowv = NA, Colv = NA, col=paste("gray",1:99,sep=""))


## Sizes of the targets not captured across all samples
### Overview
unlistedHist=unlist(histSizes)
hist(unlistedHist,breaks=150, xlab="Read size (bp)",col="darkgreen",border="white", main="")
abline(v=50, col="darkred")
abline(v=100, col="darkred")
abline(v=120, col="darkred")
abline(v=mean(unlistedHist), col="blue")

### Zoom
hist(unlistedHist,breaks=1000,xlab="Read size (bp)", xlim=c(0,1000),col="darkgreen",border="white", main="", axes=F)
abline(v=50, col="darkred")
abline(v=100, col="darkred")
abline(v=120, col="darkred")
abline(v=mean(unlistedHist), col="blue")
axis(2)
axis(1,at=c(0,50,100,120,round(mean(unlistedHist)),
            seq(200,1000,200)),
     labels = c(0,50,100,120,paste0("~",round(mean(unlistedHist))),
                seq(200,1000,200)),
     las=2)

### Size of targeted regions not captured (Summary)
kable(t(summary(unlistedHist)), format="markdown")

### Information in terms of numbers

df=data.frame(Intersection=length(intersectedRegions), Union=length(allRegions), scaffolds=length(allScaffolds))
colnames(df)<-c("**Number of targeted regions missing throusgh all samples**",
                "**Number of the different targeted regions missing though all the sampless**",
                "**Number of different scaffolds**")
rownames(df)=c("Value")
kable(t(df),format="markdown")


heatmap(as.matrix(avgCoveragePerTargetPerSample[,1:100]), Rowv = NA,Colv = NA)

sizes=unlist(sizes)
tableSizes=table(sizes)
tableSizes=data.frame(tableSizes)
maxFreq=max(as.numeric(tableSizes$Freq))
maxSize=max(as.numeric(as.character(tableSizes$sizes)))


## Depth per site from targeted regions
hist(sizes, breaks = 100, col="darkgreen",border="white")
hist(sizes, breaks = 300, col="darkgreen",border="white", xlim=c(0,1000))


## ANGSD Depth calculations

sampleLabels=unlist(read.table("label.samples.txt", stringsAsFactors = F))
angsdAllSampleFile="empirical/angsd-depth/all.depthSample"
angsdAllGlobalFile="empirical/angsd-depth/all.depthGlobal"
angsdTargetedSampleFile="empirical/angsd-depth/filtered.depthSample"
angsdTargetedGlobalFile="empirical/angsd-depth/filtered.depthGlobal"
angsdUntargetedSampleFile="empirical/angsd-depth/untargeted-captured.depthSample"
angsdUntargetedGlobalFile="empirical/angsd-depth/untargeted-captured.depthGlobal"

allSample=data.frame(log10(t(read.table(angsdAllSampleFile))))
allGlobal=log10(scan(angsdAllGlobalFile))
rownames(allSample)=paste0(1:nrow(allSample), "x")
colnames(allSample)=sampleLabels

targetedSample=data.frame(log10(t(read.table(angsdTargetedSampleFile))))
targetedGlobal=log10(scan(angsdTargetedGlobalFile))
rownames(targetedSample)=paste0(1:nrow(targetedSample), "x")
colnames(targetedSample)=sampleLabels

untargetedSample=data.frame(log10(t(read.table(angsdUntargetedSampleFile))))
untargetedGlobal=log10(scan(angsdUntargetedGlobalFile))
rownames(untargetedSample)=paste0(1:nrow(untargetedSample), "x")
colnames(untargetedSample)=sampleLabels

## Depth of coverage per sample (all)

heatmap.2(as.matrix(allSample), Rowv=NULL,Colv=NULL,
          key.title = "",
          key.xlab="number of sites (log10)",
          key.ylab = "Frequency",
          scale="none",
          margins=c(3,1), # ("margin.Y", "margin.X")
          trace='none',
          symkey=FALSE,
          symbreaks=FALSE,
          dendrogram='none',
          density.info='histogram',
          denscol="black",
          keysize=0.5,
          #( "bottom.margin", "left.margin", "top.margin", "left.margin" )
          key.par=list(mar=c(5,0,3,0)),
          # lmat -- added 2 lattice sections (5 and 6) for padding
          lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(2.5, 5), lwid=c(1, 10, 1)
)

## Depth of coverage per sample (Targeted)

heatmap.2(as.matrix(targetedSample), Rowv=NULL,Colv=NULL,
          key.title = "",
          key.xlab="number of sites (log10)",
          key.ylab = "Frequency",
          scale="none",
          margins=c(3,1), # ("margin.Y", "margin.X")
          trace='none',
          symkey=FALSE,
          symbreaks=FALSE,
          dendrogram='none',
          density.info='histogram',
          denscol="black",
          keysize=0.5,
          #( "bottom.margin", "left.margin", "top.margin", "left.margin" )
          key.par=list(mar=c(5,0,3,0)),
          # lmat -- added 2 lattice sections (5 and 6) for padding
          lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(2.5, 5), lwid=c(1, 10, 1))

## Depth of coverage per sample (Untargeted)

heatmap.2(as.matrix(untargetedSample),Rowv=NULL,Colv=NULL,
          key.title = "",
          key.xlab="number of sites (log10)",
          key.ylab = "Frequency",
          scale="none",
          margins=c(3,1), # ("margin.Y", "margin.X")
          trace='none',
          symkey=FALSE,
          symbreaks=FALSE,
          dendrogram='none',
          density.info='histogram',
          denscol="black",
          keysize=0.5,
          #( "bottom.margin", "left.margin", "top.margin", "left.margin" )
          key.par=list(mar=c(5,0,3,0)),
          # lmat -- added 2 lattice sections (5 and 6) for padding
          lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(2.5, 5), lwid=c(1, 10, 1)
)


#### General coverage levels (Log-scale)
ymaxlim=10
xmaxlim=max(length(allGlobal),length(targetedGlobal),length(untargetedGlobal))
selectedColors=c("darkred","darkgreen","darkblue")
plot(
  allGlobal,
  type="l",
  lwd=2,
  xlab="Depth",ylab="Frequency (Log10)",
  xlim=c(0,xmaxlim),
  ylim=c(0,ymaxlim),
  col=selectedColors[1]
)
points(targetedGlobal,	type="l",lwd=2,col=selectedColors[2])
points(untargetedGlobal, type="l", lwd=2, col=selectedColors[3])
dataTypeLabels=c("All", "Targeted", "Untargeted")
legend("bottomright", legend=dataTypeLabels, col=selectedColors, lty=1, lwd=4)


#### Coverage level per sample
par(mfrow=c(3,1))
plot(allSample[,1], type="l",col="darkred",lwd=2,xlab="Depth", ylab="Frequency (Log10)",ylim=c(2,9), main="All")
for (i in 2:ncol(allSample)){
  points(allSample[, i],type='l',lwd=3,col=colors[i])
}

plot(targetedSample[,1], type="l",col="darkred",lwd=2,xlab="Depth", ylab="Frequency (Log10)",ylim=c(2,9), main="Targeted")
for (i in 2:ncol(targetedSample)){
  points(targetedSample[, i],type='l',lwd=3,col=colors[i])
}

plot(untargetedSample[,1], type="l",col="darkred",lwd=2,xlab="Depth", ylab="Frequency (Log10)",ylim=c(2,9), main="Untargeted")
for (i in 2:ncol(untargetedSample)){
  points(untargetedSample[, i],type='l',lwd=3,col=colors[i])
}
