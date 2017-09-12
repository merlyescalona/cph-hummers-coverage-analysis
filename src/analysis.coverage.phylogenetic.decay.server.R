packages<-c("ape","geiger","apTreeshape","ggplot2","gplots","RColorBrewer","knitr","phangorn","futile.logger","phytools")
for(pkg in packages ){
  suppressMessages(library(pkg,character.only=TRUE,quietly=TRUE))
}
sampleColors = colorRampPalette(brewer.pal(9, "Set1"))(48)

treefile="/home/merly/git/cph-visit/coverage-analysis/capture-phylotenetic-decay/trees/RAxML_bipartitions.2011.concat"
tree=read.nexus(treefile)
endNode=23; startNode=49
nextStep=endNode
path=c(); sumpath=c()
found=F
nextStep=which(tree$edge[,2]==nextStep)
nextStart=tree$edge[nextStep,1]
while(! found){
  path=c(path,nextStart)
  sumpath=c(sumpath,tree$edge.length[nextStep])
  if (startNode==tree$edge[nextStep,1]){
    found=T
  }
  nextStep=which(tree$edge[,2]==nextStart)
  nextStart=tree$edge[nextStep,1]
}

colors=rep("black", length(tree$Nnode))
colors[tree$tip.label=="Aan"]="blue"
colors[tree$tip.label=="Cpe"]="red"
tree$node.states=rep("S", length(tree$Nnode))
tree$node.states[tree$tip.label=="Aan"]="R"
tree$node.states[tree$tip.label=="Cpe"]="R"

png("/media/merly/Baymax/research/cph-visit/img/inferred.phylogeny.png", 800,600)
plotTree(tree,mar=c(4.1,4.1,1.1,1.1), lwd=3, ftype="i", ylim=c(0,48))
axis(1); axis(2, at=1:48, labels=1:48, las=2)
# nodelabels(); tiplabels()
abline(v=sum(sumpath), lty="dashed", col="gray70")
abline(v=0, lty="dashed", col="gray70")
par(fg="black")
lastPP<-get("last_plot.phylo",env=.PlotPhyloEnv)
text(lastPP$xx[1:length(tree$tip.label)],
     lastPP$yy[1:length(tree$tip.label)],
     tree$tip.label,col=colors,
     pos=4,offset=0.3)
dev.off()


# Distance Matrix from inferred phylogeny
# treefile="/home/merly/git/cph-visit/coverage-analysis/capture-phylotenetic-decay/trees/RAxML_bipartitions.2011.concat"
# tree=read.nexus(treefile)
distMatrix=cophenetic(tree)
dCpe=distMatrix["Cpe",]
dAan=distMatrix["Aan",]

distMatrix[upper.tri(distMatrix)]=NA

rowColors=rep("black", nrow(distMatrix))
columnColors=rep("black", nrow(distMatrix))

rowColors[colnames(distMatrix)=="Cpe"]="red"
columnColors[colnames(distMatrix)=="Cpe"]="red"
rowColors[colnames(distMatrix)=="Aan"]="blue"
columnColors[colnames(distMatrix)=="Aan"]="blue"
png("/media/merly/Baymax/research/cph-visit/img/distance.matrix.inferred.phylogeny.png", 800,800)
heatmap.2(distMatrix, Rowv = NA,Colv = NA,
          key.title = "",
          key.xlab="Pairwise distance (tips)",
          key.ylab = "Frequency",
          scale="none",
          margins=c(3,1), # ("margin.Y", "margin.X")
          trace='none',
          symkey=FALSE,
          symbreaks=FALSE,
          dendrogram='none',
          density.info='histogram',
          denscol="black",
          keysize=0.3,
          colRow = rowColors,
          colCol = columnColors,
          #( "bottom.margin", "left.margin", "top.margin", "left.margin" )
          key.par=list(mar=c(5,0,3,0)),
          # lmat -- added 2 lattice sections (5 and 6) for padding
          lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(2.5, 5), lwid=c(1, 10, 1)
)
dev.off()

## Phylogenetic distance to outgroup reference (Swift)
dCpe=dCpe[names(dCpe)!="Cpe"]
qnt=quantile(dCpe,c(0.25, 0.75))
png("/media/merly/Baymax/research/cph-visit/img/distance.to.outgroup.png", 800,600)
plot(dCpe, type="l", col="darkred", axes=F, lwd=2, xlab="Samples",ylab="Distance")
abline(h=mean(dCpe), col="gray30")
abline(h=qnt[1], lty="dashed", col="gray60")
abline(h=qnt[2], lty="dashed", col="gray60")
axis(2); axis(1,at=1:length(dCpe), labels=names(dCpe), las=2)
dev.off()
## Phylogenetic distance to ingroup reference (Aan)
dAan=dAan[names(dAan)!="Aan"]
qnt=quantile(dAan,c(0.25, 0.75))
png("/media/merly/Baymax/research/cph-visit/img/distance.to.ingroup.png", 800,600)
plot(dAan, type="l", col="darkblue", axes=F, lwd=2,xlab="Samples",ylab="Distance")
abline(h=mean(dAan), col="gray30")
abline(h=qnt[1], lty="dashed", col="gray60")
abline(h=qnt[2], lty="dashed", col="gray60")
axis(2); axis(1,at=1:length(dAan), labels=names(dAan), las=2)
dev.off()

# Depth of coverage vs. phylogenetic distance
## Mapped to outgroup
WDSwift="/media/merly/Baymax/research/cph-visit/map2Swift/"
mergedMatrix=read.table(paste0(WDSwift,"files/coverage.matrix.txt"))
mergedMatrix=mergedMatrix[,! colnames(mergedMatrix) %in% c("AnnaBGI","SwiftBGI")]
mergedMatrix=mergedMatrix/targets$Size

meanPerSample=apply(mergedMatrix,2,mean)
mediansPerSample=apply(mergedMatrix,2,median)

cpe=dCpe[! names(dCpe) %in% c("Cpe","Aan")]
aan=dAan[! names(dAan) %in% c("Cpe","Aan")]
dataMap2OutgroupDis2In=cbind(aan,mediansPerSample[names(aan)])
dataMap2OutgroupDis2Out=cbind(cpe,mediansPerSample[names(cpe)])
labelsCorrect=c("H05","H01","H08","H03","H04","H06","H02","H07")
labelsWronglyFormated=c("H5","H1","H8","H3","H4","H6","H2","H7")
namesIds=which(rownames(dataMap2OutgroupDis2In)%in%labelsCorrect)
rownames(dataMap2OutgroupDis2In)[namesIds]=labelsWronglyFormated
rownames(dataMap2OutgroupDis2Out)[namesIds]=labelsWronglyFormated
dataMap2OutgroupDis2In[labelsWronglyFormated,2]=mediansPerSample[labelsWronglyFormated]
dataMap2OutgroupDis2Out[labelsWronglyFormated,2]=mediansPerSample[labelsWronglyFormated]
quantiles=apply(mergedMatrix[,rownames(dataMap2OutgroupDis2In)],2,quantile,c(0.25,0.75))



maxDistanceX=max(dataMap2OutgroupDis2In[,1],dataMap2OutgroupDis2Out[,1])
maxCoverageY=max(dataMap2OutgroupDis2In[,2],dataMap2OutgroupDis2Out[,2])

png("/media/merly/Baymax/research/cph-visit/img/dist.vs.coverage.map2out.png", 1000,700)
layout(matrix(c(1,2),1,2,byrow=T))
plot(dataMap2OutgroupDis2In,
xlim=c(0,maxDistanceX),
ylim=c(0,0.02),
xlab="Phylogenetic distance",ylab="Depth of coverage",
pch=20, col=sampleColors, main = "Distance to ingroup")
arrows(dataMap2OutgroupDis2In[,1],quantiles[1,],dataMap2OutgroupDis2In[,1],quantiles[2,], length = 0.05,angle=90, code=3,col=sampleColors)
plot(dataMap2OutgroupDis2Out,
xlim=c(0,maxDistanceX),
ylim=c(0,0.02),
xlab="Phylogenetic distance",ylab="Depth of coverage",
pch=20, col=sampleColors, main = "Distance to outgroup")
arrows(dataMap2OutgroupDis2Out[,1],quantiles[1,],dataMap2OutgroupDis2Out[,1],quantiles[2,], length = 0.05,angle=90, code=3,col=sampleColors)
legend("topleft", legend=rownames(dataMap2OutgroupDis2Out), col=sampleColors, lty=1, lwd=2, cex=0.8)
dev.off()

## Mapped to ingroup
# General information
WDAan="/media/merly/Baymax/research/cph-visit/map2Aan/"
mergedMatrix=read.table(paste0(WDAan,"files/coverage.matrix.txt"))
mergedMatrix=mergedMatrix/targets$Size
meanPerSample=apply(mergedMatrix,2,mean)
mediansPerSample=apply(mergedMatrix,2,median)
cpe=dCpe[! names(dCpe) %in% c("Cpe","Aan")]
aan=dAan[! names(dAan) %in% c("Cpe","Aan")]
dataMap2IngroupDis2In=cbind(aan,mediansPerSample[names(aan)])
dataMap2IngroupDis2Out=cbind(cpe,mediansPerSample[names(cpe)])
labelsCorrect=c("H05","H01","H08","H03","H04","H06","H02","H07")
labelsWronglyFormated=c("H5","H1","H8","H3","H4","H6","H2","H7")
namesIds=which(rownames(dataMap2IngroupDis2In)%in%labelsCorrect)
rownames(dataMap2IngroupDis2In)[namesIds]=labelsWronglyFormated
rownames(dataMap2IngroupDis2Out)[namesIds]=labelsWronglyFormated
dataMap2IngroupDis2In[labelsWronglyFormated,2]=mediansPerSample[labelsWronglyFormated]
dataMap2IngroupDis2Out[labelsWronglyFormated,2]=mediansPerSample[labelsWronglyFormated]
quantiles=apply(mergedMatrix[,rownames(dataMap2IngroupDis2In)],2,quantile,c(0.25,0.75))

maxDistanceX=max(dataMap2IngroupDis2In[,1],dataMap2IngroupDis2Out[,1])
maxCoverageY=max(quantiles)
png("/media/merly/Baymax/research/cph-visit/img/dist.vs.coverage.map2ingroup.png", 1000,700)
layout(matrix(c(1,2),1,2,byrow=T))
plot(dataMap2IngroupDis2In,
xlim=c(0,maxDistanceX),
ylim=c(0,maxCoverageY),
xlab="Phylogenetic distance",ylab="Depth of coverage",
pch=20, col=sampleColors, main = "Distance to ingroup")
arrows(dataMap2IngroupDis2In[,1],quantiles[1,],dataMap2IngroupDis2In[,1],quantiles[2,], length = 0.05,angle=90, code=3,col=sampleColors)

plot(dataMap2IngroupDis2Out,
xlim=c(0,maxDistanceX),
ylim=c(0,maxCoverageY),
xlab="Phylogenetic distance",ylab="Depth of coverage",
pch=20, col=sampleColors, main = "Distance to outgroup")
arrows(dataMap2IngroupDis2Out[,1],quantiles[1,],dataMap2IngroupDis2Out[,1],quantiles[2,], length = 0.05,angle=90, code=3,col=sampleColors)
legend("topleft", legend=rownames(dataMap2IngroupDis2Out), col=sampleColors, lty=1, lwd=2, cex=0.8)
dev.off()
## Correlations

distance2In=dataMap2IngroupDis2In[,1]
distance2Out=dataMap2IngroupDis2Out[,1]
coverageMap2IngroupDis2In=dataMap2IngroupDis2In[,2]
coverageMap2IngroupDis2Out=dataMap2IngroupDis2Out[,2]
# SD= 0 for both cases
# coverageMap2OutgroupDis2In=dataMap2OutgroupDis2In[,2]
# coverageMap2OutgroupDis2Out=dataMap2OutgroupDis2In[,2]
cor.sp <- function(x,y) cor(rank(x),rank(y))
cor1=cor(distance2In, coverageMap2IngroupDis2In)
cor2=cor(distance2Out, coverageMap2IngroupDis2Out)
scor1=cor.sp(distance2In, coverageMap2IngroupDis2In)
scor2=cor.sp(distance2Out, coverageMap2IngroupDis2Out)
lm1=summary(lm(distance2In~coverageMap2IngroupDis2In))
lm2=summary(lm(distance2Out~coverageMap2IngroupDis2Out))

collabels=c(
  "Correlation Coefficient",
  "Corr. Squared",
  "Spearman corr."
)
rowlabels=c(
  "Distance to ingroup",
  "Distance to outgroup"
)

res=data.frame(x=c(cor1,cor2), xx=c(cor1^2,cor2^2),sp=c(scor1,scor2))
colnames(res)=collabels
rownames(res)=rowlabels
kable(res,format = "markdown")

write.table(res, "/media/merly/Baymax/research/cph-visit/files/correlation.coverage.vs.distance.txt")

#########################################################################################################
# targets not captured
matrixs=read.table("/media/merly/Baymax/research/cph-visit/map2Swift/files/coverage.matrix.txt")
matrixa=read.table("/media/merly/Baymax/research/cph-visit/map2Aan/files/coverage.matrix.txt")
captureda=matrixa
captureds=matrixs
captureda[matrixa>0]=1
captureds[matrixs>0]=1
captureda[is.na(captureda)]=0
captureds[is.na(captureds)]=0
notcaptureda=1-captureda
notcaptureds=1-captureds
notboth=notcaptureda+notcaptureds[,!colnames(notcaptureds) %in% c("AnnaBGI","SwiftBGI")]
notcapturedpersample=colSums(notboth)
notcapturedpertarget=rowSums(notboth)
png("/media/merly/Baymax/research/cph-visit/img/notcaptured.rel.numsamples.both.png", 800,600)
layout(matrix(c(1,2),2,1))
h=hist(notcapturedpertarget/2, breaks=47, plot=F)
plot(h$counts,
     type="l", lwd=2,col="darkred",
     ylim=c(0,500),
     xlab="Number of samples", ylab="Number of targets", axes=F)
axis(2); axis(1, at=1:46, labels=1:46,las=2 ,cex=0.8)
dev.off()

png("/media/merly/Baymax/research/cph-visit/img/notcaptured.rel.numsamples.png", 800,800)
layout(matrix(c(1,2),2,1))
h1=hist(rowSums(captureda), breaks=47, plot=F)
h2=hist(rowSums(captureds), breaks=47, plot=F)
h3=hist(rowSums(notboth), breaks=47, plot=F)
plot(h1$counts,
     type="l", lwd=2,col="darkred",
     xlab="Number of samples", ylab="Number of targets", axes=F)
points(h2$counts,
     type="l", lwd=2,col="darkgreen",
     xlab="Number of samples", ylab="Number of targets")
points(h3$counts,
       type="l", lwd=2,col="darkblue",
       xlab="Number of samples", ylab="Number of targets")
legend("topleft", legend=c("map2Aan","map2Swift","Both"), col=c("darkred","darkgreen","darkblue"),lty=1, lwd=4)
axis(2); axis(1, at=1:46, labels=1:46,las=2 ,cex=0.8)
plot(h1$counts,
     type="l", lwd=2,col="darkred",
     ylim=c(0,500),
     xlab="Number of samples", ylab="Number of targets", axes=F)
points(h2$counts,
       type="l", lwd=2,col="darkgreen",
       ylim=c(0,500),
       xlab="Number of samples", ylab="Number of targets")
points(h3$counts,
       type="l", lwd=2,col="darkblue",
       ylim=c(0,500),
       xlab="Number of samples", ylab="Number of targets")
axis(2); axis(1, at=1:46, labels=1:46,las=2 ,cex=0.8)
dev.off()


untargetedBedFile="/media/merly/Baymax/research/cph-visit/files/untargeted.bed"
offtargets=read.table(untargetedBedFile,stringsAsFactors = F, colClasses = c("character","numeric","numeric"))
offtargets$V4=offtargets$V3-offtargets$V2
colnames(offtargets)=c("Scaffold","Start","End","Size")
offPerScaffold=split(offtargets,offtargets$Scaffold)
totalOffTargetsPerScaffold=sapply(offPerScaffold,nrow)

png(paste0("/media/merly/Baymax/research/cph-visit/img/","6.offtarget.1.num.offtargets.per.scaffold.png"),1366,768)
plot(totalOffTargetsPerScaffold,pch=16, col="darkred", axes=F, xlab="Scaffolds", ylab="Number of targeted regions")
axis(2); axis(1)
dev.off()
# Size distribution of the target regions
png(paste0(IMG,"6.offtarget.2.freq.sizes.png"), 1366,768)
h=hist(offtargets$Size,
       breaks =500,
       axes=F,
       col="darkgreen",border="white",
       main="Targeted regions size (zoom x-axis 0-1000)",
       xlab="Size (bp)",
       xlim=c(0,1000), prob=T
)
d=density(targets$Size)
points(d,lwd=2, type="l")
axis(1); axis(2)
dev.off()
prettyNum(sum(offtargets$Size), big.mark = ",")
