packages<-c("ape","geiger","apTreeshape","ggplot2","gplots","RColorBrewer","knitr","phangorn","futile.logger","phytools","gdata", "broom")
for(pkg in packages ){
  suppressMessages(library(pkg,character.only=TRUE,quietly=TRUE))
}
colors = colorRampPalette(brewer.pal(9, "Set1"))(48)
# Helpful functions: http://stackoverflow.com/a/6463946
roundUp <- function(x) 10^ceiling(log10(x))
roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10)) {
  if(length(x) != 1) stop("'x' must be of length 1")
  10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}
#-------------------------------------------------------------------------------
# Data loading
WD2="/home/merly/git/cph-visit/coverage-analysis/empirical/data/"
WD="/media/merly/Baymax/research/cph-visit/coverage-analysis/"
gffFileAnna=paste0(WD2,"general/Calypte_anna.gene.CDS.2750.3.gff")
gffFileSwift=paste0(WD2,"general/Chaetura_pelagica.CDS.2.gff")
gffAnna=read.table(gffFileAnna, stringsAsFactors=F,
                   colClasses=c("character","character","character","numeric",
                                "numeric","character","character","character","character"
                   ))
nGenesAnna=length(unique(gffAnna$V9))
genesAnna=data.frame(Gene=as.character(unique(gffAnna$V9)))
splitGFFAnna=split(gffAnna,gffAnna$V9)
genesAnna$Size=sapply(splitGFFAnna, function(x){ sum(x$V5-x$V4)})
gffSwift=read.table(gffFileSwift, stringsAsFactors=F,
                    colClasses=c("character","character","character","numeric",
                                 "numeric","character","character","character","character"
                    ))
nGenesSwift=length(unique(gffSwift$V9))
genesSwift=data.frame(Gene=as.character(unique(gffSwift$V9)))
splitGFFSwift=split(gffSwift,gffSwift$V9)
genesSwift$Size=sapply(splitGFFSwift, function(x){sum(x$V5-x$V4)})
#-------------------------------------------------------------------------------
samplesFilenameAnna=paste0(WD2,"anna/samples.txt")
samplesFilenameSwift=paste0(WD2,"swift/samples.txt")
sampleNamesAnna=unlist(read.table(samplesFilenameAnna, stringsAsFactors = F))
sampleNamesSwift=unlist(read.table(samplesFilenameSwift, stringsAsFactors = F))
#-------------------------------------------------------------------------------
coverageMatrixFileAnna=paste0(WD2,"anna/coverage.matrix.per.gene.txt")
coverageMatrixFileSwift=paste0(WD2,"swift/coverage.matrix.per.gene.txt")
coverageMatrixAnna=read.table(coverageMatrixFileAnna, header=T)/genesAnna$Size
coverageMatrixSwift=read.table(coverageMatrixFileSwift, header=T)/genesSwift$Size
rownames(coverageMatrixAnna)=genesAnna$Gene
rownames(coverageMatrixSwift)=genesSwift$Gene
#-------------------------------------------------------------------------------
geneList=unlist(read.table(paste0(WD2,"general/g.trees.to.filter.seqs.2.txt"), stringsAsFactors=F))
distancesAnna=matrix(rep(0,length(geneList)*length(sampleNamesSwift)), length(geneList), length(sampleNamesSwift))
distancesSwift=matrix(rep(0,length(geneList)*length(sampleNamesSwift)), length(geneList), length(sampleNamesSwift))
rownames(distancesAnna)=geneList
colnames(distancesAnna)=sampleNamesSwift
rownames(distancesSwift)=geneList
colnames(distancesSwift)=sampleNamesSwift
#-------------------------------------------------------------------------------
# Need to filter coverage per gene with the filtered genes
#-------------------------------------------------------------------------------

gtreefiles <- list.files(
  path = paste0(WD,"trees/seqs.2/"),
  pattern="*.x.aa.aln.backTranslated")
gtreefiles=paste0(WD,"trees/seqs.2/",gtreefiles)
labelsCorrect=c("H05","H01","H08","H03","H04","H06","H02","H07")
labelsWronglyFormated=c("H5","H1","H8","H3","H4","H6","H2","H7")

for (i in 1:length(gtreefiles)) {
  print(i)
  msa=as.matrix(read.FASTA(gzfile(gtreefiles[i])))
  hammDist=as.matrix(dist.dna(msa, model="raw"))
  indexAnna=which(startsWith(rownames(hammDist),"Aan"))
  indexSwift=which(startsWith(rownames(hammDist),"Cpe"))
  namesIds=which(rownames(hammDist) %in% labelsCorrect)
  rownames(hammDist)[c(indexAnna, indexSwift, namesIds)]=c("AnnaBGI","SwiftBGI",labelsWronglyFormated)
  colnames(hammDist)[c(indexAnna, indexSwift, namesIds)]=c("AnnaBGI","SwiftBGI",labelsWronglyFormated)
  gene=strsplit(basename(gtreefiles[i]),"\\.")[[1]][1]
  distancesAnna[gene,colnames(distancesAnna)]=hammDist[indexAnna,colnames(distancesAnna)]
  distancesSwift[gene,colnames(distancesSwift)]=hammDist[indexSwift,colnames(distancesSwift)]
}

write.table(distancesAnna,file=paste0(WD,"files/distance.matrix.per.gene.anna.txt"), col.names=T,row.names=T)
write.table(distancesSwift,file=paste0(WD,"files/distance.matrix.per.gene.swift.txt"), col.names=T,row.names=T)

distancesAnna=read.table(file=paste0(WD2,"general/distance.matrix.per.gene.anna.txt"), header=T)
distancesSwift=read.table(file=paste0(WD2,"general/distance.matrix.per.gene.swift.txt"), header=T)
birds48filename=paste0(WD2,"general/48birds_ortholog.list.chi.anna.cpe.hum.finch")
birds=read.table(birds48filename, header=T,
    colClasses=c("character","numeric","character","character",
      "character","character","character"))

birds.2=birds[birds$anna!="-",]
birds.3=birds.2[birds.2$swift!="-",]
birds.4=unique(birds.3)
birds.5=birds.4[,c("anna","swift")]
birds.5$swift=paste0("Parent=",birds.5$swift,";")
birds.5$anna=paste0("Parent=Aan_", birds.5$anna,";")
annaGenes=unique(gffAnna$V9)
birds.6=birds.5[birds.5$anna %in% annaGenes, ]
birds.6$swift.plain=unlist(strsplit(unlist(lapply(strsplit(birds.6$swift,"_"), function(x){x[[2]]})), ";"))
birds.6$anna.plain=unlist(strsplit(unlist(lapply(strsplit(birds.6$anna,"_"), function(x){x[[2]]})), ";"))
smallbird=birds.6[birds.6$anna.plain %in% intersect(birds.6$anna.plain,geneList),]

distancesAnna=distancesAnna[smallbird$anna.plain,]
distancesSwift=distancesSwift[smallbird$anna.plain,]
indexNAAnn=which(!is.na(distancesAnna[,2]))
indexNASwift=which(!is.na(distancesSwift[,2]))
indexNA=intersect(indexNA, indexNASwift)
distancesAnna=distancesAnna[indexNA,]
distancesSwift=distancesSwift[indexNA,]
rownames(distancesAnna)=smallbird$anna[indexNA]
rownames(distancesSwift)=smallbird$swift[indexNA]

cmswift=coverageMatrixSwift[smallbird$swift[indexNA],]
cmanna=coverageMatrixAnna[smallbird$anna[indexNA],]
cmanna=cbind(cmanna,AnnaBGI=rep(0,nrow(cmanna)), SwiftBGI=rep(0,nrow(cmanna)))

smallbird=smallbird[indexNA,]

pointsCorrSwift=data.frame(dataset=c(),gene=c(),sample=c(),distance=c(),coverage=c())
for(geneIndex in 1:length(smallbird$swift)){
  gene=smallbird$swift[geneIndex]
  for (sampleIndex in 1:length(sampleNamesSwift)){
    sample=sampleNamesSwift[sampleIndex]
    pointsCorrSwift=rbind(pointsCorrSwift,
                     cbind(
                       "Swift",
                       gene,sample,
                       distancesSwift[gene,sample],
                       cmswift[gene,sample])
    )
  }
}
p1=as.numeric(as.character(pointsCorrSwift[,4]))
p2=as.numeric(as.character(pointsCorrSwift[,5]))
# removing distance outliers with paired coverage
qntSwift=quantile(p1,probs=c(0.25,0.75))
HSwift=1.5*IQR(p1)
# removing coverage outliers with paired distance
qntSwift2=quantile(p2,probs=c(0.25,0.75))
HSwift2=1.5*IQR(p2)
finalDistanceSwift=p1
finalDistanceSwift[p1 < (qntSwift[1] - HSwift)]=NA
finalDistanceSwift[p1 > (qntSwift[2] + HSwift)]=NA
finalDistanceSwift[p2 < (qntSwift2[1] - HSwift2)]=NA
finalDistanceSwift[p2 > (qntSwift2[2] + HSwift2)]=NA
indicesNA=is.na(finalDistanceSwift)
finalDistanceSwift=finalDistanceSwift[!indicesNA]
finalCoverageSwift=p2[!indicesNA]

pointsCorrAnna=c()
for(geneIndex in 1:length(smallbird$anna)){
  gene=smallbird$anna[geneIndex]
  for (sampleIndex in 1:length(sampleNamesAnna)){
    sample=sampleNamesAnna[sampleIndex]
    pointsCorrAnna=rbind(pointsCorrAnna,
                     cbind(
                       "Anna",
                       gene,sample,
                       distancesAnna[gene,sample],
                       cmanna[gene,sample])
    )
  }
}
p11=as.numeric(pointsCorrAnna[,4])
p21=as.numeric(pointsCorrAnna[,5])

qntAnna=quantile(p11,probs=c(0.25,0.75))
qntAnna2=quantile(p21,probs=c(0.25,0.75))
HAnna=1.5*IQR(p11)
HAnna2=1.5*IQR(p21)
finalDistanceAnna=p11
finalDistanceAnna[p11 < (qntAnna[1] - HAnna)]=NA
finalDistanceAnna[p11 > (qntAnna[2] + HAnna)]=NA
finalDistanceAnna[p21 < (qntAnna2[1] - HAnna2)]=NA
finalDistanceAnna[p21 > (qntAnna2[2] + HAnna2)]=NA
indicesNA=is.na(finalDistanceAnna)
finalDistanceAnna=finalDistanceAnna[!indicesNA]
finalCoverageAnna=p21[!indicesNA]


lmAnna=lm(finalCoverageSwift~finalDistanceSwift)
lmSwift=lm(finalCoverageAnna~finalDistanceAnna)

sink(file=paste0(WD2,"anna/lm.dist.cov.phylo.txt"))
summary(lmAnna)
sink(NULL)
sink(file=paste0(WD2,"swift/lm.dist.cov.phylo.txt"))
summary(lmSwift)
sink(NULL)

png(paste0(WD2,"anna/phylo.distance.coverage.genes.lm.4.png"), 900,800)
layout(matrix(c(1,2,3,4), 2,2, byrow=T))
plot(lmAnna, pch=16)
dev.off()
png(paste0(WD2,"swift/phylo.distance.coverage.genes.lm.4.png"), 900,800)
layout(matrix(c(1,2,3,4), 2,2, byrow=T))
plot(lmSwift, pch=16)
dev.off()
#####################################################################
#####################################################################
gene=smallbird$swift[1]
sample=1
maxY=roundUpNice(max(cmswift))
maxX=roundUpNice(max(distancesSwift))sd
png(paste0(WD2,"swift/phylo.distance.coverage.genes.png"), 900,800)
plot(
  distancesSwift[gene,sampleNamesSwift[sample]],
  cmswift[gene,sampleNamesSwift[sample]],
  col=colors[sample], pch=16,
  ylim=c(0,maxY),xlim=c(0,maxX),
  xlab="Distance", ylab="Coverage",
  main="map2Swift",
  axes=F, lwd=2
)
for(geneIndex in 1:length(smallbird$swift)){
  gene=smallbird$swift[geneIndex]
  for (sampleIndex in 2:length(sampleNamesSwift)){
    sample=sampleNamesSwift[sampleIndex]
    points(
      distancesSwift[gene,sample],
      cmswift[gene,sample],
      pch=16,col=colors[sampleIndex], lwd=2
    )
  }
}
rect(
  0,0,
  (qntSwift[2] + HSwift),
  (qntSwift2[2] + HSwift2),
  lwd=2,lty="solid"
)
abline(lmSwift, lwd=2, col="red")
axis(2,las=2)
axis(1,las=2)
legend("topright", legend=sampleNamesSwift,col=colors, ncol=6, pch=20)
legend("topleft", legend=c("Outliers","LM"),col=c("black","red"),lwd=2,lty=c("solid","solid"))
dev.off()

png(paste0(WD2,"anna/phylo.distance.coverage.genes.png"), 900,800)
gene=smallbird$anna[1]
sample=1
maxY=roundUpNice(max(cmanna))
maxX=roundUpNice(max(distancesAnna))
plot(
  distancesAnna[gene,sampleNamesAnna[sample]],
  cmanna[gene,sampleNamesAnna[sample]],
  col=colors[sample], pch=15,
  ylim=c(0,maxY),xlim=c(0,maxX),
  xlab="Distance", ylab="Coverage",
  main="map2anna",
  axes=F, lwd=2
)
for(geneIndex in 1:length(smallbird$anna)){
  gene=smallbird$anna[geneIndex]
  for (sampleIndex in 2:length(sampleNamesAnna)){
    sample=sampleNamesAnna[sampleIndex]
    points(
      distancesAnna[gene,sample],
      cmanna[gene,sample],
      pch=15,col=colors[sampleIndex], lwd=3
    )
  }
}

legend("topright", legend=sampleNamesSwift,col=colors, ncol=6, pch=20)
legend("topleft", legend=c("Outliers","LM"),col=c("black","red"),lwd=2,lty=c("solid","solid"))
rect(
  0,0,
  (qntAnna[2] + HAnna),
  (qntAnna2[2] + HAnna2),
  lwd=2,lty="solid"
)
abline(lmAnna, lwd=2, col="red")

axis(2,las=2)
axis(1,las=2)
dev.off()


pointsCorrAnna=data.frame(pointsCorrAnna, stringsAsFactors=F)
colnames(pointsCorrAnna)=c("Dataset","Gene","Sample","Distance","Coverage")
pointsCorrSwift=data.frame(pointsCorrSwift, stringsAsFactors=F)
colnames(pointsCorrSwift)=c("Dataset","Gene","Sample","Distance","Coverage")
cAnna=pointsCorrAnna[pointsCorrAnna$Sample %in% c("H21", "H47"), ]
cSwift=pointsCorrSwift[pointsCorrSwift$Sample %in% c("H21", "H47"), ]
bAnna=pointsCorrAnna[pointsCorrAnna$Sample %in% c("H40", "H46"), ]
bSwift=pointsCorrSwift[pointsCorrSwift$Sample %in% c("H40", "H46"), ]
mAnna=pointsCorrAnna[pointsCorrAnna$Sample %in% c("H8", "H24"), ]
mSwift=pointsCorrSwift[pointsCorrSwift$Sample %in% c("H8", "H24"), ]
hAnna=pointsCorrAnna[pointsCorrAnna$Sample %in% c("H2", "H16"), ]
hSwift=pointsCorrSwift[pointsCorrSwift$Sample %in% c("H2", "H16"), ]

cAnna$Distance=as.numeric(as.character(cAnna$Distance)); cAnna$Coverage=as.numeric(as.character(cAnna$Coverage))
cSwift$Distance=as.numeric(as.character(cSwift$Distance)); cSwift$Coverage=as.numeric(as.character(cSwift$Coverage))
bAnna$Distance=as.numeric(as.character(bAnna$Distance)); bAnna$Coverage=as.numeric(as.character(bAnna$Coverage))
bSwift$Distance=as.numeric(as.character(bSwift$Distance)); bSwift$Coverage=as.numeric(as.character(bSwift$Coverage))
mAnna$Distance=as.numeric(as.character(mAnna$Distance)); mAnna$Coverage=as.numeric(as.character(mAnna$Coverage))
mSwift$Distance=as.numeric(as.character(mSwift$Distance)); mSwift$Coverage=as.numeric(as.character(mSwift$Coverage))
hAnna$Distance=as.numeric(as.character(hAnna$Distance)); hAnna$Coverage=as.numeric(as.character(hAnna$Coverage))
hSwift$Distance=as.numeric(as.character(hSwift$Distance)); hSwift$Coverage=as.numeric(as.character(hSwift$Coverage))

cAnnalmH21=lm(cAnna[cAnna$Sample=="H21",]$Coverage~cAnna[cAnna$Sample=="H21",]$Distance)
cAnnalmH47=lm(cAnna[cAnna$Sample=="H47",]$Coverage~cAnna[cAnna$Sample=="H47",]$Distance)
cSwiftlmH21=lm(cSwift[cSwift$Sample=="H21",]$Coverage~cSwift[cSwift$Sample=="H21",]$Distance)
cSwiftlmH47=lm(cSwift[cSwift$Sample=="H47",]$Coverage~cSwift[cSwift$Sample=="H47",]$Distance)
bAnnalmH40=lm(bAnna[bAnna$Sample=="H40",]$Coverage~bAnna[bAnna$Sample=="H40",]$Distance)
bAnnalmH46=lm(bAnna[bAnna$Sample=="H46",]$Coverage~bAnna[bAnna$Sample=="H46",]$Distance)
bSwiftlmH40=lm(bSwift[bSwift$Sample=="H40",]$Coverage~bSwift[bSwift$Sample=="H40",]$Distance)
bSwiftlmH46=lm(bSwift[bSwift$Sample=="H46",]$Coverage~bSwift[bSwift$Sample=="H46",]$Distance)
mAnnalmH8=lm(mAnna[mAnna$Sample=="H8",]$Coverage~mAnna[mAnna$Sample=="H8",]$Distance)
mAnnalmH24=lm(mAnna[mAnna$Sample=="H24",]$Coverage~mAnna[mAnna$Sample=="H24",]$Distance)
mSwiftlmH8=lm(mSwift[mSwift$Sample=="H8",]$Coverage~mSwift[mSwift$Sample=="H8",]$Distance)
mSwiftlmH24=lm(mSwift[mSwift$Sample=="H24",]$Coverage~mSwift[mSwift$Sample=="H24",]$Distance)
hAnnalmH2=lm(hAnna[hAnna$Sample=="H2",]$Coverage~hAnna[hAnna$Sample=="H2",]$Distance)
hAnnalmH16=lm(hAnna[hAnna$Sample=="H16",]$Coverage~hAnna[hAnna$Sample=="H16",]$Distance)
hSwiftlmH2=lm(hSwift[hSwift$Sample=="H2",]$Coverage~hSwift[hSwift$Sample=="H2",]$Distance)
hSwiftlmH16=lm(hSwift[hSwift$Sample=="H16",]$Coverage~hSwift[hSwift$Sample=="H16",]$Distance)

pvalues=data.frame(
    Coquettes=c(glance(cAnnalmH21)$p.value,glance(cAnnalmH47)$p.value,glance(cSwiftlmH21)$p.value,glance(cSwiftlmH47)$p.value),
    Brilliants=c(glance(bAnnalmH40)$p.value, glance(bAnnalmH46)$p.value, glance(bSwiftlmH40)$p.value, glance(bSwiftlmH46)$p.value),
    Mangoes=c(glance(mAnnalmH8)$p.value, glance(mAnnalmH24)$p.value, glance(mSwiftlmH8)$p.value, glance(mSwiftlmH24)$p.value),
    Hermits=c(glance(hAnnalmH2)$p.value, glance(hAnnalmH16)$p.value, glance(hSwiftlmH2)$p.value, glance(hSwiftlmH16)$p.value )
)
rownames(pvalues)=c(
  "Low Coverage - map2Anna",
  "High Coverage - map2Anna",
  "Low Coverage - map2Swift",
  "High Coverage - map2Swift"
)
write.table(pvalues, "/home/merly/git/cph-visit/coverage-analysis/empirical/data/general/phylo.coverage.lm.pvalues.txt",
  row.names=T,
  col.names=T
)

maxX=roundUpNice(max(cAnna$Distance,cSwift$Distance,bAnna$Distance,bSwift$Distance,mAnna$Distance,mSwift$Distance,hAnna$Distance,hSwift$Distance))
maxY=roundUpNice(max(cAnna$Coverage, cSwift$Coverage, bAnna$Coverage, bSwift$Coverage, mAnna$Coverage, mSwift$Coverage, hAnna$Coverage, hSwift$Coverage))

png(paste0("/home/merly/git/cph-visit/coverage-analysis/empirical/data/general/corr.high.low.coverage.phylodist.vs.cov.coquettes.png"), 768,364)
layout(matrix(c(1,2),1,2, byrow=T))
plot(cAnna[cAnna$Sample=="H21",]$Distance,cAnna[cAnna$Sample=="H21",]$Coverage,col=colors[which(sampleNamesSwift=="H21")],pch=16,xlim=c(0,maxX), ylim=c(0,maxY),  lwd=2,xlab="Distance", ylab="Coverage",main=" C -  Anna ")
points(cAnna[cAnna$Sample=="H47",]$Distance,cAnna[cAnna$Sample=="H47",]$Coverage, col=colors[which(sampleNamesSwift=="H47")],pch=15, lwd=2)
abline(cAnnalmH21,lwd=2,col=colors[which(sampleNamesSwift=="H21")],lty="solid")
abline(cAnnalmH47,lwd=2,col=colors[which(sampleNamesSwift=="H47")],lty="solid")
legend("topright",legend=c("H21 - H","H47 - L"),pch=c(16,15),col=c(col=colors[which(sampleNamesSwift=="H21")],col=colors[which(sampleNamesSwift=="H47")]))
legend("topleft", horiz=T,legend=c("H21 - lm","H47- lm"), lwd=2,col=c(col=colors[which(sampleNamesSwift=="H21")],col=colors[which(sampleNamesSwift=="H47")]))
plot(cSwift[cSwift$Sample=="H21",]$Distance,cSwift[cSwift$Sample=="H21",]$Coverage,col=colors[which(sampleNamesSwift=="H21")],pch=16,xlim=c(0,maxX), ylim=c(0,maxY),  lwd=2,xlab="Distance", ylab="Coverage",main=" C - Swift")
points(cSwift[cSwift$Sample=="H47",]$Distance,cSwift[cSwift$Sample=="H47",]$Coverage, col=colors[which(sampleNamesSwift=="H47")],pch=15, lwd=2)
abline(cSwiftlmH21,lwd=2,col=colors[which(sampleNamesSwift=="H21")],lty="solid")
abline(cSwiftlmH47,lwd=2,col=colors[which(sampleNamesSwift=="H47")],lty="solid")
legend("topright",legend=c("H21 - H","H47 - L"),pch=c(16,15),col=c(col=colors[which(sampleNamesSwift=="H21")],col=colors[which(sampleNamesSwift=="H47")]))
legend("topleft", horiz=T,legend=c("H21 - lm","H47- lm"), lwd=2,col=c(col=colors[which(sampleNamesSwift=="H21")],col=colors[which(sampleNamesSwift=="H47")]))
dev.off()
png(paste0("/home/merly/git/cph-visit/coverage-analysis/empirical/data/general/corr.high.low.coverage.phylodist.vs.cov.brilliants.png"), 768,364)
layout(matrix(c(1,2),1,2, byrow=T))
plot(bAnna[bAnna$Sample=="H40",]$Distance,bAnna[bAnna$Sample=="H40",]$Coverage,col=colors[which(sampleNamesSwift=="H40")],pch=16,xlim=c(0,maxX), ylim=c(0,maxY),  lwd=2,xlab="Distance", ylab="Coverage",main=" B -  Anna ")
points(bAnna[bAnna$Sample=="H46",]$Distance,bAnna[bAnna$Sample=="H46",]$Coverage, col=colors[which(sampleNamesSwift=="H46")],pch=15, lwd=2)
abline(bAnnalmH40,lwd=2,col=colors[which(sampleNamesSwift=="H40")],lty="solid")
abline(bAnnalmH46,lwd=2,col=colors[which(sampleNamesSwift=="H46")],lty="solid")
legend("topright",legend=c("H40 - H","H46 - L"),pch=c(16,15),col=c(col=colors[which(sampleNamesSwift=="H40")],col=colors[which(sampleNamesSwift=="H46")]))
legend("topleft", horiz=T,legend=c("H40 - lm","H46- lm"), lwd=2,col=c(col=colors[which(sampleNamesSwift=="H40")],col=colors[which(sampleNamesSwift=="H46")]))
plot(bSwift[bSwift$Sample=="H40",]$Distance,bSwift[bSwift$Sample=="H40",]$Coverage,col=colors[which(sampleNamesSwift=="H40")],pch=16,xlim=c(0,maxX), ylim=c(0,maxY),  lwd=2,xlab="Distance", ylab="Coverage",main=" B - Swift")
points(bSwift[bSwift$Sample=="H46",]$Distance,bSwift[bSwift$Sample=="H46",]$Coverage, col=colors[which(sampleNamesSwift=="H46")],pch=15, lwd=2)
abline(bSwiftlmH40,lwd=2,col=colors[which(sampleNamesSwift=="H40")],lty="solid")
abline(bSwiftlmH46,lwd=2,col=colors[which(sampleNamesSwift=="H46")],lty="solid")
legend("topright",legend=c("H40 - H","H46 - L"),pch=c(16,15),col=c(col=colors[which(sampleNamesSwift=="H40")],col=colors[which(sampleNamesSwift=="H46")]))
legend("topleft", horiz=T,legend=c("H40 - lm","H46- lm"), lwd=2,col=c(col=colors[which(sampleNamesSwift=="H40")],col=colors[which(sampleNamesSwift=="H46")]))
dev.off()
png(paste0("/home/merly/git/cph-visit/coverage-analysis/empirical/data/general/corr.high.low.coverage.phylodist.vs.cov.mangoes.png"), 768,364)
layout(matrix(c(1,2),1,2, byrow=T))
plot(mAnna[mAnna$Sample=="H8",]$Distance,mAnna[mAnna$Sample=="H8",]$Coverage,col=colors[which(sampleNamesSwift=="H8")],pch=16,xlim=c(0,maxX), ylim=c(0,maxY),  lwd=2,xlab="Distance", ylab="Coverage",main=" M -  Anna ")
points(mAnna[mAnna$Sample=="H24",]$Distance,mAnna[mAnna$Sample=="H24",]$Coverage, col=colors[which(sampleNamesSwift=="H24")],pch=15, lwd=2)
abline(mAnnalmH8,lwd=2,col=colors[which(sampleNamesSwift=="H8")],lty="solid")
abline(mAnnalmH24,lwd=2,col=colors[which(sampleNamesSwift=="H24")],lty="solid")
legend("topright",legend=c("H8 - H","H24 - L"),pch=c(16,15),col=c(col=colors[which(sampleNamesSwift=="H8")],col=colors[which(sampleNamesSwift=="H24")]))
legend("topleft", horiz=T,legend=c("H8 - lm","H24- lm"), lwd=2,col=c(col=colors[which(sampleNamesSwift=="H8")],col=colors[which(sampleNamesSwift=="H24")]))
plot(mSwift[mSwift$Sample=="H8",]$Distance,mSwift[mSwift$Sample=="H8",]$Coverage,col=colors[which(sampleNamesSwift=="H8")],pch=16,xlim=c(0,maxX), ylim=c(0,maxY),  lwd=2,xlab="Distance", ylab="Coverage",main=" M - Swift")
points(mSwift[mSwift$Sample=="H24",]$Distance,mSwift[mSwift$Sample=="H24",]$Coverage, col=colors[which(sampleNamesSwift=="H24")],pch=15, lwd=2)
abline(mSwiftlmH8,lwd=2,col=colors[which(sampleNamesSwift=="H8")],lty="solid")
abline(mSwiftlmH24,lwd=2,col=colors[which(sampleNamesSwift=="H24")],lty="solid")
legend("topright",legend=c("H8 - H","H24 - L"),pch=c(16,15),col=c(col=colors[which(sampleNamesSwift=="H8")],col=colors[which(sampleNamesSwift=="H24")]))
legend("topleft", horiz=T,legend=c("H8 - lm","H24- lm"), lwd=2,col=c(col=colors[which(sampleNamesSwift=="H8")],col=colors[which(sampleNamesSwift=="H24")]))
dev.off()
png(paste0("/home/merly/git/cph-visit/coverage-analysis/empirical/data/general/corr.high.low.coverage.phylodist.vs.cov.hermits.png"), 768,364)
layout(matrix(c(1,2),1,2, byrow=T))
plot(hAnna[hAnna$Sample=="H2",]$Distance,hAnna[hAnna$Sample=="H2",]$Coverage,col=colors[which(sampleNamesSwift=="H2")],pch=16,xlim=c(0,maxX), ylim=c(0,maxY),  lwd=2,xlab="Distance", ylab="Coverage",main=" H - Anna  ")
points(hAnna[hAnna$Sample=="H16",]$Distance,hAnna[hAnna$Sample=="H16",]$Coverage, col=colors[which(sampleNamesSwift=="H16")],pch=15, lwd=2)
abline(hAnnalmH2,lwd=2,col=colors[which(sampleNamesSwift=="H2")],lty="solid")
abline(hAnnalmH16,lwd=2,col=colors[which(sampleNamesSwift=="H16")],lty="solid")
legend("topright",legend=c("H2 - H","H16 - L"),pch=c(16,15),col=c(col=colors[which(sampleNamesSwift=="H2")],col=colors[which(sampleNamesSwift=="H16")]))
legend("topleft", horiz=T,legend=c("H2 - lm","H16- lm"), lwd=2,col=c(col=colors[which(sampleNamesSwift=="H2")],col=colors[which(sampleNamesSwift=="H16")]))
plot(hSwift[hSwift$Sample=="H2",]$Distance,hSwift[hSwift$Sample=="H2",]$Coverage,col=colors[which(sampleNamesSwift=="H2")],pch=16,xlim=c(0,maxX), ylim=c(0,maxY),  lwd=2,xlab="Distance", ylab="Coverage",main=" H - Swift")
points(hSwift[hSwift$Sample=="H16",]$Distance,hSwift[hSwift$Sample=="H16",]$Coverage, col=colors[which(sampleNamesSwift=="H16")],pch=15, lwd=2)
abline(hSwiftlmH2,lwd=2,col=colors[which(sampleNamesSwift=="H2")],lty="solid")
abline(hSwiftlmH16,lwd=2,col=colors[which(sampleNamesSwift=="H16")],lty="solid")
legend("topright",legend=c("H2 - H","H16 - L"),pch=c(16,15),col=c(col=colors[which(sampleNamesSwift=="H2")],col=colors[which(sampleNamesSwift=="H16")]))
legend("topleft", horiz=T,legend=c("H2 - lm","H16- lm"), lwd=2,col=c(col=colors[which(sampleNamesSwift=="H2")],col=colors[which(sampleNamesSwift=="H16")]))
dev.off()
