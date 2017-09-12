packages<-c("ape","geiger","apTreeshape","ggplot2","gplots","RColorBrewer","knitr","phangorn","futile.logger","phytools")
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
WD="/media/merly/Baymax/research/cph-visit/coverage-analysis/"
gffFileAnna=paste0(WD,"files/Calypte_anna.gene.CDS.2750.3.gff")
gffFileSwift=paste0(WD,"files/Chaetura_pelagica.CDS.2.gff")
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
samplesFilenameAnna=paste0(WD,"anna/files/samples.txt")
samplesFilenameSwift=paste0(WD,"swift/files/samples.txt")
sampleNamesAnna=unlist(read.table(samplesFilenameAnna, stringsAsFactors = F))
sampleNamesSwift=unlist(read.table(samplesFilenameSwift, stringsAsFactors = F))
#-------------------------------------------------------------------------------
targetFileAnna=paste0(WD2,"files/targets.anna.3.bed")
targetFileSwift=paste0(WD2,"files/targets.swift.2.bed")
targetsAnna=read.table(targetFileAnna, stringsAsFactors=F)
targetsSwift=read.table(targetFileSwift, stringsAsFactors=F)
targetsAnna$V4=targetsAnna$V3-targetsAnna$V2
targetsSwift$V4=targetsSwift$V3-targetsSwift$V2
colnames(targetsAnna)=c("Scaffold","Start","End","Size")
colnames(targetsSwift)=c("Scaffold","Start","End","Size")
rownames(targetsAnna)=paste0(targetsAnna$Scaffold,
  rep("-", nrow(targetsAnna)),
  targetsAnna$Start,
  rep("-", nrow(targetsAnna)),
  targetsAnna$End)
rownames(targetsSwift)=paste0(targetsSwift$Scaffold,
  rep("-", nrow(targetsSwift)),
  targetsSwift$Start,
  rep("-", nrow(targetsSwift)),
  targetsSwift$End)


coverageMatrixFileAnna=paste0(WD,"anna/files/coverage.matrix.per.target.txt")
coverageMatrixFileSwift=paste0(WD,"swift/files/coverage.matrix.per.target.txt")
coverageMatrixAnna=read.table(coverageMatrixFileAnna, header=T)/genesAnna$Size
coverageMatrixSwift=read.table(coverageMatrixFileSwift, header=T)/genesSwift$Size
rownames(coverageMatrixAnna)=genesAnna$Gene
rownames(coverageMatrixSwift)=genesSwift$Gene
#-------------------------------------------------------------------------------

distancesAnna=matrix(rep(0,length(genesAnna$Gene)*length(sampleNamesSwift)), length(genesAnna$Gene), length(sampleNamesSwift))
distancesSwift=matrix(rep(0,length(genesSwift$Gene)*length(sampleNamesSwift)), length(genesSwift$Gene), length(sampleNamesSwift))
rownames(distancesAnna)=genesAnna$Gene
colnames(distancesAnna)=sampleNamesSwift
rownames(distancesSwift)=genesSwift$Gene
colnames(distancesSwift)=sampleNamesSwift
filesWO <- list.files(
  path = paste0(WD,"bedtools2/nohist/"),
  pattern="*.gz$")
filesWO=paste0(WD,"bedtools2/nohist/",filesWO)
mergedMatrix=matrix(0, nrow=nrow(targets),ncol = length(sampleNames))
mergedMatrixGenes=matrix(0, nrow=nrow(genes),ncol = length(sampleNames))


gtreefiles <- list.files(
  path = paste0(WD,"trees/gtrees.2/"),
  pattern="^R")
  gtreefiles=paste0(WD,"trees/gtrees.2/",gtreefiles)
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
  colnames(data)=c("scaffold","start","end","gene","pos","cov1","size","breadth", "targetname")
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
