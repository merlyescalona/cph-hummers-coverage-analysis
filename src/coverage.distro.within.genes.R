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
coverageMatrixFileAnna=paste0(WD,"anna/files/coverage.matrix.per.gene.txt")
coverageMatrixFileSwift=paste0(WD,"swift/files/coverage.matrix.per.gene.txt")
coverageMatrixAnna=read.table(coverageMatrixFileAnna, header=T)/genesAnna$Size
coverageMatrixSwift=read.table(coverageMatrixFileSwift, header=T)/genesSwift$Size
rownames(coverageMatrixAnna)=genesAnna$Gene
rownames(coverageMatrixSwift)=genesSwift$Gene
#-------------------------------------------------------------------------------

coveragesGenesAnna=unlist(coverageMatrixAnna)
coveragesGenesSwift=unlist(coverageMatrixSwift)
d1=density(coveragesGenesAnna)
d2=density(coveragesGenesSwift)
maxY=roundUpNice(max(d1$y,d2$y))
plot(d1, type="l", col="darkred", lwd=2, ylim=c(0,maxY))
points(d2, type="l", col="darkblue", lwd=2)
