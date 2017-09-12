# This has been developed to work in triploid.uvigo.es
# Folder cph-bam-coverage
currentDirectory="/home/merly/cph-bam-coverage/"
setwd(currentDirectory)
genomeCovDir=paste0(currentDirectory,"genomecov")
targetsFilename=paste0(currentDirectory,"files/","targets.bed")
samplesFilename=paste0(currentDirectory,"files/","samples.txt")
scaffoldsFilename=paste0(currentDirectory,"files/","scaffolds.txt")

targets=read.table(
    targetsFilename,
    stringsAsFactors = F
  )
scaffolds=unlist(read.table(
  scaffoldsFilename,
  stringsAsFactors = F
))
samples=unlist(read.table(
  samplesFilename,
  stringsAsFactors = F
))

shortlist=samples[c(2:11,13:14)]


targetsPerScaffold=split(targets,targets$V1)

# Sinlg esample analysis
sample="H09"
sampleFilename=paste0(
  rep(genomeCovDir, length(scaffolds)),
  rep("/", length(scaffolds)),
  rep(sample, length(scaffolds)),
  rep("/", length(scaffolds)),
  rep(sample, length(scaffolds)),
  rep(".", length(scaffolds)),
  scaffolds,
  rep(".genomecov", length(scaffolds))
  )

coverageData=list()
for(i in 1:length(sampleFilename)){
  coverageData[[i]]=read.table(sampleFilename[i])
}

d1=coverageData[[1]]
avgCovScaffold=lapply(coverageData,function(x){
  mean(x$V4)
})



png(paste0(currentDirectory,"imgs/",sample,"avg.cov.scaffolds.png"), width=3000, height = 800,units = "px")
plot(
  1:length(scaffolds),
  unlist(avgCovScaffold), 
  xlab="Scaffold",ylab="Avg. Coverage",
  main="Avg, Coverage. SampleH09", 
  type='l', col="darkred", axes=F)
axis(1,at=seq(1,length(scaffolds)), labels = scaffolds,las=2, cex=0.5)
axis(2)
dev.off()
