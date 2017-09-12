# Generation of genomewide scaffolds
WD="/media/merly/Baymax/research/cph-visit/coverage-analysis/"
# WD="/home/merly/"
allgenes=paste0(WD,"files/Calypte_anna.gene.CDS.gff")
anna=read.table(allgenes, colClasses = c("character","character","character",
                                         "numeric","numeric","character",
                                         "character","character","character"))
swiftfull=paste0(WD,"files/Chaetura_pelagica.CDS.gff")
swift=read.table(swiftfull, colClasses = c("character","character","character",
                                         "numeric","numeric","character",
                                         "character","character","character"))
sAnnaAll=split(anna,anna$V1)
sSwiftAll=split(swift,swift$V1)
maxPositionPerScaffoldAnna=sapply(sAnnaAll,function(x){max(x$V5)})
maxPositionPerScaffoldSwift=sapply(sSwiftAll,function(x){max(x$V5)})

maxLen=max(length(maxPositionPerScaffoldAnna),length(maxPositionPerScaffoldSwift))
unionScaffolds=union(names(maxPositionPerScaffoldAnna), names(maxPositionPerScaffoldSwift))
df=data.frame(anna=rep(0,length(unionScaffolds)), swift=rep(0,length(unionScaffolds)))
rownames(df)=unionScaffolds

df[names(maxPositionPerScaffoldAnna),]$anna=maxPositionPerScaffoldAnna[names(maxPositionPerScaffoldAnna)]
df[names(maxPositionPerScaffoldSwift),]$swift=maxPositionPerScaffoldSwift[names(maxPositionPerScaffoldSwift)]

dfAnna=df$anna[df$anna!=0]
names(dfAnna)=rownames(df)[df$anna!=0]
finalBed=cbind(names(dfAnna),rep(1,length(dfAnna)),dfAnna)
write.table(finalBed,paste0(WD,"files/scaffolds.anna.bed"), row.names = F,col.names = F,quote = F, sep="\t")

dfSwift=df$swift[df$swift!=0]
names(dfSwift)=rownames(df)[df$swift!=0]
finalBed=cbind(names(dfSwift),rep(1,length(dfSwift)),dfSwift)
write.table(finalBed,paste0(WD,"files/scaffolds.swift.bed"), row.names = F,col.names = F,quote = F, sep="\t")
