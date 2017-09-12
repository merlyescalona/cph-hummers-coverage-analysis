# WD="/home/merly/files/"
WD="/media/merly/Baymax/research/cph-visit/coverage-analysis/files/"
# WD="/home/merly/files/"
birds48filename=paste0(WD,"48birds_ortholog.list.chi.anna.cpe.hum.finch")
swiftgff=paste0(WD,"Chaetura_pelagica.CDS.gff")
annagff=paste0(WD,"Calypte_anna.gene.CDS.2750.gff")
birds=read.table(birds48filename, header=T,
    colClasses=c("character","numeric","character","character",
      "character","character","character"))
swift=read.table(swiftgff,
    colClasses=c("character","character","character","numeric",
      "numeric","character","character","character","character"))
anna=read.table(annagff,
    colClasses=c("character","character","character","numeric",
      "numeric","character","character","character","character"))

birds.2=birds[birds$anna!="-",]
birds.3=birds.2[birds.2$swift!="-",]
birds.4=unique(birds.3)
birds.5=birds.4[,c("anna","swift")]
birds.5$swift=paste0("Parent=",birds.5$swift,";")
birds.5$anna=paste0("Parent=Aan_", birds.5$anna,";")
annaGenes=unique(anna$V9[birds.5$anna %in% anna$V9])
birds.6=birds.5[birds.5$anna %in% annaGenes, ]
swiftGenes=birds.6$swift

swift.2=swift[swift$V3=="CDS",]
swift.3=swift.2[swift.2$V9 %in% birds.6$swift,]
swift.3$size=swift.3$V5-swift.3$V4
swift.4=swift.3[swift.3$size>0,]
anna.2=anna[(anna$V5-anna$V4)>0,]
swift.4=swift.4[,1:(ncol(swift.4)-1)]
anna.3=anna.2[anna.2$V9 %in% unique(birds.6$anna),]


write.table(swift.4[,c(1,4,5)],
  file=paste0(WD,"targets.swift.2.bed"),
  row.names = F,col.names = F, quote = F, sep = "\t")
write.table(swift.4,
  file=paste0(WD,"Chaetura_pelagica.CDS.2.gff"),
  row.names = F,col.names = F, quote = F, sep = "\t")
write.table(anna.3[,c(1,4,5)],
  file=paste0(WD,"targets.anna.3.bed"),
  row.names = F,col.names = F, quote = F, sep = "\t")
write.table(anna.3,
  file=paste0(WD,"Calypte_anna.gene.CDS.2750.3.gff"),
  row.names = F,col.names = F, quote = F, sep = "\t")
write.table(birds.6,
  file=paste0(WD,"Calypte_anna.Chaetura_pelagica.ids"),
  row.names = F,col.names = F, quote = F, sep = "\t")
