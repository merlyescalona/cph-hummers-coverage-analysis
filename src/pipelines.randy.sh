################################################################################
# Coverteing GFF to BED files
GFF="$HOME/files/Calypte_anna.gene.CDS.2750.gff"
cat $GFF | awk -v OFS='\t' '{print $1,$4,$5}' > $HOME/files/targets.filtered.bed
# filtering 0bp length targets
targetsFilteredBed="$HOME/files/targets.filtered.bed"
targetsBed="$HOME/files/targets.bed"
cat $targetsBed | awk '($3-$2)>0 {print $0}' > $HOME/files/targets.filtered.bed
cat $targetsBed | awk '{print $1}' | uniq > $HOME/files/scaffolds.txt
cat $targetsBed | awk '{print $1}' | uniq -c | awk -v OFS="\t" '{print $2,$1}'> $HOME/files/num.targets.per.scaffold.txt
################################################################################
for bamfile in $(find $HOME/map2Aan/originals -name "*.bam"); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  echo $tag >>$HOME/map2Aan/files/samples.txt
done
for bamfile in $(find $HOME/map2Swift/originals -name "*.bam"); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  echo $tag >>$HOME/map2Swift/files/samples.txt
done
################################################################################################################################################################################################################################################
# Getting targeted and untargeted bam files
################################################################################################################################################################################################################################################

for bamfile in $(find $HOME/map2Aan/originals -name "*.bam" ); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  echo "nohup /home/merly/src/bedtools2/bin/bedtools intersect -abam $bamfile -b $targetsFilteredBed > $HOME/map2Aan/targeted/${tag}.targeted.bam &"
  echo "nohup /home/merly/src/bedtools2/bin/bedtools intersect -v -abam $bamfile -b $targetsFilteredBed > $HOME/map2Aan/untargeted/${tag}.untargeted.bam &"
done
<<MAP2AANintersect
MAP2AANintersect


for bamfile in $(find $HOME/map2Swift/originals -name "*.bam" ); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  echo "nohup /home/merly/src/bedtools2/bin/bedtools intersect -abam $bamfile -b $targetsFilteredBed > $HOME/map2Swift/targeted/${tag}.targeted.bam &"
  echo "nohup /home/merly/src/bedtools2/bin/bedtools intersect -v -abam $bamfile -b $targetsFilteredBed > $HOME/map2Swift/untargeted/${tag}.untargeted.bam &"
done
<<MAP2SWIFTintersect
MAP2SWIFTintersect
################################################################################################################################################################################################################################################
# Coverage por targeted regions
################################################################################################################################################################################################################################################

for bamfile in $(find $HOME/map2Aan/originals -name "*.bam"); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  echo "nohup /home/merly/src/bedtools2/bin/bedtools coverage -hist -abam $bamfile -b $targetsFilteredBed | gzip > $HOME/map2Aan/bedtools/cov/${tag}.cov.gz &"
done
<<BEDTOOLSHISTMAPTOANNA
BEDTOOLSHISTMAPTOANNA
for tag in $(cat $HOME/map2Aan/samples.txt); do
  echo "nohup zcat $HOME/map2Aan/bedtools/cov/$tag.cov.gz | grep -v ^all | gzip > $HOME/map2Aan/bedtools/nohist/$tag.nohist.gz &"
  echo "nohup zcat $HOME/map2Aan/bedtools/cov/$tag.cov.gz | grep ^all | gzip > $HOME/map2Aan/bedtools/hist/$tag.hist.gz &"
done
<<BEDTOOLSHISTsplitMAPTOANNA
BEDTOOLSHISTsplitMAPTOANNA

################################################################################################################################################################################################################################################
# Working now with data mapped to the outgroup
for bamfile in $(find $HOME/map2Swift/originals -name "*.bam" ); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  echo $tag >> samples.txt
  #echo $tag
  echo "nohup /home/merly/src/bedtools2/bin/bedtools intersect -abam $bamfile -b $HOME/files/targets.filtered.bed > $HOME/map2Swift/targeted/${tag}.targeted.bam &"
  echo "nohup /home/merly/src/bedtools2/bin/bedtools intersect -v -abam $bamfile -b $HOME/files/targets.filtered.bed > $HOME/map2Swift/untargeted/${tag}.untargeted.bam &"
done
<<BEDTOOLSINTERSECTSWIFT
BEDTOOLSINTERSECTSWIFT
################################################################################################################################################################################################################################################
# map2Swift bedtools coverage
for bamfile in $(find $HOME/map2Swift/originals -name "*.bam"); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  echo "nohup /home/merly/src/bedtools2/bin/bedtools coverage -hist -abam $bamfile -b $HOME/files/targets.filtered.bed | gzip > $HOME/map2Swift/bedtools/cov/${tag}.cov.gz &"
done

<<BEDTOOLSHISTMAPTOSWIFT
BEDTOOLSHISTMAPTOSWIFT

################################################################################################################################################################################################################################################
# map2Swift hist split
for tag in $(cat $HOME/map2Swift/samples.txt); do
  echo "nohup zcat $HOME/map2Swift/bedtools/cov/$tag.cov.gz | grep -v ^all | gzip > $HOME/map2Swift/bedtools/nohist/$tag.nohist.gz &"
  echo "nohup zcat $HOME/map2Swift/bedtools/cov/$tag.cov.gz | grep ^all | gzip > $HOME/map2Swift/bedtools/hist/$tag.hist.gz &"
done

<<BEDTOOLSHISTMAPTOSWIFTgrep
BEDTOOLSHISTMAPTOSWIFTgrep

################################################################################################################################################################################################################################################
# files for angsd
find $HOME/map2Swift/originals -name "*.bam" >  $HOME/map2Swift/files/map2Swift.originals.bamfiles.txt
find $HOME/map2Swift/targeted -name "*.bam"  > $HOME/map2Swift/files/map2Swift.targeted.bamfiles.txt
find $HOME/map2Swift/untargeted -name "*.bam"  > $HOME/map2Swift/files/map2Swift.untargeted.bamfiles.txt
find $HOME/map2Aan/originals -name "*.bam"  > $HOME/map2Aan/files/map2Aan.originals.bamfiles.txt
find $HOME/map2Aan/targeted -name "*.bam" > $HOME/map2Aan/files/map2Aan.targeted.bamfiles.txt
find $HOME/map2Aan/untargeted -name "*.bam"  > $HOME/map2Aan/files/map2Aan.untargeted.bamfiles.txt

nohup  $HOME/src/angsd/angsd -bam $HOME/map2Swift/files/map2Swift.originals.bamfiles.txt -out $HOME/map2Swift/files/angsd.map2Swift.originals -doCounts 1 -doDepth 1 -maxDepth 1000 &
nohup  $HOME/src/angsd/angsd -bam $HOME/map2Swift/files/map2Swift.targeted.bamfiles.txt -out $HOME/map2Swift/files/angsd.map2Swift.targeted -doCounts 1 -doDepth 1 -maxDepth 1000 &
nohup  $HOME/src/angsd/angsd -bam $HOME/map2Swift/files/map2Swift.untargeted.bamfiles.txt -out $HOME/map2Swift/files/angsd.map2Swift.untargeted -doCounts 1 -doDepth 1 -maxDepth 1000 &
nohup  $HOME/src/angsd/angsd -bam $HOME/map2Aan/files/map2Aan.originals.bamfiles.txt -out $HOME/map2Aan/files/angsd.map2Aan.originals -doCounts 1 -doDepth 1 -maxDepth 1000 &
nohup  $HOME/src/angsd/angsd -bam $HOME/map2Aan/files/map2Aan.targeted.bamfiles.txt -out $HOME/map2Aan/files/angsd.map2Aan.targeted -doCounts 1 -doDepth 1 -maxDepth 1000 &
nohup  $HOME/src/angsd/angsd -bam $HOME/map2Aan/files/map2Aan.untargeted.bamfiles.txt -out $HOME/map2Aan/files/angsd.map2Aan.untargeted -doCounts 1 -doDepth 1 -maxDepth 1000 &


################################################################################################################################################################################################################################################
# Getting the untargeted positions
for bamfile in $(find $HOME/map2Aan/originals -name "*.bam"); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  echo "nohup /home/merly/src/bedtools2/bin/bedtools bamtobed -i $HOME/map2Aan/untargeted/$tag.untargeted.bam > $HOME/map2Aan/bed/$tag.untargeted.bed &"
done
<<BEDMergeMap2Aan
BEDMergeMap2Aan

for bamfile in $(find $HOME/map2Swift/originals -name "*.bam"); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  echo "nohup /home/merly/src/bedtools2/bin/bedtools bamtobed -i $HOME/map2Swift/untargeted/$tag.untargeted.bam > $HOME/map2Swift/bed/$tag.untargeted.bed &"
done
<<BEDMergeMap2Swift
BEDMergeMap2Swift
####################################################################################################################################################
####################################################################################################################################################
for bamfile in $(find $HOME/map2Aan/originals -name "*.bam"); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  echo "nohup sort -k1,1 -k2,2n $HOME/map2Aan/bed/$tag.untargeted.bed > $HOME/map2Aan/bed/$tag.untargeted.sorted.bed &"
done
<<BEDMergeMap2Aan
BEDMergeMap2Aan

for bamfile in $(find $HOME/map2Swift/originals -name "*.bam"); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  echo "nohup sort -k1,1 -k2,2n $HOME/map2Swift/bed/$tag.untargeted.bed > $HOME/map2Swift/bed/$tag.untargeted.sorted.bed &"
done
<<BEDMergeMap2Swift
BEDMergeMap2Swift

for bamfile in $(find $HOME/map2Aan/originals -name "*.bam"); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  echo "nohup bedtools merge -i $HOME/map2Aan/bed/$tag.untargeted.sorted.bed > $HOME/map2Aan/bed/$tag.untargeted.merged.bed &"
done
<<BEDMergeMap2Aan
BEDMergeMap2Aan

for bamfile in $(find $HOME/map2Swift/originals -name "*.bam"); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  echo "nohup bedtools merge -i $HOME/map2Swift/bed/$tag.untargeted.sorted.bed > $HOME/map2Swift/bed/$tag.untargeted.merged.bed &"
done
<<BEDMergeMap2Swidt
BEDMergeMap2Swidt


####################################################################################################################################
for bamfile in $(find $HOME/map2Swift/originals -name "*.bam"); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  echo "nohup /home/merly/src/bedtools2/bin/bedtools bamtobed -i $HOME/map2Swift/untargeted/$tag.untargeted.bam > $HOME/map2Swift/bed/$tag.untargeted.bed &"
done
<<BEDbamtobedMap2Swift
BEDbamtobedMap2Swift
####################################################################################################################################
cat $HOME/map2Aan/bed/*.untargeted.merged.bed > /home/merly/map2Aan/files/untargeted.all.bed
nohup /home/merly/src/bedtools2/bin/bedtools merge -i /home/merly/map2Aan/files/untargeted.all.bed > /home/merly/map2Aan/files/untargeted.merged.bed &
cat $HOME/map2Swift/bed/*.untargeted.merged.bed > /home/merly/map2Swift/files/untargeted.all.bed
nohup /home/merly/src/bedtools2/bin/bedtools merge -i /home/merly/map2Swift/files/untargeted.all.bed > /home/merly/map2Swift/files/untargeted.merged.bed &
####################################################################################################################################
## REMAPPING
################################################################################################################################################################################################################################################
# map2Swift hist split
mkdir map2Aan/cov2 map2Aan/his2 map2Aan/nohist2 map2Swift/cov2 map2Swift/his2 map2Swift/nohist2
cat $HOME/map2Aan/files/untargeted.merged.bed | grep -v ^Error > $HOME/map2Aan/files/untargeted.bed
cat $HOME/map2Swift/files/untargeted.merged.bed | grep -v ^Error > $HOME/map2Swift/files/untargeted.bed

for bamfile in $(find $HOME/map2Aan/originals -name "*.bam"); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  echo "nohup /home/merly/src/bedtools2/bin/bedtools coverage -hist -abam $bamfile -b $HOME/map2Aan/files/untargeted.bed | gzip > $HOME/map2Aan/bedtools/cov2/${tag}.cov.gz &"
done
<<CYCLEDCOVMAP2AAN
CYCLEDCOVMAP2AAN
for bamfile in $(find $HOME/map2Aan/originals -name "*.bam"); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  echo "nohup zcat $HOME/map2Aan/bedtools/cov2/$tag.cov.gz | grep -v ^all | gzip > $HOME/map2Aan/bedtools/nohist2/$tag.nohist.gz &"
  echo "nohup zcat $HOME/map2Aan/bedtools/cov2/$tag.cov.gz | grep ^all | gzip > $HOME/map2Aan/bedtools/hist2/$tag.hist.gz &"
done
<<CYCLED2COVMAP2AAN
CYCLED2COVMAP2AAN
for bamfile in $(find $HOME/map2Swift/originals -name "*.bam"); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  echo "nohup /home/merly/src/bedtools2/bin/bedtools coverage -hist -abam $bamfile -b $HOME/files/untargeted.bed | gzip > $HOME/map2Swift/bedtools/cov2/${tag}.cov.gz &"
done
<<CYCLEDCOVMAP2swift
CYCLEDCOVMAP2swift
for bamfile in $(find $HOME/map2Swift/originals -name "*.bam"); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  echo "nohup zcat $HOME/map2Swift/bedtools/cov2/$tag.cov.gz | grep -v ^all | gzip > $HOME/map2Swift/bedtools/nohist2/$tag.nohist.gz &"
  echo "nohup zcat $HOME/map2Swift/bedtools/cov2/$tag.cov.gz | grep ^all | gzip > $HOME/map2Swift/bedtools/hist2/$tag.hist.gz &"
done
<<CYCLED2COVMAP2swift
CYCLED2COVMAP2swift



# merging both offtarget files (map2aan and map2Swift)
cat /home/merly/map2Aan/files/untargeted.bed /home/merly/map2Swift/files/untargeted.bed > /home/merly/files/untargeted.all.1.bed
nohup /home/merly/src/bedtools2/bin/bedtools merge -i /home/merly/files/untargeted.all.1.bed | grep -v ^Error > /home/merly/files/untargeted.bed &


# Half of the things i did were incorrect.
# tusns out original target files does not match the map2Swift dataset.
# Chaetura_pelagica.CDS.gff
# Rscript target.coversion.R
# This files extracts the targets from the 48birds_ortholog.list.chi.anna.cpe.hum.finch and the Swift gff
