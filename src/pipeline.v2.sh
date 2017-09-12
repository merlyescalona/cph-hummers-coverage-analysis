cd# Getting coverage per BAM file
# both ANGSD and bedtools
# 1. all
# 2. filtered by targets
# 3. filtered by untargeted
WD="/media/merly/Luigi/research/cph-visit/coverage-analysis"
out="/media/merly/Baymax/research/cph-visit/"
targets="$WD/files/Calypte_anna.gene.CDS.2750.gff"
targetbed="/media/merly/Luigi/research/cph-visit/coverage-analysis/files/targets.bed"
cd $WD
mkdir -p targeted untargeted
mkdir -p genomecov/all genomecov/targeted genomecov/untargeted
mkdir -p samtools/all samtools/targeted samtools/untargeted
mkdir -p angsd bed bedtools targetZero

# indexing bam files
for bamfile in $(find $WD/originals -name "*.bam"); do
  echo $tag
  samtools index $bamfile
done

# Transform GFF file to BED
# cols: SCAFFOLD, START,END,SIZE
zcat $targets | awk -v OFS='\t' '{print $1,$4,$5, $5-$4}' > $WD/files/targets.bed
zcat /home/merly/files/Calypte_anna.gene.CDS.2750.gff | awk -v OFS='\t' '{print $1,$4,$5}' > $WD/files/targets.bed

# Generate intersected bam
# BEDTOOLS - ALL - TARGETED
# dowloading H22.
for bamfile in $(find $(pwd)/originals -name "*.bam" ); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  echo $tag
  bedtools intersect -abam $bamfile -b $targetbed > targeted/${tag}.targeted.bam &
  bedtools intersect -v -abam $bamfile -b $targetbed > untargeted/${tag}.untargeted.bam
done

# BEDTOOLS - ALL - TARGETED
for bamfile in $(find $WD/originals -name "*.bam"); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  bedtools genomecov -ibam $bamfile > genomecov/all/${tag}.all.cov &
  bedtools genomecov -ibam targeted/${tag}.targeted.bam > genomecov/targeted/${tag}.targeted.cov &
  bedtools genomecov -ibam untargeted/${tag}.untargeted.bam > genomecov/untargeted/${tagqq}.untargeted.cov
done

# Forgot to say to genome cov that I wanted the indices to start in 1
# so I'll take now the files and add 1 to the position it represents
for bamfile in $(find $WD/originals -name "*.bam"); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  mv genomecov/all/$tag.all.cov genomecov/all/$tag.all.cov.bak
  mv genomecov/targeted/$tag.targeted.cov genomecov/targeted/$tag.targeted.cov.bak
  mv genomecov/untargeted/$tag.untargeted.cov genomecov/untargeted/$tag.untargeted.cov.bak
  cat genomecov/all/$tag.all.cov.bak | awk -v OFS='\t'  '{ pos=$2; newpos=pos+1; print $1,newpos,$3,$4,$5 }'  > genomecov/all/$tag.all.cov
  cat genomecov/targeted/$tag.targeted.cov.bak | awk -v OFS='\t'  '{ pos=$2; newpos=pos+1; print $1,newpos,$3,$4,$5 }' > genomecov/targeted/$tag.targeted.cov
  cat genomecov/untargeted/$tag.untargeted.cov.bak | awk -v OFS='\t'  '{ pos=$2; newpos=pos+1; print $1,newpos,$3,$4,$5 }' > genomecov/untargeted/$tag.untargeted.cov
done

for bamfile in $(find $WD/originals -name "*.bam"); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  echo -e "\r$tag - bedtools    "
  bedtools coverage -hist -abam $bamfile -b $targets | gzip > "$out/bedtools/cov/${tag}.cov.gz"
  echo -e "\r$tag - nohist      "
  zcat $out/bedtools/cov/$tag.cov.gz | grep -v ^all | gzip > "$out/bedtools/nohist/${tag}.no.hist.cov.gz"
  echo -e "\r$tag - hist        "
  zcat $out/bedtools/cov/$tag.cov.gz| grep ^all | gzip > "$out/bedtools/hist/${tag}.hist.cov.gz"
done

# Checking targeted regions where all bases have 0 coverage
for bamfile in $(find $WD/originals -name "*.bam"); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  bedtools coverage -hist -abam $bamfile -b $targets  -hist | gzip > "$out/bedtools/${tag}.cov.gz" &
  echo -e "\r$tag                                 "
  samtools depth $bamfile | gzip > $out/samtools/all/$tag.depth.txt.gz &
  echo -e "\r$tag-targeted"
  samtools depth -a targeted/$tag.targeted.bam | gzip > $out/samtools/targeted/$tag.depth.txt.gz &
  echo -e "\r$tag-untargeted"
  samtools depth -a untargeted/$tag.untargeted.bam | gzip > $out/samtools/untargeted/$tag.depth.txt.gz
done

<<comments
Until here, I have coverage per positions from 3 files: all, targeted and untargeted
I'll need to get:
3. get size of the untargeted regions
    - bam to bed of untargeted:
    so i'll be able to get the size of the regions and calculate the coverage
2. extract coverage of the untargeted regions
3. extract coverage of the untargeted regions per sample
4. extract coverage of the targeted regions
5. extract coverage of the targeted regions per sample
6. extract coverage per sample - targeted and captured
7. check those regions that are not captured by any sample:
    - get new covered size (number of bases covered with the targets)
    - extract all the coverage related information again
comments

# Extract untargeted region size
for bamfile in $(find $HOME/map2Aan/originals -name "*.bam"); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  echo "nohup bedtools bamtobed -i map2Aan/untargeted/$tag.untargeted.bam > bed/$tag.untargeted.bed &"
  echo "nohup sort -k1,1 -k2,2n map2Aan/bed/$tag.untargeted.bed > bed/$tag.untargeted.sorted.bed &"º
  echo "nohup bedtools merge -i map2Aan/bed/$tag.untargeted.sorted.bed > bed/$tag.untargeted.merged.bed &"
done

for bamfile in $(find $HOME/map2Swift/originals -name "*.bam"); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  echo "nohup bedtools bamtobed -i map2Swift/untargeted/$tag.untargeted.bam > bed/$tag.untargeted.bed &"
  echo "nohup sort -k1,1 -k2,2n map2Swift/bed/$tag.untargeted.bed > bed/$tag.untargeted.sorted.bed &"º
  echo "nohup bedtools merge -i map2Swift/bed/$tag.untargeted.sorted.bed > bed/$tag.untargeted.merged.bed &"
done

# merge all bedfiles to have a common sense of all the untargeted regions
bedfiles=""
for bamfile in $(find $WD/originals -name "*.bam"); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  bedtools merge -i bed/$tag.untargeted.merged.bed > bed/$tag.untargeted.merged.sorted.bed
  bedfiles="$bedfiles -i bed/$tag.untargeted.merged.sorted.bed"
done

bedtools merge $bedfiles > bed/untargeted.all.merged.bed
bedtools merge bed/untargeted.all.merged.bed > bed/untargeted.all.merged.sorted.bed




<<TOREAD
https://davetang.org/muse/2013/09/07/creating-a-coverage-plot-in-r/
TOREAD

# Getting all file paths for all bam files
find $WD/originals/ -name "*.bam"  > $WD/files/bamfiles.originals.txt
find $WD/targeted-captured/ -name "*.bam"  > $WD/files/bamfiles.targeted.txt
find $WD/untargeted-captured/ -name "*.bam"  > $WD/files/bamfiles.untargeted.txt

angsd -bam $WD/files/bamfiles.originals.txt -out files/originals -doCounts 1 -doDepth 1 -maxDepth 1000
angsd -bam $WD/files/bamfiles.targeted.txt -out files/targeted -doCounts 1 -doDepth 1 -maxDepth 1000
angsd -bam $WD/files/bamfiles.untargeted.txt -out files/untargeted -doCounts 1 -doDepth 1 -maxDepth 1000
