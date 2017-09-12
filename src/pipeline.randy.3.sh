################################################################################################################################################################
# Coverteing GFF to BED files
gffAnna="$HOME/files/Calypte_anna.gene.CDS.2750.3.gff"
gffSwift="$HOME/files/Chaetura_pelagica.CDS.2.gff"
# filtering 0bp length targets
# --------------------------------------------------------------------------------------------------------------------------------------------------------------
targetSwift="$HOME/files/targets.swift.2.bed"
targetAnna="$HOME/files/targets.anna.3.bed"
originalsSwift="$HOME/swift/originals"
originalsAnna="$HOME/anna/originals"
offbedanna="$HOME/files/scaffolds.anna.bed"
offbedswift="$HOME/files/scaffolds.swift.bed"
################################################################################################################################################################
mkdir -p anna/files2 anna/bedtools2/cov anna/bedtools2/hist anna/bedtools2/nohist anna/img2
mkdir -p swift/files2 swift/bedtools2/cov swift/bedtools2/hist swift/bedtools2/nohist swift/img2
mkdir -p anna/bedtools2/cov2 anna/bedtools2/hist2 anna/bedtools2/nohist2
mkdir -p swift/bedtools2/cov2 swift/bedtools2/hist2 swift/bedtools2/nohist2
################################################################################################################################################################
# Off-target regions
# generated a total bed with possible position from all the scaffolds in the GFFs.
# I cant filter both files because they are different genoms and positions do not match
nohup /home/merly/src/bedtools2/bin/bedtools subtract -a $offbedanna -b $gffAnna > "$HOME/files/offtarget.anna.bed" &
nohup /home/merly/src/bedtools2/bin/bedtools subtract -a $offbedswift -b $gffSwift > "$HOME/files/offtarget.swift.bed" &
################################################################################################################################################################
for bamfile in $(find $originalsAnna -name "*.bam"); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  echo $tag >>$HOME/anna/files/samples.txton
done
for bamfile in $(find $originalsSwift -name "*.bam"); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  echo $tag >>$HOME/swift/files/samples.txt
done

# swift bedtools coverage
for bamfile in $(find $originalsAnna -name "*.bam" |  tail -n+31   ); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  nohup /home/merly/src/bedtools2/bin/bedtools coverage -hist -abam $bamfile -b $gffAnna| gzip > $HOME/anna/bedtools2/cov/${tag}.cov.gz &
done

for bamfile in $(find $originalsSwift -name "*.bam" | tail -n+32 | head -16 ); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  nohup /home/merly/src/bedtools2/bin/bedtools coverage -hist -abam $bamfile -b $gffSwift | gzip > $HOME/swift/bedtools2/cov/${tag}.cov.gz &
done
################################################################################################################################################################################################################################################
#  hist split
for tag in $(cat $HOME/anna/files/samples.txt  ); do
  nohup zcat $HOME/anna/bedtools2/cov/$tag.cov.gz | grep -v ^all | gzip > $HOME/anna/bedtools2/nohist/$tag.nohist.gz &
  nohup zcat $HOME/anna/bedtools2/cov/$tag.cov.gz | grep ^all | gzip > $HOME/anna/bedtools2/hist/$tag.hist.gz &
  sleep 3
done
for tag in $(cat $HOME/swift/files/samples.txt  ); do
  echo $tag
  nohup zcat $HOME/swift/bedtools2/cov/$tag.cov.gz | grep -v ^all | gzip > $HOME/swift/bedtools2/nohist/$tag.nohist.gz &
  nohup zcat $HOME/swift/bedtools2/cov/$tag.cov.gz | grep ^all | gzip > $HOME/swift/bedtools2/hist/$tag.hist.gz &
  sleep 1
done
################################################################################################################################################################################################################################################
# Getting targeted and untargeted bam files
################################################################################################################################################################################################################################################
ls -lah $HOME/anna/offtarget2/
for bamfile in $(find $originalsAnna -name "*.bam" |  tail -n+41 | head -10); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  # nohup /home/merly/src/bedtools2/bin/bedtools intersect -abam $bamfile -b $targetAnna > $HOME/anna/ontarget/${tag}.ontarget.bam &
  nohup /home/merly/src/bedtools2/bin/bedtools intersect -abam $bamfile -b "$HOME/files/offtarget.anna.bed" > $HOME/anna/offtarget/${tag}.offtarget.bam &
done
ls $HOME/swift/offtarget2/
for bamfile in $(find $originalsSwift -name "*.bam" | head -1); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  # nohup /home/merly/src/bedtools2/bin/bedtools intersect -abam $bamfile -b $targetSwift > $HOME/swift/ontarget/${tag}.ontarget.bam &
  nohup /home/merly/src/bedtools2/bin/bedtools intersect -abam $bamfile -b "$HOME/files/offtarget.swift.bed" > $HOME/swift/offtarget/${tag}.offtarget.bam &
done


for bamfile in $(find $originalsAnna -name "*.bam" |  tail -n+41 | head -20); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  nohup /home/merly/src/bedtools2/bin/bedtools coverage -hist -abam $bamfile -b "$HOME/files/offtarget.anna.bed" | gzip > $HOME/anna/bedtools2/cov/${tag}.cov.gz &
done
# next 32
for bamfile in $(find $originalsSwift -name "*.bam" | tail -n+42   | head -6); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  nohup /home/merly/src/bedtools2/bin/bedtools coverage -hist -abam $bamfile -b "$HOME/files/offtarget.swift.bed" | gzip > $HOME/swift/bedtools2/cov/${tag}.cov.gz &
done
################################################################################################################################################################################################################################################
#  hist split
for tag in $(cat $HOME/anna/files/samples.txt  ); do
  echo $tag
  nohup zcat $HOME/anna/bedtools2/cov/$tag.cov.gz | grep -v ^all | gzip > $HOME/anna/bedtools2/nohist/$tag.nohist.gz &
  nohup zcat $HOME/anna/bedtools2/cov/$tag.cov.gz | grep ^all | gzip > $HOME/anna/bedtools2/hist/$tag.hist.gz &
  sleep 1
done
for tag in $(cat $HOME/swift/files/samples.txt   ); do
  echo $tag
  nohup zcat $HOME/swift/bedtools2/cov/$tag.cov.gz | grep -v ^all | gzip > $HOME/swift/bedtools2/nohist/$tag.nohist.gz &
  nohup zcat $HOME/swift/bedtools2/cov/$tag.cov.gz | grep ^all | gzip > $HOME/swift/bedtools2/hist/$tag.hist.gz &
  sleep 1
done





# Half of the things i did were incorrect.
# tusns out original target files does not match the swift dataset.
# Chaetura_pelagica.CDS.gff
# Rscript target.coversion.R
# This files extracts the targets from the 48birds_ortholog.list.chi.anna.cpe.hum.finch and the Swift gff

#################################################################################################################################################################################
# files for angsd
#################################################################################################################################################################################
find $HOME/swift/originals -name "*.bam"  > $HOME/swift/files/swift.originals.bamfiles.txt
find $HOME/swift/ontarget -name "*.bam"  > $HOME/swift/files/swift.ontarget.bamfiles.txt
find $HOME/swift/offtarget -name "*.bam"  > $HOME/swift/files/swift.offtarget.bamfiles.txt
find $HOME/swift/offtarget2 -name "*.bam"  > $HOME/swift/files/swift.offtarget2.bamfiles.txt
find $HOME/anna/originals -name "*.bam"  > $HOME/anna/files/anna.originals.bamfiles.txt
find $HOME/anna/ontarget -name "*.bam" > $HOME/anna/files/anna.ontarget.bamfiles.txt
find $HOME/anna/offtarget -name "*.bam"  > $HOME/anna/files/anna.offtarget.bamfiles.txt
find $HOME/anna/offtarget2 -name "*.bam"  > $HOME/anna/files/anna.offtarget2.bamfiles.txt

# nohup  $HOME/src/angsd/angsd -bam  $HOME/swift/files/swift.originals.bamfiles.txt -out $HOME/swift/files/swift.originals.angsd. -doCounts 1 -doDepth 1 -maxDepth 1000 &
nohup  $HOME/src/angsd/angsd -bam  $HOME/swift/files/swift.ontarget.bamfiles.txt -out $HOME/swift/files2/swift.ontarget.angsd -doCounts 1 -doDepth 1 -maxDepth 1000 &
# nohup  $HOME/src/angsd/angsd -bam  $HOME/swift/files/swift.offtarget.bamfiles.txt -out $HOME/swift/files/swift.offtarget.angsd. -doCounts 1 -doDepth 1 -maxDepth 1000 &
nohup  $HOME/src/angsd/angsd -bam  $HOME/swift/files/swift.offtarget2.bamfiles.txt -out $HOME/swift/files2/swift.offtarget2.angsd -doCounts 1 -doDepth 1 -maxDepth 1000 &
# nohup  $HOME/src/angsd/angsd -bam  $HOME/anna/files/anna.originals.bamfiles.txt -out $HOME/anna/files/anna.originals.angsd. -doCounts 1 -doDepth 1 -maxDepth 1000 &
nohup  $HOME/src/angsd/angsd -bam  $HOME/anna/files/anna.ontarget.bamfiles.txt -out $HOME/anna/files2/anna.ontarget.angsd -doCounts 1 -doDepth 1 -maxDepth 1000 &
# nohup  $HOME/src/angsd/angsd -bam  $HOME/anna/files/anna.offtarget.bamfiles.txt -out $HOME/anna/files/anna.offtarget.angsd. -doCounts 1 -doDepth 1 -maxDepth 1000 &
nohup  $HOME/src/angsd/angsd -bam  $HOME/anna/files/anna.offtarget2.bamfiles.txt -out $HOME/anna/files2/anna.offtarget2.angsd -doCounts 1 -doDepth 1 -maxDepth 1000 &



#################################################################################################################################################################################
#################################################################################################################################################################################
#################################################################################################################################################################################
#################################################################################################################################################################################
cat Calypte_anna.gene.CDS.2750.3.gff | cut -f9 | uniq | tr "_" " " | awk '{print $2}' | tr ";" " " > g.trees.to.filter.txt
WD="/media/merly/Baymax/research/cph-visit/coverage-analysis"
for gtreefile in $(cat $WD/files/g.trees.to.filter.txt); do
  cp "$WD/trees/gtrees/$gtreefile" "$WD/trees/gtrees.2/"
done
#################################################################################################################################################################################
# Gettign numbre of possible reads per sample
#################################################################################################################################################################################
for item in $(find /home/fonseca/WORK/HUMMERS/1-CLEANmerged -name "*.fastq.gz"); do
  tag=$(basename $item .fastq.gz | tr "_" " " | awk '{print $3}')
  sample=$(basename $item .fastq.gz | tr "_" " " | awk '{print $2}')
  numreads=$(wc -l $item | awk '{print $1}')
  echo -e "$sample\t$tag"
  echo -e "$sample\t$tag\t$numreads" >> $HOME/files/raw.numreads.txt
done
#################################################################################################################################################################################
# getting number of mapped reads
#################################################################################################################################################################################
<<MAPPEDREADS
offtargetAnna="$HOME/anna/offtarget"
for bamfile in $(find $offtargetAnna -name "*.bam"); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  numreads=$(samtools view -c -F 4 $bamfile)
  echo -e "$tag"
  echo -e "$tag\t$numreads" >> $HOME/files/mapped.anna.numreads.offtarget.txt
done

offtargetSwift="$HOME/swift/offtarget"
for bamfile in $(find $offtargetSwift -name "*.bam"); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  numreads=$(samtools view -c -F 4 $bamfile)
  echo -e "$tag"
  echo -e "$tag\t$numreads" >> $HOME/files/mapped.swift.numreads.offtarget.txt
done
MAPPEDREADS

nohup bash mapped.anna.reads.sh &
nohup bash mapped.swift.reads.sh &
