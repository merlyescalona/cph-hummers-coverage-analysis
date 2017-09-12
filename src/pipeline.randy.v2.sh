################################################################################################################################################################
# Coverteing GFF to BED files
gffAnna="$HOME/files/Calypte_anna.gene.CDS.2750.2.gff"
gffSwift="$HOME/files/Chaetura_pelagica.CDS.2.gff"
# filtering 0bp length targets
# --------------------------------------------------------------------------------------------------------------------------------------------------------------
targetSwift="$HOME/files/targets.swift.2.bed"
targetAnna="$HOME/files/targets.anna.2.bed"
originalsSwift="$HOME/swift/originals"
originalsAnna="$HOME/anna/originals"
offbedanna="$HOME/files/scaffolds.anna.bed"
offbedswift="$HOME/files/scaffolds.swift.bed"
################################################################################################################################################################
mkdir -p anna/untargeted anna/files anna/targeted anna/bed anna/bedtools/cov anna/bedtools/hist anna/bedtools/nohist anna/img
mkdir -p swift/untargeted swift/files swift/targeted swift/bed swift/bedtools/cov swift/bedtools/hist swift/bedtools/nohist swift/img
mkdir -p anna/bedtools/cov2 anna/bedtools/hist2 anna/bedtools/nohist2 anna/bedtools/cov3 anna/bedtools/hist3 anna/bedtools/nohist3
mkdir -p swift/bedtools/cov2 swift/bedtools/hist2 swift/bedtools/nohist2 swift/bedtools/cov3 swift/bedtools/hist3 swift/bedtools/nohist3
################################################################################################################################################################
# Off-target regions
# generated a total bed with possible position from all the scaffolds in the GFFs.
# I cant filter both files because they are different genoms and positions do not match
nohup /home/merly/src/bedtools2/bin/bedtools subtract -a $offbedanna -b $gffAnna > "$HOME/files/offtarget.anna.bed" &> nohup.bedtools.substract.anna.out &
nohup /home/merly/src/bedtools2/bin/bedtools subtract -a $offbedswift -b $gffSwift > "$HOME/files/offtarget.swift.bed" &> nohup.bedtools.substract.sw.out &
################################################################################################################################################################
for bamfile in $(find $originalsAnna -name "*.bam"); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  echo $tag >>$HOME/anna/files/samples.txt
done
for bamfile in $(find $originalsSwift -name "*.bam"); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  echo $tag >>$HOME/swift/files/samples.txt
done

# swift bedtools coverage
for bamfile in $(find $originalsAnna -name "*.bam" | tail -n+41  ); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  nohup /home/merly/src/bedtools2/bin/bedtools coverage -hist -abam $bamfile -b $gffAnna| gzip > $HOME/anna/bedtools/cov/${tag}.cov.gz &
done
for bamfile in $(find $originalsSwift -name "*.bam" ); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  nohup /home/merly/src/bedtools2/bin/bedtools coverage -hist -abam $bamfile -b $gffSwift | gzip > $HOME/swift/bedtools/cov/${tag}.cov.gz &
done
################################################################################################################################################################################################################################################
#  hist split
for tag in $(cat $HOME/anna/files/samples.txt  ); do
  nohup zcat $HOME/anna/bedtools/cov/$tag.cov.gz | grep -v ^all | gzip > $HOME/anna/bedtools/nohist/$tag.nohist.gz &
  nohup zcat $HOME/anna/bedtools/cov/$tag.cov.gz | grep ^all | gzip > $HOME/anna/bedtools/hist/$tag.hist.gz &
done
for tag in $(cat $HOME/swift/files/samples.txt  ); do
  nohup zcat $HOME/swift/bedtools/cov/$tag.cov.gz | grep -v ^all | gzip > $HOME/swift/bedtools/nohist/$tag.nohist.gz &
  nohup zcat $HOME/swift/bedtools/cov/$tag.cov.gz | grep ^all | gzip > $HOME/swift/bedtools/hist/$tag.hist.gz &
done
################################################################################################################################################################################################################################################
# Getting targeted and untargeted bam files
################################################################################################################################################################################################################################################
ls $HOME/anna/offtarget2/
for bamfile in $(find $originalsAnna -name "*.bam" |  tail -n+41 | head ); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  # nohup /home/merly/src/bedtools2/bin/bedtools intersect -abam $bamfile -b $targetAnna > $HOME/anna/ontarget/${tag}.ontarget.bam &
  nohup /home/merly/src/bedtools2/bin/bedtools intersect -abam $bamfile -b "$HOME/files/offtarget.anna.bed" > $HOME/anna/offtarget2/${tag}.offtarget.bam &
done
for bamfile in $(find $originalsSwift -name "*.bam" | tail -n+42 | head -5 ); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  # nohup /home/merly/src/bedtools2/bin/bedtools intersect -abam $bamfile -b $targetSwift > $HOME/swift/ontarget/${tag}.ontarget.bam &
  nohup /home/merly/src/bedtools2/bin/bedtools intersect -abam $bamfile -b "$HOME/files/offtarget.swift.bed" > $HOME/swift/offtarget2/${tag}.offtarget.bam &
done


for bamfile in $(find $originalsAnna -name "*.bam" |  tail -n+41 ); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  nohup /home/merly/src/bedtools2/bin/bedtools coverage -hist -abam $bamfile -b "$HOME/files/offtarget.anna.bed" | gzip > $HOME/anna/bedtools/cov2/${tag}.cov.gz &
done
# next 32
for bamfile in $(find $originalsSwift -name "*.bam"  | tail -n+42 | head -6); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  nohup /home/merly/src/bedtools2/bin/bedtools coverage -hist -abam $bamfile -b "$HOME/files/offtarget.swift.bed" | gzip > $HOME/swift/bedtools/cov2/${tag}.cov.gz &
done
################################################################################################################################################################################################################################################
#  hist split
for tag in $(cat $HOME/anna/files/samples.txt  ); do
  nohup zcat $HOME/anna/bedtools/cov2/$tag.cov.gz | grep -v ^all | gzip > $HOME/anna/bedtools/nohist2/$tag.nohist.gz &
  nohup zcat $HOME/anna/bedtools/cov2/$tag.cov.gz | grep ^all | gzip > $HOME/anna/bedtools/hist2/$tag.hist.gz &
done
for tag in $(cat $HOME/swift/files/samples.txt  ); do
   zcat $HOME/swift/bedtools/cov2/$tag.cov.gz | grep -v ^all | gzip > $HOME/swift/bedtools/nohist2/$tag.nohist.gz
   zcat $HOME/swift/bedtools/cov2/$tag.cov.gz | grep ^all | gzip > $HOME/swift/bedtools/hist2/$tag.hist.gz
done



################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################
# New off targets

for bamfile in $(find $originalsAnna -name "*.bam" |  tail -n+41 | head); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  nohup /home/merly/src/bedtools2/bin/bedtools coverage -hist -abam $bamfile -b "$HOME/files/offtarget.anna.bed" | gzip > $HOME/anna/bedtools/cov3/${tag}.cov.gz &
done
# next 32
for bamfile in $(find $originalsSwift -name "*.bam"  | tail -n+42 | head -6 ); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  nohup /home/merly/src/bedtools2/bin/bedtools coverage -hist -abam $bamfile -b "$HOME/files/offtarget.swift.bed" | gzip > $HOME/swift/bedtools/cov3/${tag}.cov.gz &
done
for tag in $(cat $HOME/anna/files/samples.txt  ); do
  nohup zcat $HOME/anna/bedtools/cov3/$tag.cov.gz | grep -v ^all | gzip > $HOME/anna/bedtools/nohist3/$tag.nohist.gz &
  nohup zcat $HOME/anna/bedtools/cov3/$tag.cov.gz | grep ^all | gzip > $HOME/anna/bedtools/hist3/$tag.hist.gz &
  sleep 2
done
for tag in $(cat $HOME/swift/files/samples.txt  ); do
  nohup zcat $HOME/swift/bedtools/cov3/$tag.cov.gz | grep ^all | gzip > $HOME/swift/bedtools/hist3/$tag.hist.gz &
  nohup zcat $HOME/swift/bedtools/cov3/$tag.cov.gz | grep -v ^all | gzip > $HOME/swift/bedtools/nohist3/$tag.nohist.gz &
  sleep 5
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
# nohup  $HOME/src/angsd/angsd -bam  $HOME/swift/files/swift.ontarget.bamfiles.txt -out $HOME/swift/files/swift.ontarget.angsd. -doCounts 1 -doDepth 1 -maxDepth 1000 &
# nohup  $HOME/src/angsd/angsd -bam  $HOME/swift/files/swift.offtarget.bamfiles.txt -out $HOME/swift/files/swift.offtarget.angsd. -doCounts 1 -doDepth 1 -maxDepth 1000 &
nohup  $HOME/src/angsd/angsd -bam  $HOME/swift/files/swift.offtarget2.bamfiles.txt -out $HOME/swift/files/swift.offtarget2.angsd. -doCounts 1 -doDepth 1 -maxDepth 1000 &
# nohup  $HOME/src/angsd/angsd -bam  $HOME/anna/files/anna.originals.bamfiles.txt -out $HOME/anna/files/anna.originals.angsd. -doCounts 1 -doDepth 1 -maxDepth 1000 &
# nohup  $HOME/src/angsd/angsd -bam  $HOME/anna/files/anna.ontarget.bamfiles.txt -out $HOME/anna/files/anna.ontarget.angsd. -doCounts 1 -doDepth 1 -maxDepth 1000 &
# nohup  $HOME/src/angsd/angsd -bam  $HOME/anna/files/anna.offtarget.bamfiles.txt -out $HOME/anna/files/anna.offtarget.angsd. -doCounts 1 -doDepth 1 -maxDepth 1000 &
nohup  $HOME/src/angsd/angsd -bam  $HOME/anna/files/anna.offtarget2.bamfiles.txt -out $HOME/anna/files/anna.offtarget2.angsd. -doCounts 1 -doDepth 1 -maxDepth 1000 &
