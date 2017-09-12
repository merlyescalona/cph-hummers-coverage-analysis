################################################################################
# (c) 2015-2017 Merly Escalona <merlyescalona@uvigo.es>
# Phylogenomics Lab. University of Vigo.
# CPH Visit 2017
# ------------------------------------------------------------------------------
# Description:
# Script to extract and analyze depth of coverage information from empirical
# data that will help to model the coverage used in the simulations
################################################################################
#!/bin/bash -l
################################################################################
# Data files:
<<DATAFILES
H09_map2Anna_single-b1HiS_pair.clean0.sort.bam
H10_map2Anna_single-b1HiS_pair.clean0.sort.bam
H11_map2Anna_merge_pair.clean0.sort.bam
H12_map2Anna_merge_pair.clean0.sort.bam
H13_map2Anna_single-b1HiS_pair.clean0.sort.bam
H14_map2Anna_single-b1HiS_pair.clean0.sort.bam
H15_map2Anna_single-b1HiS_pair.clean0.sort.bam
H16_map2Anna_single-b1HiS_pair.clean0.sort.bam
H17_map2Anna_single-b1HiS_pair.clean0.sort.bam
H18_map2Anna_single-b1HiS_pair.clean0.sort.bam
H19_map2Anna_single-b1HiS_pair.clean0.sort.bam
H1_map2Anna_merge_pair.clean0.sort.bam
H20_map2Anna_single-b1HiS_pair.clean0.sort.bam
H21_map2Anna_merge_pair.clean0.sort.bam
H22_map2Anna_single-b1HiS_pair.clean0.sort.bam
H23_map2Anna_single-b1HiS_pair.clean0.sort.bam
H24_map2Anna_single-b1HiS_pair.clean0.sort.bam
H25_map2Anna_single-b1HiS_pair.clean0.sort.bam
H26_map2Anna_single-b1HiS_pair.clean0.sort.bam
H27_map2Anna_merge_pair.clean0.sort.bam
H28_map2Anna_merge_pair.clean0.sort.bam
H29_map2Anna_single-b1HiS_pair.clean0.sort.bam
H2_map2Anna_merge_pair.clean0.sort.bam
H30_map2Anna_single-b1HiS_pair.clean0.sort.bam
H31_map2Anna_merge_pair.clean0.sort.bam
H33_map2Anna_single-b1HiS_pair.clean0.sort.bam
H34_map2Anna_single-b1HiS_pair.clean0.sort.bam
H35_map2Anna_single-b1HiS_pair.clean0.sort.bam
H36_map2Anna_single-b1HiS_pair.clean0.sort.bam
H38_map2Anna_single-b1HiS_pair.clean0.sort.bam
H39_map2Anna_single-b1HiS_pair.clean0.sort.bam
H3_map2Anna_merge_pair.clean0.sort.bam
H40_map2Anna_single-b1HiS_pair.clean0.sort.bam
H41_map2Anna_single-b1HiS_pair.clean0.sort.bam
H42_map2Anna_single-b1HiS_pair.clean0.sort.bam
H43_map2Anna_single-b1HiS_pair.clean0.sort.bam
H44_map2Anna_single-b1HiS_pair.clean0.sort.bam
H45_map2Anna_single-b1HiS_pair.clean0.sort.bam
H46_map2Anna_merge_pair.clean0.sort.bam
H47_map2Anna_merge_pair.clean0.sort.bam
H48_map2Anna_single-b1HiS_pair.clean0.sort.bam
H4_map2Anna_merge_pair.clean0.sort.bam
H5_map2Anna_merge_pair.clean0.sort.bam
H6_map2Anna_merge_pair.clean0.sort.bam
H7_map2Anna_merge_pair.clean0.sort.bam
H8_map2Anna_merge_pair.clean0.sort.bam
DATAFILES

# There are 46 files - Im assuming those are 46 samples.
module load bio/samtools/1.2.0
for bamfile in $(find $HOME/cph-bam-coverage/bam/ -name "*.bam"); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  echo $tag
  samtools depth -a $bamfile > "cov/${tag}.sorted.cov"
done
module load bio/bedtols/2.22.0
# Following this process:
# http://www.gettinggeneticsdone.com/2014/03/visualize-coverage-exome-targeted-ngs-bedtools.html
# Bedtools protocol followed:
# https://github.com/arq5x/bedtools-protocols/blob/master/bedtools.md#ap2b-assessing-coverage-in-exome-capture-experiments
captureRegion="$HOME/cph-bam-coverage/Calypte_anna.gene.CDS.2750.gff"
for bamfile in $(find $HOME/cph-bam-coverage/bam/ -name "*.bam"); do
  tag=$(echo $(basename $bamfile) | tr "_" " " | awk '{print $1}')
  echo $tag
  bedtools coverage -hist -abam $bamfile -b $captureRegion> "bedtools/${tag}.cov"
done

for bamfile in $(find $HOME/cph-bam-coverage/bedtools/ -name "*.cov"); do
  tag=$(basename $bamfile .cov)
  echo $tag
  cat "${tag}.cov" | grep "^all" > "$HOME/cph-bam-coverage/bedtools/${tag}.filtered.cov"
done


for bamfile in $(find $HOME/cph-bam-coverage/bedtools/ -name "*.cov" | grep -v filtered); do
  tag=$(basename $bamfile .cov)
  echo $tag
  nTargets=$(cat $bamfile | grep -v ^all | awk '$10>0 {print $1,$4,$5}' | uniq | wc -l | awk '{print $1}')
  echo -e "$tag\t$nTargets" >> numTargets.per.sample.txt
done

for bamfile in $(find $HOME/cph-bam-coverage/bedtools/ -name "*.cov" | grep -v filtered); do
  tag=$(basename $bamfile .cov)
  echo $tag
  cat $bamfile | grep -v ^all | awk '$10==0 {print $1,$4,$5}' | uniq > "$HOME/cph-bam-coverage/targetZero/$tag.target.zero.txt"
done


# trying to see how to check this information
wc -l Calypte_anna.gene.CDS.2750.gff
24047 Calypte_anna.gene.CDS.2750.gff
