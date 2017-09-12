#! /bin/bash
#$ -wd /home/merly/cph-bam-coverage/angsd-depth
#$ -o /home/merly/cph-bam-coverage/output/angsd.untargeted.coverage.o
#$ -e /home/merly/cph-bam-coverage/output/angsd.untargeted.coverage.e
#$ -N angsd-untarget

module load bio/angsd/0.902 bio/bedtools/2.22.0

bedtools intersect -a /home/merly/cph-bam-coverage/originals/H3_map2Anna_merge_pair.clean0.sort.bam -b /home/merly/cph-bam-coverage/bam-targeted-captured/H3.filtered.bam > /home/merly/cph-bam-coverage/bam-untargeted-captured/H3.untargeted.captured.bam

angsd -bam /home/merly/cph-bam-coverage/untargeted.captured.txt -doDepth 1 -out untargeted-captured -doCounts 1

module unload bio/angsd/0.902  bio/bedtools/2.22.0
bam-untargeted-captured/H4.untargeted.captured.bam'
bam-untargeted-captured/H4.untargeted.captured.bam
