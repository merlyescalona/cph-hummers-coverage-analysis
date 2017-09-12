#! /bin/bash
#$ -wd /home/merly/cph-bam-coverage
#$ -o /home/merly/cph-bam-coverage/output/bedtools.intersect.o
#$ -e /home/merly/cph-bam-coverage/output/bedtools.intersect.e
#$ -N bedtools.i

module load bio/bedtools/2.22.0

echo "bedtools intersect -a $(cat /home/merly/cph-bam-coverage/files/original/original.${SGE_TASK_ID}) -b $(cat /home/merly/cph-bam-coverage/files/targeted/targeted.${SGE_TASK_ID})  > $(cat /home/merly/cph-bam-coverage/files/untargeted/untargeted.${SGE_TASK_ID})"
bedtools intersect -a $(cat /home/merly/cph-bam-coverage/files/original/original.${SGE_TASK_ID}) -b $(cat /home/merly/cph-bam-coverage/files/targeted/targeted.${SGE_TASK_ID})  > $(cat /home/merly/cph-bam-coverage/files/untargeted/untargeted.${SGE_TASK_ID})

module unload bio/bedtools/2.22.0
