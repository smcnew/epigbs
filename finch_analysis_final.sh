#!/bin/bash
#
#SBATCH --job-name=methyl-dackel-zebra-finch
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --output=output.%j.methyldackl-zebra-finch
#SBATCH --time=72:00:00

module purge
module load samtools
module load bwa
module load python/2.7.11

export OMP_NUM_THREADS=$SLURM_JOB_CPUS_PER_NODE
export WORK_DIR=/scratch/general/nfs1/u0813967/finch
export TMPDIR=/scratch/general/lustre/u0813967/tmp/SLURM_JOB_ID
export PATH=/uufs/chpc.utah.edu/common/home/u0813967/bwa-meth-master/:$PATH

echo "Job started on `date`"
echo "hostname:`hostname`"

mkdir -p $TMPDIR
cd $WORK_DIR

bwameth.py index zfinch_ref.fa 

mkdir -p methyld_methylkit
mkdir -p methyld_merged
mkdir -p mbias 
finch=(1P_1 10_W_1 4_P_168 9_W_167 9_W_1 5_P_1 1P_2 4_P_N 9_W_165 9_W_2 5_P_2)

for u in "${finch[@]}"; do
bwameth.py --reference zfinch_ref.fa ${u}F.fastq ${u}R.fastq -t 6 | samtools view -b - > ${u}bwa-meth.bam
samtools sort ${u}bwa-meth.bam -o ${u}bwameth.sorted.bam
samtools index ${u}bwameth.sorted.bam
/uufs/chpc.utah.edu/sys/installdir/methyldackel/07302018/bin/MethylDackel extract --methylKit zfinch_ref.fa ${u}bwameth.sorted.bam -o methyld_methylkit/${u} --nOT 10,10,10,10 --nOB 10,10,10,10 --nCTOT 10,10,10,10 --nCTOB 10,10,10,10 --keepDupes
/uufs/chpc.utah.edu/sys/installdir/methyldackel/07302018/bin/MethylDackel extract --mergeContext --minDepth 10 --maxVariantFrac 0.25 zfinch_ref.fa ${u}bwameth.sorted.bam -o methyld_merged/${u} --nOT 10,10,10,10 --nOB 10,10,10,10 --nCTOT 10,10,10,10 --nCTOB 10,10,10,10 --keepDupes
done


# making plots of methylation bias 
#for u in "${finch[@]}"; do 
#/uufs/chpc.utah.edu/sys/installdir/methyldackel/07302018/bin/MethylDackel mbias zfinch_ref.fa ${u}bwameth.sorted.bam mbias/${u}
#done
#clean up
rm -rf $TMPDIR

echo "Job ended on `date`"

# Notes: 
# Script to align zebra finch reads to the reference and call variants. 
# Required: reference genome (zfinch_ref.fa, from refseq), as well as F and R fastqs from each bird. 
# First, index genome and then align using bwameth. Inspect reads for methylation bias, and then trim accordingly 
# using methyldackl. In this case 10 bp was trimmed from each end for each bird, based on inspection of bias plots. 
