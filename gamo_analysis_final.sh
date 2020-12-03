#!/bin/bash
#
#SBATCH --job-name=biscuit_gamo_ref
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --output=output.%j.biscuit_gamo_ref
#SBATCH --time=72:00:00

module purge
module load samtools
module load python/2.7.11
module load bwa

export OMP_NUM_THREADS=$SLURM_JOB_CPUS_PER_NODE
export WORK_DIR=/scratch/general/nfs1/u0813967/gamo_fastq
export TMPDIR=/scratch/general/lustre/u0813967/tmp/SLURM_JOB_ID
export PATH=/uufs/chpc.utah.edu/common/home/u0813967/htslib/:$PATH
export PATH=/uufs/chpc.utah.edu/common/home/u0813967/bwa-meth-master/:$PATH

echo "Job started on `date`"
echo "hostname:`hostname`"

mkdir -p $TMPDIR
cd $WORK_DIR
mkdir -p methyld_methylkit
mkdir -p methyld_merged
mkdir -p mbias

gamos=(DC050 DC014 DC055 DC004 DC102 DC052 DC041 DC062 DC037 DC059 DC067 DC031 DC033 DC045 DC077 DC108 DC039 \
 DC111 DC203 DC196 DC248 DC199 DC185 DC167 DC193 DC213 DC230 DC209 DC220 DC226 DC180 DC189 DC236 DC514 DC489 \
 DC487 DC546 DC497 DC510 DC518 DC465 DC503 DC418 DC537 DC685 DC689 DC669 DC676 DC673 DC656 DC666 DC715 DC694 \
DC661 DC650 DC653 DC687 DC671 DCYY DC508 DC686 DC216 DC524 DC505 DC438 DC477 DC539 DC512 DC693 DC645 DC680)

# 5' trim will be 10 
# Create 3' trim 
gamos2=(50 50 25 50 10 10 50 25 50 50 50 40 50 50 50 60 30 30 10 10 30 10 10 40 40 \
40 10 10 10 10 10 10 40 50 20 100 30 100 50 40 50 40 100 10 50 50 50 40 10 \
40 50 40 50 50 10 10 10 10 20 10 10 40 50 20 10 10 40 50 10 10 10)


bwameth.py index gamo_ref.fa

for u in "${gamos[@]}"; do
bwameth.py --reference gamo_ref.fa ${u}F.fastq ${u}R.fastq -t 6 | samtools view -b - > ${u}bwa-meth.bam
samtools sort ${u}bwa-meth.bam -o ${u}bwameth.sorted.bam
samtools index ${u}bwameth.sorted.bam
/uufs/chpc.utah.edu/sys/installdir/methyldackel/07302018/bin/MethylDackel mbias gamo_ref.fa ${u}bwameth.sorted.bam mbias/${u}
done



for u in ${!gamos[*]}; do

/uufs/chpc.utah.edu/sys/installdir/methyldackel/07302018/bin/MethylDackel extract --methylKit gamo_ref.fa ${gamos[$u]}bwameth.sorted.bam -o methyld_methylkit/${gamos[$u]} --nOT 10,${gamos2[$u]},10,${gamos2[$u]} --nOB 10,${gamos2[$u]},10,${gamos2[$u]} --nCTOT 10,${gamos2[$u]},10,${gamos2[$u]} --nCTOB 10,${gamos2[$u]},10,${gamos2[$u]} --keepDupes
/uufs/chpc.utah.edu/sys/installdir/methyldackel/07302018/bin/MethylDackel extract --mergeContext --minDepth 10 --maxVariantFrac 0.25 gamo_ref.fa ${gamos[$u]}bwameth.sorted.bam -o methyld_merged/${gamos[$u]} --nOT 10,${gamos2[$u]},10,${gamos2[$u]} --nOB 10,${gamos2[$u]},10,${gamos2[$u]} --nCTOT 10,${gamos2[$u]},10,${gamos2[$u]} --nCTOB 10,${gamos2[$u]},10,${gamos2[$u]} --keepDupes

done

#clean up
rm -rf $TMPDIR

echo "Job ended on `date`"

# Notes: 
# Script to align and call methylation variants for mockingbird samples using bwameth and methyldackl
# Requires: mockingbird reference genome (gamo_ref.fa), as well as F and R fastqs for all samples. 
# First loop aligns sampls to the reference and then creats mbias plots to inspect alignments for methylation bias
# Second loop extracts methylation information for each position, trimming alignments based on mbias plots. 
