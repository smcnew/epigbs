#!/bin/bash
#
#SBATCH --job-name=grepgamo
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --output=output.%j.py-test
#SBATCH --time=72:00:00
#SBATCH --account=clayton
#SBATCH --partition=kingspeak
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sabrina.mcnew@gmail.com

module purge

export OMP_NUM_THREADS=$SLURM_JOB_CPUS_PER_NODE
export WORK_DIR=/scratch/general/lustre/u0813967/sabri/biscuit/gamo

echo "Job started on `date`"
echo "hostname:`hostname`"

cd $WORK_DIR

zcat R1-gamo.fastq.gz > R1-gamo.fastq
zcat R2-gamo.fastq.gz > R2-gamo.fastq

gamos=(DC050 DC014 DC055 DC004 DC102 DC052 DC041 DC062 DC037 DC059 DC067 DC031 DC033 DC045 DC077 DC108 DC039 \
 DC111 DC203 DC196 DC248 DC199 DC185 DC167 DC193 DC213 DC230 DC209 DC220 DC226 DC180 DC189 DC236 DC514 DC489 \
 DC487 DC546 DC497 DC510 DC518 DC465 DC503 DC418 DC537 DC685 DC689 DC669 DC676 DC673 DC656 DC666 DC715 DC694 \
DC661 DC650 DC653 DC687 DC671 DCYY DC508 DC686 DC216 DC524 DC505 DC438 DC477 DC539 DC512 DC693 DC645 DC680)

for u in "${gamos[@]}"; do
cat R1-gamo.fastq | grep --no-group-separator -A3 ${u} > ${u}F.fastq; \
cat R2-gamo.fastq | grep --no-group-separator -A3 ${u} > ${u}R.fastq; \
done

finch=(1P_1 10_W_1 4_P_168 9_W_167 9_W_1 5_P_1 1P_2 4_P_N 9_W_165 9_W_2 5_P_2)

for u in "${finch[@]}"; do
cat R1-zfinch.fastq | grep --no-group-separator -A3 ${u} > $u/${u}F.txt; \
cat R2-zfinch.fastq | grep --no-group-separator -A3 ${u} > $u/${u}R.txt; \
done

rm R1-gamo.fastq
rm R2-gamo.fastq


echo "Job ended on `date`"
