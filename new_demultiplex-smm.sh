#!/bin/bash
#
#SBATCH --job-name=py-test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=output.%j.py-test
#SBATCH --time=72:00:00
#SBATCH --account=clayton
#SBATCH --partition=notchpeak
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sabrina.mcnew@gmail.com

### SLURM 1 processor / single node python test to run for 72 hours.
### multiple cpu cores can be requested for MPI4Py projects by adjusting --ntasks and adding mpirun to the command line
##### For example, to use 6 cores would need "--ntasks=6" and "mpirun -n 6".  Example mpirun is shown in comment below


module purge
module load python/2.7.11
export OMP_NUM_THREADS=$SLURM_JOB_CPUS_PER_NODE
export WORK_DIR=/scratch/general/lustre/u0813967/sabri/demult
export OUTPUT_DIR=/scratch/general/lustre/u0813967/files/demult2
export TMPDIR =/scratch/general/lustre/u0813967/tmp/SLURM_JOB_ID

echo "Job started on `date`"
echo "hostname:`hostname`"

# Create the OUTPUT directory if it does not exist
if [ ! -d $OUTPUT_DIR ] 
then 
     echo " Creating the directory $OUTPUT_DIR ..."
     mkdir -p $OUTPUT_DIR
fi

#WRC: 032/26/2018
#create alternative tmp dir to take care of small/tmp

mkdir -p $TMPDIR
cd $WORK_DIR


python demultiplex_MHR.py --r1_in $WORK_DIR/Mp_TKD171202028_HHY7FCCXY_L3_1_val_1.fq.gz \
                          --r2_in $WORK_DIR/Mp_TKD171202028_HHY7FCCXY_L3_2_val_2.fq.gz \
                          --barcodes $WORK_DIR/MP_barcodesFCL3.txt \
                          --output-dir $OUTPUT_DIR

#clean up 
rm -rf $TMPDIR

echo "Job ended on `date`"

#notes: 
#12 JUNE: input files changed to take the trim-galored files
#Nov 3, re-running analysis. val_1 and val_2 are outputs of trimgalore 
#Had to copy MP_barcodes...txt.complete into work directory again hope it's right..
#May 21: doing this again. 
