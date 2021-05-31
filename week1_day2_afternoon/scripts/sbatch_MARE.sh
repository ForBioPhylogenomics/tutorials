#!/bin/sh

#SBATCH --job-name=MARE_Reduction
#SBATCH --account=nn9408k
#SBATCH --output=MARE_Reduction_slurm-%j.txt
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=normal 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=t.h.struck@nhm.uio.no

## Set up job environment:
set -o errexit # exit on errors
set -o nounset  # Treat any unset variables as an error

module purge  # Reset the modules to the system default

## Create work dir
backdir=$(pwd)
workdir=$USERWORK/$SLURM_JOB_ID
mkdir -p $workdir

## Copy input files to the work directory and move to work dir:
cp -r Matrix_original_supermatrix.fas Matrix_original_supermatrix_partition_MARE.txt MARE_v0.1.2-rc/MARE $workdir
cd $workdir

./MARE Matrix_original_supermatrix_partition_MARE.txt Matrix_original_supermatrix.fas -t 1.5 -d 1 -m -s > Matrix_original_supermatrix_t1.5_d1.0.log
mv results Matrix_original_supermatrix_t1.5_d1.0_results
./MARE Matrix_original_supermatrix_partition_MARE.txt Matrix_original_supermatrix.fas -t 1 -d 2 -m > Matrix_original_supermatrix_t1.0_d2.0.log
mv results Matrix_original_supermatrix_t1.0_d2.0_results
./MARE Matrix_original_supermatrix_partition_MARE.txt Matrix_original_supermatrix.fas -t 1 -d 3 -m > Matrix_original_supermatrix_t1.0_d3.0.log
mv results Matrix_original_supermatrix_t1.0_d3.0_results

## Make sure the results are copied back to the submit directory:
cp -r Matrix_original_supermatrix_t1.5_d1.0_results Matrix_original_supermatrix_t1.0_d2.0_results Matrix_original_supermatrix_t1.0_d3.0_results $backdir
rm -rf $workdir

#### submit job as "sbatch this_script_name" ####

