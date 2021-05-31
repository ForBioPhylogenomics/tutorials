#!/bin/sh

#SBATCH --job-name=BaCoCa
#SBATCH --account=nn9408k
#SBATCH --output=BaCoCa_slurm-%j.txt
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=normal 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=t.h.struck@nhm.uio.no

## Set up job environment:
set -o errexit # exit on errors
set -o nounset  # Treat any unset variables as an error

module purge  # Reset the modules to the system default
module load Perl/5.32.0-GCCcore-10.2.0

## Create work dir
backdir=$(pwd)
workdir=$USERWORK/$SLURM_JOB_ID
mkdir -p $workdir

## Copy input files to the work directory and move to work dir:
cp Matrix_original_supermatrix.fas Matrix_original_supermatrix_partition.txt ../Programs/BaCoCa/BaCoCa.v1.105.pl $workdir
cd $workdir

sed '/>/s/-/_/' < Matrix_original_supermatrix.fas >> Matrix_original_supermatrix_mod.fas
sed 's/^DNA, //' < Matrix_original_supermatrix_partition.txt | sed 's/-/ - /' >> Matrix_original_supermatrix_partition_mod.txt
perl BaCoCa.v1.105.pl -i Matrix_original_supermatrix_mod.fas -p Matrix_original_supermatrix_partition_mod.txt


## Make sure the results are copied back to the submit directory:
cp -r Matrix_original_supermatrix_partition_mod.txt Matrix_original_supermatrix_mod.fas BaCoCa_Results/ $backdir
rm -rf $workdir

#### submit job as "sbatch this_script_name" ####

