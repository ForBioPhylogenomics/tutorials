#!/bin/sh

#SBATCH --job-name=AliGROOVE
#SBATCH --account=nn9408k
#SBATCH --output=AliGROOVE_slurm-%j.txt
#SBATCH --time=16:00:00
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
cp -r *.fas *.treefile ../Programs/aligroove_1_07/* $workdir
cd $workdir

sed '/>/s/-/_/' < Matrix_original_supermatrix.fas >> Matrix_original_supermatrix.fas_mod
sed 's/-/_/g' < Matrix_original_supermatrix.fas.treefile >> Matrix_original_supermatrix.fas.treefile.tre
perl AliGROOVE_v.1.07.pl -i Matrix_original_supermatrix.fas_mod -z Matrix_original_supermatrix.fas.treefile.tre -N
#done

## Make sure the results are copied back to the submit directory:
cp -r *.svg *.txt $backdir
rm -rf $workdir

#### submit job as "sbatch this_script_name" ####
