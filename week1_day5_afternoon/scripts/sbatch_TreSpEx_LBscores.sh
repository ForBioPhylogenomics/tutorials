#!/bin/sh

#SBATCH --job-name=TreSpEx_LBscores
#SBATCH --account=nn9458k
#SBATCH --output=TreSpEx_LBscores_slurm-%j.txt
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=normal 

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
cp -r Matrix_original_supermatrix.fas *.treefile ../Programs/TreSpEx/TreSpEx.v1.2_SAGA.pl $workdir
cd $workdir

grep "^>" < Matrix_original_supermatrix.fas | sed "s/>//" > Taxa_Matrix_original_supermatrix.txt
ls *.treefile > Treefiles_LBscores.txt

perl TreSpEx.v1.2_SAGA.pl -fun e -ipt Treefiles_LBscores.txt -tf Taxa_Matrix_original_supermatrix.txt 

## Make sure the results are copied back to the submit directory:
cp -r *.txt *.log LBscore_Results $backdir
rm -rf $workdir

#### submit job as "sbatch this_script_name" ####
