#!/bin/sh

#SBATCH --job-name=Supermatrix_Cleaned_tree
#SBATCH --account=nn9458k
#SBATCH --output=Supermatrix_Cleaned_tree_slurm-%j.txt
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16 
#SBATCH --mem-per-cpu=4G 

set -o errexit # exit on errors
set -o nounset  # Treat any unset variables as an error

module purge  # Reset the modules to the system default
module load IQ-TREE/1.6.12-foss-2018b 

## Create work dir
backdir=$(pwd)
workdir=$USERWORK/$SLURM_JOB_ID
mkdir -p $workdir

## Copy input files to the work directory and move to work dir:
cp Matrix_Concatenated_Cleaned_supermatrix.phy $workdir
cd $workdir

## Generate the phylogenetic tree of the aligned nuc sequences:
iqtree -s Matrix_Concatenated_Cleaned_supermatrix.phy -m TEST -bb 1000 -wbt -nt AUTO

## Make sure the results are copied back to the submit directory:
cp * $backdir
rm -rf $workdir

#### submit job as "sbatch this_script_name" ####


