#!/bin/sh

#SBATCH --job-name=Concatenation
#SBATCH --account=nn9458k
#SBATCH --output=Concatenation_slurm-%j.txt
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G 


## Set up job environment:
set -o errexit # exit on errors
set -o nounset  # Treat any unset variables as an error

module purge  # Reset the modules to the system default
module load Perl/5.32.0-GCCcore-10.2.0

## Create work dir
mkdir Supermatrix
backdir=$(pwd)/Supermatrix
workdir=$USERWORK/$SLURM_JOB_ID
mkdir -p $workdir

## Copy input files to the work directory and move to work dir:
cp FASconCAT-G_v1.05.pl locus_*.fas $workdir
cd $workdir

#Finding the contamination using a barcoding approach
perl FASconCAT-G_v1.05.pl -s -p -p -n -l
mv FcC_info.xls Matrix_Concatenated_Para_info.xls
mv FcC_supermatrix.fas Matrix_Concatenated_Para_supermatrix.fas
mv FcC_supermatrix.nex Matrix_Concatenated_Para_supermatrix.nex
mv FcC_supermatrix.phy Matrix_Concatenated_Para_supermatrix.phy
mv FcC_supermatrix_partition.nex Matrix_Concatenated_Para_supermatrix_partition.nex
mv FcC_supermatrix_partition.txt Matrix_Concatenated_Para_supermatrix_partition.txt

## Make sure the results are copied back to the submit directory:
cp -r Matrix_Concatenated_Para* $backdir
rm -rf $workdir

#### submit job as "sbatch this_script_name" ####
