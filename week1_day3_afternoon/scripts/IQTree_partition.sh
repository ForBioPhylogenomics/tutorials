#!/bin/bash

# Job name:
#SBATCH --job-name=PartitionFinding
#
# Project:
#SBATCH --account=nn9458k
#
# Wall time limit:
#SBATCH --time=04:00:00
#
# Other parameters:
#SBATCH --mem-per-cpu=3G
#SBATCH	--ntasks=16
#SBATCH --output=IQTreePartition_slurm-%j.txt
## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

module --quiet purge  # Reset the modules to the system default
module load IQ-TREE/2.1.2-foss-2020a
module list

## Create work dir
backdir=$(pwd)
workdir=$USERWORK/$SLURM_JOB_ID
mkdir -p $workdir

## Copy input files to the work directory and move to work dir:
cp Matrix_Concatenated_supermatrix* $workdir
cd $workdir

## Do some work:
iqtree2-mpi -s Matrix_Concatenated_supermatrix.fas -p Matrix_Concatenated_supermatrix_partition.txt -m MFP+MERGE --prefix no_rcluster

iqtree2-mpi -s Matrix_Concatenated_supermatrix.fas -p Matrix_Concatenated_supermatrix_partition.txt -m MFP+MERGE -rcluster 10 --prefix rcluster

## Make sure the results are copied back to the submit directory:
cp * $backdir
rm -rf $workdir

#### submit job as "sbatch this_script_name" ####


