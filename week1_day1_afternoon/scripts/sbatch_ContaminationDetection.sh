#!/bin/sh

#SBATCH --job-name=Contamination
#SBATCH --account=nn9458k
#SBATCH --output=Contamination_slurm-%j.txt
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G 

## Set up job environment:
set -o errexit # exit on errors
set -o nounset  # Treat any unset variables as an error

module --quiet purge
module load BLAST+/2.13.0-gompi-2022a 

## Create work dir
backdir=$(pwd)
workdir=$USERWORK/$SLURM_JOB_ID
mkdir -p $workdir

## Copy input files to the work directory and move to work dir:
cp ContaminatingTaxaDetectionPipeline.sh Argentina_sp_contaminated.fasta Barcode_Fish18S.fasta $workdir
cd $workdir

#Finding the contamination using a barcoding approach
sh ContaminatingTaxaDetectionPipeline.sh Argentina_sp_contaminated.fasta Barcode_Fish18S.fasta 1e-20 1 95

## Make sure the results are copied back to the submit directory:
cp -r Results $backdir
rm -rf $workdir

#### submit job as "sbatch this_script_name" ####
