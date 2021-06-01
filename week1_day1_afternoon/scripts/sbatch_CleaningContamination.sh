#!/bin/sh

#SBATCH --job-name=CleaningContamination
#SBATCH --account=nn9458k
#SBATCH --output=CleaningContamination_slurm-%j.txt
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G 


## Set up job environment:
set -o errexit # exit on errors
set -o nounset  # Treat any unset variables as an error

module purge  # Reset the modules to the system default
module load BLAST+/2.11.0-gompi-2020b
set -o errexit # exit on errors

## Create work dir
backdir=$(pwd)/Results
workdir=$USERWORK/$SLURM_JOB_ID
mkdir -p $workdir

## Copy input files to the work directory and move to work dir:
cp CleaningContamination.sh Argentina_in_alignments.fasta neg_*.fasta pos_*.fasta  $workdir
cd $workdir

#Finding the contamination using a barcoding approach
sh CleaningContamination.sh Argentina_in_alignments.fasta 10

## Make sure the results are copied back to the submit directory:
cp -r Argentina_in_alignments.fasta_* $backdir
rm -rf $workdir

#### submit job as "sbatch this_script_name" ####
