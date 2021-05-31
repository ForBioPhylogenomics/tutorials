#!/bin/sh

#SBATCH --job-name=Contamination
#SBATCH --account=nn9408k
#SBATCH --output=Contamination_slurm-%j.txt
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=torsths@nhm.uio.no


## Set up job environment:
set -o errexit # exit on errors
set -o nounset  # Treat any unset variables as an error

#module purge  # Reset the modules to the system default
module load BLAST+/2.11.0-gompi-2020b

## Create work dir
backdir=$(pwd)
workdir=$USERWORK/$SLURM_JOB_ID
mkdir -p $workdir

## Copy input files to the work directory and move to work dir:
cp ContaminatingTaxaDetectionPipeline.sh Argentina_sp_contaminated.fasta Barcode_Fish18S.fasta RetrieveTaxonomyDatabase.sh RetrieveTaxonPath.sh $workdir
cd $workdir

#Finding the contamination using a barcoding approach
sh ContaminatingTaxaDetectionPipeline.sh Argentina_sp_contaminated.fasta Barcode_Fish18S.fasta 1e-20 1 95

## Make sure the results are copied back to the submit directory:
cp -r Results $backdir
rm -rf $workdir

#### submit job as "sbatch this_script_name" ####
