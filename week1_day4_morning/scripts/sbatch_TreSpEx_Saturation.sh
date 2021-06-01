#!/bin/sh

#SBATCH --job-name=TreSpEx_Saturation
#SBATCH --account=nn9458k
#SBATCH --output=TreSpEx_Saturation_slurm-%j.txt
#SBATCH --time=1:00:00
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
cp -r *.phy *.treefile ../Programs/TreSpEx/TreSpEx.v1.2_SAGA.pl $workdir

cd $workdir

ls *.treefile > Treefiles_Saturation.txt
sed "s/fas.treefile/phy/" < Treefiles_Saturation.txt | sed "s/locus/FcC_locus/" > Alignments_Saturation.txt

perl TreSpEx.v1.2_SAGA.pl -fun g -ipt Treefiles_Saturation.txt -ipa Alignments_Saturation.txt

## Make sure the results are copied back to the submit directory:
cp -r Treefiles_Saturation.txt Alignments_Saturation.txt Correlation_Results/ $backdir
rm -rf $workdir

#### submit job as "sbatch this_script_name" ####
