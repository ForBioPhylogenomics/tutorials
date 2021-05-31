#!/bin/sh

#SBATCH --job-name=TreSpEx_AveBoot
#SBATCH --account=nn9408k
#SBATCH --output=TreSpEx_AveBoot_slurm-%j.txt
#SBATCH --time=1:00:00
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
cp -r *.treefile ../Day2Morning/TreSpEx.v1.2_SAGA.pl $workdir
cd $workdir

# run Paralogy screening - bootstrap screening
ls *.treefile > Treefiles_Paralogy.txt

perl TreSpEx.v1.2_SAGA.pl -fun f -ipt Treefiles_Paralogy.txt
  
## Make sure the results are copied back to the submit directory:
cp -r Treefiles_Paralogy.txt Average_BS_Results/Average_BS_perPartition.txt $backdir
rm -rf $workdir

#### submit job as "sbatch this_script_name" ####

