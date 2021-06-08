#!/bin/sh

#SBATCH --job-name=TreSpEx_Paralogy
#SBATCH --account=nn9458k
#SBATCH --output=TreSpEx_Paralogy_slurm-%j.txt
#SBATCH --time=16:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=normal 

## Set up job environment:
set -o errexit # exit on errors
set -o nounset  # Treat any unset variables as an error

module purge  # Reset the modules to the system default
module load Perl/5.32.0-GCCcore-10.2.0
module load BLAST+/2.11.0-gompi-2020b

## Create work dir
mkdir Results
backdir=$(pwd)/Results
workdir=$USERWORK/$SLURM_JOB_ID
mkdir -p $workdir

## Copy input files to the work directory and move to work dir:
cp -r SingleGenes/Results/*.treefile SingleGenes/*.phy TreSpEx.v1.2_SAGA.pl blast $workdir
cd $workdir

# run Paralogy screening - bootstrap screening
ls *.treefile > Treefiles_Paralogy.txt
sed "s/fas.treefile/phy/" < Treefiles_Paralogy.txt | sed "s/locus/FcC_locus/" > Alignments_Paralogy.txt

perl TreSpEx.v1.2_SAGA.pl -fun a -ipt Treefiles_Paralogy.txt -gts N -lowbs 95 -upbs 100 -possc 0 -poslb 0 -possb 0 

# run Paralogy screening - blast search

perl TreSpEx.v1.2_SAGA.pl -fun c -ppf PotentialParalogsBootstrap.txt -ipt Treefiles_Paralogy.txt -ipa Alignments_Paralogy.txt -db1 Bos_taurus -db2 Branchiostoma_floridae -ediff 5 -ltp 0.05 -utp 0.95 -evalue 1e-20 

# run Paralogy screening - pruning of sequences

cp blast_search/blast_comparisons/Certain_paralogy/Certain_PotentialParalogsBootstrap.txt .

perl TreSpEx.v1.2_SAGA.pl -fun d -ppf Certain_PotentialParalogsBootstrap.txt -ipt Treefiles_Paralogy.txt -ipa Alignments_Paralogy.txt 
  
## Make sure the results are copied back to the submit directory:
cp -r *.csv *.log *.txt blast_search FilesNotPruned FilesPruned LBfactor* $backdir
rm -rf $workdir

#### submit job as "sbatch this_script_name" ####

