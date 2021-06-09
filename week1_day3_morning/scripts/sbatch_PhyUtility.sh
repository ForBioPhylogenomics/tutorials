#!/bin/sh

#SBATCH --job-name=PhyUtility
#SBATCH --account=nn9408k
#SBATCH --output=PhyUtility_slurm-%j.txt
#SBATCH --time=16:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8G
#SBATCH --partition=normal 

## Set up job environment:
set -o errexit # exit on errors
set -o nounset  # Treat any unset variables as an error

module purge  # Reset the modules to the system default
module load Java/11.0.2

## Create work dir
backdir=$(pwd)
workdir=$USERWORK/$SLURM_JOB_ID
mkdir -p $workdir

## Copy input files to the work directory and move to work dir:
cp Matrix_original_supermatrix.fas.treefile Matrix_original_supermatrix.fas.ufboot $workdir
cd $workdir

# run PhyUtility to reroot to the Xenopus root
java -jar /cluster/projects/nn9458k/phylogenomics/week1/Programs/phyutility/phyutility.jar -rr -in Matrix_original_supermatrix.fas.treefile -out Matrix_original_supermatrix.fas.treefile_rooted.tre -names Xenopus_tropicalis
java -jar /cluster/projects/nn9458k/phylogenomics/week1/Programs/phyutility/phyutility.jar -rr -in Matrix_original_supermatrix.fas.ufboot -out Matrix_original_supermatrix.fas.ufboot_rooted.tre -names Xenopus_tropicalis

#calculate leaf stability indices
java -jar /cluster/projects/nn9458k/phylogenomics/week1/Programs/phyutility/phyutility.jar -ls -in Matrix_original_supermatrix.fas.ufboot_rooted.tre > PhyUtility_results.Matrix_original_supermatrix.fas.ufboot_rooted.txt

# calculate lineage movement (BAF) for Argentina
java -jar /cluster/projects/nn9458k/phylogenomics/week1/Programs/phyutility/phyutility.jar -lm -in Matrix_original_supermatrix.fas.ufboot_rooted.tre -tree Matrix_original_supermatrix.fas.treefile_rooted.tre -out Matrix_original_supermatrix.fas.treefile_BAF_Argentina.tre -names Teleost_Euteleost_Argentiniformes_Argentinidae_Argentina_sp

## Make sure the results are copied back to the submit directory:
cp -r *.txt *.tre phyutility.log $backdir
rm -rf $workdir

#### submit job as "sbatch this_script_name" ####

