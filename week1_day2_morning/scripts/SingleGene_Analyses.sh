#!/bin/bash

mkdir Results
FOLDER=$(pwd)/Results
for FILE in locus*.fas
do
echo "#!/bin/sh" >> iqtree_${FILE}.sh
echo "" >> iqtree_${FILE}.sh
echo "### Here change the job name and replace the species in accordance with species analysed ###" >> iqtree_${FILE}.sh
echo "#SBATCH --job-name=FORBIO_SingleGene_${FILE}" >> iqtree_${FILE}.sh
echo "" >> iqtree_${FILE}.sh
echo "# Leave this unchanged" >> iqtree_${FILE}.sh
echo "#SBATCH --account=nn9408k" >> iqtree_${FILE}.sh
echo "" >> iqtree_${FILE}.sh
echo "### Here change the job name and replace the species in accordance with species analysed ###" >> iqtree_${FILE}.sh
echo "#SBATCH --output=FORBIO_SingleGene_${FILE}_slurm-%j.txt" >> iqtree_${FILE}.sh
echo "" >> iqtree_${FILE}.sh
echo "### Adjust memory if needed ###" >> iqtree_${FILE}.sh
echo "#SBATCH --mem-per-cpu=8G" >> iqtree_${FILE}.sh
echo "" >> iqtree_${FILE}.sh
echo "### Adjust the number of tasks if needed ###" >> iqtree_${FILE}.sh
echo "#SBATCH --ntasks=1" >> iqtree_${FILE}.sh
echo "" >> iqtree_${FILE}.sh
echo "### Adjust run wall time if needed ###" >> iqtree_${FILE}.sh
echo "#SBATCH --time=2-00:00:00" >> iqtree_${FILE}.sh
echo "" >> iqtree_${FILE}.sh
echo "# Leave this unchanged" >> iqtree_${FILE}.sh
echo "#SBATCH --mail-type=FAIL" >> iqtree_${FILE}.sh
echo "#SBATCH --mail-user=t.h.struck@nhm.uio.no" >> iqtree_${FILE}.sh
echo "" >> iqtree_${FILE}.sh

echo "set -o errexit # exit on errors" >> iqtree_${FILE}.sh
echo "set -o nounset  # Treat any unset variables as an error" >> iqtree_${FILE}.sh
echo "" >> iqtree_${FILE}.sh

echo "module --quiet purge  # Reset the modules to the system default" >> iqtree_${FILE}.sh
echo "export LMOD_DISABLE_SAME_NAME_AUTOSWAP=no" >> iqtree_${FILE}.sh
echo "module load IQ-TREE/1.6.12-foss-2018b" >> iqtree_${FILE}.sh
echo "" >> iqtree_${FILE}.sh

echo "## Create work dir" >> iqtree_${FILE}.sh
echo "workdir=\$USERWORK/\$SLURM_JOB_ID" >> iqtree_${FILE}.sh
echo "mkdir -p \$workdir" >> iqtree_${FILE}.sh
echo "" >> iqtree_${FILE}.sh

echo "## Copy input files to the work directory and move to work dir:" >> iqtree_${FILE}.sh
echo "cp ${FILE} \$workdir" >> iqtree_${FILE}.sh
echo "cd \$workdir" >> iqtree_${FILE}.sh
echo "" >> iqtree_${FILE}.sh

echo "## Generate the phylogenetic tree of the aligned nuc sequences:" >> iqtree_${FILE}.sh
echo "iqtree -s ${FILE} -m TEST -bb 1000 -wbt -nt AUTO" >> iqtree_${FILE}.sh
echo "" >> iqtree_${FILE}.sh

echo "## Make sure the results are copied back to the submit directory:" >> iqtree_${FILE}.sh
echo "cp * ${FOLDER}" >> iqtree_${FILE}.sh
echo "rm -rf ${workdir}" >> iqtree_${FILE}.sh

sbatch iqtree_${FILE}.sh
done