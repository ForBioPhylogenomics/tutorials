#!/bin/bash

#SBATCH --job-name=red
#SBATCH --account=nn9458k
#SBATCH --time=4:0:0
#SBATCH --mem-per-cpu=20000M
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

#set -o errexit  # Recommended for easier debugging

module purge

#might not need --ignore-cache anymore
module --ignore-cache load Red/2015-05-22-GCC-7.3.0-2.30 Biopython/1.72-foss-2018b-Python-3.6.6

mkdir -p masked_assemblies

start_dir=$PWD

#running from $USERWORK since this creates quite a but of files, overloading the number of files quota on the course folder
cd $USERWORK


/cluster/projects/nn9458k/phylogenomics/week2/src/redmask.py -i ${start_dir}/../../week2/data/cichlids/original/${1}.fasta -o ${start_dir}/masked_assemblies/${1}






