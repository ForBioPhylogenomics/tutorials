#!/bin/bash

#SBATCH --job-name=cactus
#SBATCH --account=nn9458k
#SBATCH --time=25:50:0
#SBATCH --partition=bigmem
#SBATCH --mem-per-cpu=32G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64

#set -o errexit  # Recommended for easier debugging

## Load your modules
#module purge   # Recommended for reproducibility

module --force purge


start_dir=$PWD

#if the softmasked genome assemblies are not present, copy them
rsync ../../week2/data/cichlids/softmasked/* masked_assemblies

#set up input
cat nj_tree > cactus_setup.txt

for i in $(ls masked_assemblies/*.softmasked.fa); do
	j=$(basename $i)
	j=${j%.softmasked.fa}
	printf "%s /data/%s\n" $j $i >> cactus_setup.txt
done

cd $USERWORK

mkdir -p cactus
cd cactus

#singularity pull --name cactus-v1.3.0.sif docker://quay.io/comparative-genomics-toolkit/cactus:v1.3.0

rsync -ravz ${start_dir}/masked_assemblies .
cp ${start_dir}/cactus_setup.txt .

rsync  /cluster/projects/nn9458k/phylogenomics/week2/src/cactus-v1.3.0.sif . 

singularity exec -B $(pwd):/data cactus-v1.3.0.sif cactus  --maxCores 64 /data/jobStore /data/cactus_setup.txt /data/cichlids_all.hal --binariesMode local > ${start_dir}/cactus_1.out 2> ${start_dir}/cactus_1.err

rsync cichlids_all.hal ${start_dir}

#echo "Validation of the HAL file:" >  ${start_dir}/halValidation.txt
singularity exec -B $(pwd):/data cactus-v1.3.0.sif halValidate /data/cichlids_all.hal >  ${start_dir}/halValidation.all.txt

#echo "Statistics of the HAL file:" 
singularity exec -B $(pwd):/data cactus-v1.3.0.sif   halStats /data/cichlids_all.hal > ${start_dir}/halStats.all.txt


#singularity exec -B $(pwd):/data cactus-v1.3.0.sif hal2maf --refGenome orenil --onlyOrthologs --noAncestors --maxBlockLen 1000000 --maxRefGap 1000 cichlids_chr5.hal cichlids_chr5.maf 
