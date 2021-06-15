#!/bin/bash

#SBATCH --job-name=cactus
#SBATCH --account=nn9458k
#SBATCH --time=24:0:0
#SBATCH --mem-per-cpu=10G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10

#set -o errexit  # Recommended for easier debugging

## Load your modules

module --force purge

start_dir=$PWD


mkdir -p masked_assemblies
cd masked_assemblies

#copying the fasta files containing only chr5 from Nile tilapia or sequences mapping to it
rsync /cluster/projects/nn9458k/phylogenomics/week2/data/cichlids/softmasked/*chr5.fa .
cd ..

#set up input
cat nj_tree > cactus_setup.chr5.txt

for i in $(ls masked_assemblies/*.chr5.fa); do
        j=$(basename $i)
        j=${j%.softmasked.chr5.fa}
        printf "%s /data/%s\n" $j $i >> cactus_setup.chr5.txt
done


cd $USERWORK

mkdir -p cactus_chr5
cd cactus_chr5

#singularity pull --name cactus-v1.3.0.sif docker://quay.io/comparative-genomics-toolkit/cactus:v1.3.0

mkdir -p masked_assemblies
rsync -razv ${start_dir}/masked_assemblies/*chr5.fa  masked_assemblies
cp ${start_dir}/cactus_setup.chr5.txt .

rsync  /cluster/projects/nn9458k/phylogenomics/week2/src/cactus-v1.3.0.sif . 

singularity exec -B $(pwd):/data cactus-v1.3.0.sif cactus  --maxCores 10 /data/jobStore /data/cactus_setup.chr5.txt /data/cichlids_chr5.hal --binariesMode local > ${start_dir}/cactus_chr5_1.out 2> ${start_dir}/cactus_chr5_1.err

#echo "Validation of the HAL file:" >  ${start_dir}/halValidation.txt
singularity exec -B $(pwd):/data cactus-v1.3.0.sif halValidate /data/cichlids_chr5.hal >  ${start_dir}/halValidation.chr5.txt

#echo "Statistics of the HAL file:" 
singularity exec -B $(pwd):/data cactus-v1.3.0.sif   halStats /data/cichlids_chr5.hal > ${start_dir}/halStats.chr5.txt


singularity exec -B $(pwd):/data cactus-v1.3.0.sif hal2maf --refGenome orenil --onlyOrthologs --noAncestors --maxBlockLen 1000000 --maxRefGap 1000 /data/cichlids_chr5.hal /data/cichlids_chr5.maf 

rsync cichlids_chr5.hal ${start_dir}
rsync cichlids_chr5.maf ${start_dir}
