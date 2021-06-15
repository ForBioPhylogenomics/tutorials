#!/bin/bash
#SBATCH --job-name=mash
#SBATCH --account=nn9458k
#SBATCH --time=1:0:0
#SBATCH --mem-per-cpu=4500M
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10


module --force purge

module --ignore-cache load StdEnv Mash/2.3-GCC-10.2.0

mkdir -p mash

cd mash

#soft link in original genome assemblies (not softmasked, but that doesn't really matter)

ln -s ../../../week2/data/cichlids/original/*fasta .


#the sed at the end is to remove .fasta so we end up with just species short names
mash triangle -p 10 neobri.fasta  neogra.fasta  neomar.fasta  neooli.fasta neopul.fasta  orenil.fasta 2> tri.err | sed "s/\.fasta//g" > infile 

