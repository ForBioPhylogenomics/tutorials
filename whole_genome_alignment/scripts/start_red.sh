#list all the fasta files of the genome assemblies
for i in $(ls ../../week2/data/cichlids/original/*fasta); do
	j=${i%.fasta} #remove the suffix .fasta
	l=${j##*/} #remove everything before and including /, ending with just the species name
	#j="${basename $i .fasta}"
	echo $l
	sbatch run_red.sh $l
done
