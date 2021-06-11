# [DATA](data)
The folder "_data_" contains the smaller data files for the afternoon session of day 2. However, we would advice you to use the data provided on SAGA as it will make the data transfer to your folder much faster.

# [SCRIPTS](scripts)
The folder "_scripts_" contains all the scripts to be run on SAGA. However, we would advice you that you transfer them from SAGA as described in the exercise below.

# [LECTURE](Lecture)
The folder "_Lecture_" contains the lecture from this session.
* [Missingness & Phylogenetic signal](https://github.com/ForBioPhylogenomics/tutorials/blob/main/week1_day2_afternoon/Lecture/Day2_02_Missingness_Signal.pdf)

# EXERCISE
1. For getting started, copy all data from the folder "_Day2Afternon_" to your folder in the project area in SAGA
	
	```
	cd /cluster/projects/nn9458k/phylogenomics/$YOURNAME
	cp -r ../week1/Day2Afternoon .
	cd Day2Afternoon
	```
	
2. From the concatenation this morning we have a concatenated dataset in fasta format as well as a corresponding partition file. However, here we will proceed with the original orthologous files from the study, so that it is only a subset of the files they used and not something new. These files are already provided in the folder. However, the FASconCAT out is a bit different from the input format needed for MARE. Therefore, we first must change it.
	
	```
	sed "s/DNA,/charset/" < Matrix_original_supermatrix_partition.txt | sed "s/$/ ;/" > Matrix_original_supermatrix_partition_MARE.txt
	```
	
3. We need the program MARE.
	
	```
	wget http://software.zfmk.de/MARE_v0.1.2-rc.zip
	unzip MARE_v0.1.2-rc.zip
	cd MARE_v0.1.2-rc/
	make
	cd ../
	```
	
4. For the matrix reduction, we will run MARE with different settings for -d and -t to test their influence on the exclusion of taxa and loci
	
	```
	sbatch sbatch_MARE.sh
	```
	
5. Next we want to see what effect this has on the tree reconstruction
	
	
	* Copy the file "_sbatch_Supermatrix_MARE_tree.sh_" to the three results folders
	
	* Change them so that they fit to appropriate fasta file (check the MARE manual on the Github folder for this; what name has the reduced supermatrix?)
	
	* Submit the sbatch files

6. Next will calculate the average bootstrap support for each orthologous loci included in the original dataset
	
	* Copy the .treefile files of the first 100 loci from the excerise of this morning to the folder "_Day2Afternoon_"
	
	* "_sbatch sbatch_TreSpEx_AveBoot.sh_" from the folder "_Day2Afternoon_"
	
7. Let's take a look at the values
	
	```
	cat Average_BS_perPartition.txt
	cut -f2 Average_BS_perPartition.txt | sort
	```
	
8. Now we extract all files, which have an average bootstrap support above 70
	
	```
	awk -F" " '{if($2>70)print$1}' < Average_BS_perPartition.txt | sed "s/fas.treefile/fas/" > Average_BS_above70.txt
	mkdir Above70
	while read LINE; do cp ../Day2Morning/SingleGenes/$LINE Above70; done < Average_BS_above70.txt
	```
	
9. Now we need to concatenate these again and run a tree reconstruction of the new supernatrix
	
	```
	cd Above70
	```
	
	* Copy _FASconCAT-G_ to this folder as well as the "_sbatch_Concatenation.sh_" we used yesterday (see yesterday's exercise for this) 
	* Modify the .sh-file to suite your needs now
	
	
	```
	sbatch sbatch_Concatenation.sh
	```
	
10. When it is done, run a tree reconstruction again on the supermatrix as explained above in point 5 (the concatenated files should be in the folder "_Supermatrix_"); I would suggest to rename "_sbatch_Supermatrix_MARE_tree.sh_" to "_sbatch_Supermatrix_AveBoot_tree.sh_"

# [RESULTS](Results)
The folder "_Results_" contains the most important results from this session.
