# [DATA](data)
The folder "_data_" contains the smaller data files for the morning session of day 2. However, we would advice you to use the data provided on SAGA as it will make the data transfer to your folder much faster. 

# [SCRIPTS](scripts)
The folder "_scripts_" contains all the scripts to be run on SAGA. However, we would advice you that you transfer them from SAGA as described in the exercise below.

# [LECTURE](Lecture)
The folder "_Lecture_" contains the lecture from this session.
* [Paralogy](https://github.com/ForBioPhylogenomics/tutorials/blob/main/week1_day2_morning/Lecture/Day2_01_Paralogy.pdf)

# EXERCISE
1. For the TreSpEx analyses we first need to install some modules in your home directory as they are not included in the availables ones on SAGA.<br>

	```
	module load Perl/5.32.0-GCCcore-10.2.0
	export PERL_CPANM_HOME=/tmp/cpanm_$USER # this time the variable is a variable of the system and can be used as is.
	cpanm Statistics::LineFit
	cpanm --force Statistics::Test::WilcoxonRankSum
	export PERL5LIB=/cluster/home/$USER/perl5/lib/perl5:/cluster/home/$USER/perl5/:$PERL5LIB # this time the variable is a variable of the system and can be used as is.
	nano .bashrc
	copy the line from e) into the text editor at the end
	close the file and save it with the same name
	```
	
2. For getting started, copy all data from the folder "_Day2Morning_" to your folder in the project area in SAGA.<nr>
	
	```
	cd /cluster/projects/nn9458k/phylogenomics/$YOURNAME
	cp -r ../week1/Day2Morning .
	cd Day2Morning
	```
	
3. First we need to get a tree with bootstrap values for each individual loci.<br>
	
	```
	cd SingleGenes
	sh SingleGene_Analyses.sh
	```
	
4. For the actual paralogy screening, we need the program TreSpEx and the blast folder within it. We will also need a relaxed phylip file format for the alignment input.<br>
	
	```
	perl /cluster/projects/nn9458k/phylogenomics/week1/Programs/FASconCAT-G/FASconCAT-G_v1.05.pl -o -a -p -p -s
	cd ../
	cp -r /cluster/projects/nn9458k/phylogenomics/week1/Programs/TreSpEx/TreSpEx.v1.2_SAGA.pl /cluster/projects/nn9458k/phylogenomics/week1/Programs/TreSpEx/blast .
	```
	
5. For the paralogy screening and cleaning you will run TreSpEx on the trees with bootstrap values.<br>
	
	```
	sbatch sbatch_TreSpEx_Paralogy.sh
	```
	
6. The alignment files have been sorted into the folders "_FilesPruned_" and "_FilesNotPruned_" within the "_Results_" folder. They have to be concatenated into one file now.<br>
	
	```
	mkdir SingleGenesPara
	cp Results/FilesNotPruned/FcC_* Results/FilesPruned/FcC_* SingleGenesPara/
	cd SingleGenesPara
	cp ../../../week1/Programs/FASconCAT-G/FASconCAT-G_v1.05.pl .
	cp ../../Day1Afternoon/SingleGenes/sbatch_Concatenation.sh .
	```
	
	* Modify the "_.sh_" file to fit your needs for this analysis
	
	```
	sbatch sbatch_Concatenation.sh
	```
	
7. Next you will run a tree reconstruction to assess effect this paralogy pruning approach.<br>
	
	```
	cd Supermatrix
	cp ../../sbatch_Supermatrix_Cleaned_tree.sh .
	```
	
	* Check if the name of the supermatrix is fitting for your supermatrix. If not please change it accordingly.
	
	```
	sbatch sbatch_Supermatrix_Cleaned_tree.sh
	```

# [RESULTS](Results)
The folder "_Results_" contains the most important results from this session.
