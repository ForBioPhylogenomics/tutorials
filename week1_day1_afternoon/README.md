# [DATA](data)
The folder "_data_" contains the smaller data files for the afternoon session of day 1. However, we would advice you to use the data provided on SAGA as it will make the data transfer to your folder much faster. SAGA is also the only place, where the larger data files are.

# [SCRIPTS](scripts)
The folder "_scripts_" contains all the scripts to be run on SAGA. However, we would advice you that you transfer them from SAGA as described in the exercise below.

# [LECTURES](Lectures)
The folder "_Lectures_" contains all the lectures from this session.
* [Dataset](https://github.com/ForBioPhylogenomics/tutorials/blob/main/week1_day1_afternoon/Lectures/Day1_06_Dataset.pdf)
* [Contamination](https://github.com/ForBioPhylogenomics/tutorials/blob/main/week1_day1_afternoon/Lectures/Day1_07_Contamination.pdf)

# EXERCISE
1. For getting started, copy all data from the folder "_Day1Afternoon_" to your folder in the project area in SAGA<br>

	```
	cd /cluster/projects/nn9458k/phylogenomics/
	mkdir $YOURNAME # You need to replace the variable $YOURNAME with, for example, your first or last name. This will be the folder you are working in the next two weeks.
	cd $YOURNAME # You need to replace the variable $YOURNAME with the name you chose in the line above
	cp -r ../week1/Day1Afternoon .
	cd Day1Afternoon
	cp .ncbirc ~
	```

2. Run the following script to find the contaminations in your dataset<br>

	```
	sbatch --get-user-env sbatch_ContaminationDetection.sh
	```
	
3. Go to the "_Results_" folder and compare "_Argentina_sp_contaminated.fasta_nr_RNA_BLAST_Matches.txt_" with "_Argentina_sp_contaminated.fasta_Taxa_found.txt_" and find out which contaminations have been found. Which datasets would be needed having in mind that we added _Protodrilus symbioticus_ articfically as a contamination to the dataset?

	```
	cd Results
	cat Argentina_sp_contaminated.fasta_nr_RNA_BLAST_Matches.txt
	```
	Take a look at the output before executing the next command. Which hits belong to which dataset? (__Hint:__ Look at the different name structures for the query sequences.)
	
	```
	cat Argentina_sp_contaminated.fasta_Taxa_found.txt
	cd ..
	```
	
	Check which taxa where found for the added _Protodrilus symbioticus_-dataset and which for the original _Argentina_ sp.-dataset? Which taxa will hence be relevant for the screening of contamination given that we have the artificially added one at hand already? Do we need an additional one? 
	
4. The next step would now to find proper datasets allowing for screening against positive and negative reference datasets. For this course, we already did this for you. What could be good datasets for such a screening of the entire transcriptomic or genomic datasets?<br>

	1. Please also have in mind that we will use here only the genes for screening, which are in the alignments we selected. However, you would usually do this on the entire dataset before the orthology determination.<br>
	2. If you would set up the libraries yourself, all databases for positive references would have to have a trailing "_pos\__" and end with "_.fasta_". For the negative references, it would be "_neg\__" instead of "_pos\__". The file format has to be fasta and nucleotide sequences.<br> 

5. Run the following script to generate your reference dataset:<br>
	
	```
	sbatch sbatch_CleaningContamination.sh
	```
	
6. Normally you would have now a cleaned assembly dataset (the name of the assembly extended by "_\_pruned.fas_"), which would go into the next step determining orthologies. However, here we did it after the orthology determination and you would therefore need to clean these contaminated sequences manually from the affected alignments.<br>

	Luckily, we already did this for you. The pruned single genes are in the subfolder "_SingleGenes_".
	
7. You need  to concatenate the single genes now into a supermatrix. For this you can use the program FASconCAT-G, which is written in perl.<br>

	```
	cd SingleGenes
	cp ../../week1/Programs/FASconCAT-G/FASconCAT-G_v1.05.pl .
	sbatch sbatch_Concatenation.sh
	```
	
8. Now you can assess, what effect these contaminations had on the tree by running a phylogenetic tree reconstruction. You will use the same settings as before to ensure that the only difference is the exclusion of the sequences from a contamination.<br>

	```
	cp Supermatrix/Matrix_Concatenated_Para_supermatrix.phy  ../ConcatenatedData/
	cd ../ConcatenatedData/
	sbatch sbatch_Supermatrix_Para_tree.sh
	```
	
9. Download the final tree (ending on _.treefile_) as well as the one with the contaminated sequences to your local computer. You can find the latter tree in the folder "_Day1Morning_". This will allow you to look at them using _FigTree_.
	You can use either _WinSCP_ or any other similar program to download the trees or you can use the command _scp_. In both cases, you should nagivate first to the folder, where you want the data to be on your local computer.
	The commands below are to be executed from your local computer if you are using the command line for download.

	```
	scp $USERID@saga.sigma2.no:/cluster/projects/nn9458k/phylogenomics/week1/Day1Morning/Concatenated_Para_Conta_SupermatrixTree/*.treefile . #This is for the tree with contaminated data. BE AWARE: You need to replace the variable $USERID with your user-id on Saga.
	scp $USERID@saga.sigma2.no:/cluster/projects/nn9458k/phylogenomics/$YOURNAME/week1/Day1Afternoon/Results/*.treefile . #This is for the tree with uncontaminated data. BE AWARE: You need to replace the variable $USERID with your user-id on Saga and the variable $YOURNAME with the name of the folder you generated above.
	```
	
10. Now you have a final tree (ending on _.treefile_) and compare it with the tree including the contaminated sequences. Are there any differences between the trees? Have in mind Argentina was the taxon with deliberate contamintion in this example and two loci were affected by this.<br>
	
# [RESULTS](Results)
The folder "_Results_" contains the most important results from this session.
