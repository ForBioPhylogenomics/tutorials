1. For getting started, copy all data from the folder "Day1Afternoon" to your folder in the project area in SAGA<br>
	1. cd /cluster/projects/nn9458k/phylogenomics/<br>
	2. mkdir $YOURNAME<br>
	3. cd $YOURNAME<br>
	4. cp -r ../week1/Day1Afternoon .<br>
	5. cd Day1Afternoon<br>
	6. cp .ncbirc ~<br>
2. Run the following script to find the contaminations in your dataset<br>
	1. sbatch --get-user-env sbatch_ContaminationDetection.sh<br>
3. Compare "Argentina_sp_contaminated.fasta_nr_RNA_BLAST_Matches.txt" with "Argentina_sp_contaminated.fasta_Taxa_found.txt" and find out which contaminations have been found. Which datasets would be needed having in mind that we added Protodrilus symbioticus articfically as a contamination to the dataset?<br>
4. The next step would now to find proper datasets allowing for screening against positive and negative reference datasets. For this course, we already did this for you. What could be good datasets for such a screening of the entire transcriptomic or genomic datasets?<br>
	1. Please also have in mind that we will use here only the genes for screening, which are in the alignments we selected. However, you would usually do this on the entire dataset before the orthology determination.<br>
	2. If you would set up the libraries yourself, all databases for positive references would have to have a trailing "pos_" and end with ".fasta". For the negative references, it would be "neg_" instead of "pos_". The file format has to be fasta and nucleotide sequences.<br> 
5. Run the following script to generate your reference dataset:<br>
	1. sbatch sbatch_CleaningContamination.sh<br>
6. Normally you would have now a cleaned assembly dataset (the name of the assembly extended by "\_pruned.fas"), which would go into the next step determining orthologies. However, here we did it after the orthology determination and you would therefore need to clean these contaminated sequences manually from the affected alignments.<br>
		Luckily, we already did this for you. The pruned single genes are in the subfolder "SingleGenes"<br>
7. You need  to concatenate the single genes now into a supermatrix. For this you can use the program FASconCAT-G, which is written in perl.<br>
	1. cd SingleGenes<br>
	2. cp ../../week1/Programs/FASconCAT-G/FASconCAT-G_v1.05.pl .<br>
	3. sbatch sbatch_Concatenation.sh<br>
8. Now you can assess, what effect these contaminations had on the tree by running a phylogenetic tree reconstruction. You will use the same settings as before to ensure that the only difference is the exclusion of the sequences from a contamination.<br>
	1. cp Supermatrix/Matrix_Concatenated_Para_supermatrix.phy  ../ConcatenatedData/<br>
	2. cd ../ConcatenatedData/<br>
	3. sbatch sbatch_Supermatrix_Para_tree.sh<br>
9. Now you have a final tree (ending on .treefile) and compare it with the tree including the contaminated sequences. You can find this tree in the folder "Day1Morning". Are there any differences between the trees? Have in mind Argentina was the taxon with deliberate contamintion in this example and two loci were affected by this.<br>
	
