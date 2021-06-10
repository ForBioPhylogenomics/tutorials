# [DATA](data)
The folder "_data_" contains the smaller data files for the afternoon session of day 4. However, we would advice you to use the data provided on SAGA as it will make the data transfer to your folder much faster. 

# [SCRIPTS](scripts)
The folder "_scripts_" contains all the scripts to be run on SAGA. However, we would advice you that you transfer them from SAGA as described in the exercise below.

# [LECTURE](Lecture)
The folder "_Lecture_" contains the lecture from this session.
* [Base composition heterogeneity](https://github.com/ForBioPhylogenomics/tutorials/blob/main/week1_day4_afternoon/Lecture/Day4_02_BaseCompositionHeterogeneity.pdf)

# EXERCISE
**On your computer**

1. The good things this afternoon is that we already have the data for the RCFV analyses. They are in the "_summarized_frequencies.txt_" file, but more data are in "_BaCoCa_Results/taxon_basefrequencies_single_partions/Matrix_original_supermatrix_mod.fas_summarized_taxon_basefrequencies.txt_" and "_BaCoCa_Results/taxon_basefrequencies_all_partions/RCFV_frequencies_all_partitions.txt_". Please download them to your computer like this morning.
	
	* Open "_Matrix_original_supermatrix_mod.fas_summarized_taxon_basefrequencies.txt_" in a text editor and delete the first line
	
2. Important these three files in R studio using "_Import Database/From text (base)_"; the heading to "_yes_"; RowNames to "_Use first column_"

3. We create now density plots, histograms and heatmaps in R to explore the distribution of the data

	* Open R script "_RCFV_analyses.R_" and execute it
	
	* explore the plots, what could be a reasonable threshold? What does the heatmap show you?
	
	* Save the results as pdf files

**On SAGA**

4. Copy all data from the folder "_Day4Afternoon_" to your folder in the project area in SAGA
	
	```
	cd /cluster/projects/nn9458k/phylogenomics/$YOURNAME
	cp -r ../week1/Day4Afternoon .
	cd Day4Afternoon
	```
	
5. We will now first prune the taxa above a certain RCFV value from the supermatrix
	
	```
	cp ../Day4Morning/BaCoCa_Results/taxon_basefrequencies_single_partions/Matrix_original_supermatrix_mod.fas_summarized_taxon_basefrequencies.txt .
	awk -F"\t" '{if($12>0.003)print$2}' < Matrix_original_supermatrix_mod.fas_summarized_taxon_basefrequencies.txt > RCFV_above0003.txt
	sh CleaningTaxaSupermatrix.sh Matrix_original_supermatrix.fas RCFV_above0003.txt
	```
	
6. When it is done, run a tree reconstruction again on the supermatrix as done before

7. The next step is do the selection at the loci level, here we have to use the value the threshold as we keep the genes
	
	```
	cp ../Day4Morning/BaCoCa_Results/summarized_frequencies.txt .
	awk -F"\t" '{if($25<0.10)print$1}' < summarized_frequencies.txt | sed "s/locus/FcC_locus/" | sed "s/$/.phy/" > RCFV_locus_below010.txt
	mkdir RCFV_locus_below010
	while read LINE; do cp ../Day2Morning/SingleGenes/$LINE RCFV_locus_below010; done < RCFV_locus_below010.txt
	```
	
8. Now we need to concatenate these again and run a tree reconstruction of the new supernatrix
	
	```
	cd RCFV_locus_below010
	```
	
	* Copy FASconCAT-G to this folder as well as the "_sbatch_Concatenation.sh_" we used yesterday (see yesterday's exercise for this) 
	
	* Modify the .sh-file to suite your needs now
	
	```
	sbatch sbatch_Concatenation.sh
	```
	
9. When it is done, run a tree reconstruction again on the supermatrix as done before

