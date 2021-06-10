# [SCRIPTS](scripts)
The folder "_scripts_" contains all the scripts to be run on SAGA. However, we would advice you that you transfer them from SAGA as described in the exercise below.

# EXERCISE
1. For getting started, copy all data from the folder "_Day5Afternoon_" to your folder in the project area in SAGA
	
	```
	cd /cluster/projects/nn9458k/phylogenomics/$YOURNAME
	cp -r ../week1/Day5Afternoon .
	cd Day5Afternoon
	```
	
2. For the calculation of the LB scores we will need the treefiles of the original supermatrix and the individual loci. Additionally, we need the fasta alignment file for the original dataset of 100 loci:

	* Copy the _.treefile_ files of individual loci you used in the excerise of the morning of Day 4 to this folder
	
	* Do the same for the _.treefile_ file and the _.fas_ file from the exercise of this morning for the original supermatrix
	
	```
	sbatch sbatch_TreSpEx_LBscores.sh
	```

3. Download the _.txt_ files in folder "_LBscore_Results_" to your own computer using _scp_

**On your own computer**

4. Important these three files in R studio using "_Import Database/From text (base)_"; the heading to "_yes_"; RowNames to "_Use first column_"

5. We create now density plots, histograms, correlation plots and heatmaps in R to explore the distribution of the data

	* Open R script "_LBscore_analyses.R_" and execute it
	
	* explore the plots, what could be reasonable thresholds? What do the correlation plots tell you? What does the heatmap show you?
	
	* Save the results as pdf files

**Back on SAGA**

6. We will now first prune the taxa above a certain LB or TR score value from the supermatrix (so two times)
	
	```
	cp LBscore_Results/*.txt .
	awk -F"\t" '{if($2>20)print$1}' < LB_scores_perTaxon.txt > LB_above20.txt
	sh CleaningTaxaSupermatrix.sh Matrix_original_supermatrix.fas LB_above20.txt
	mv Matrix_original_supermatrix.fas_pruned.fas Matrix_original_supermatrix.fas_LBpruned.fas
	```
	
	* Repeat this for _TR_scores_perTaxon.txt_ (remember to adjust your code)

7. When it is done, run a tree reconstruction again on each pruned supermatrix as done before

8. We now extract all files, which have a value below your specified threshold and which shall be included; please do the following step for both values (LB score heterogeneity, Average PD); one example is given
	
	```
	awk -F"\t" '{if($3<50)print$1}' < LB_scores_summary_perPartition.txt | sed "s/locus/FcC_locus/" | sed "s/.fas.treefile/.phy/" > LBscores_below50.txt
	mkdir LBscores_below50
	while read LINE; do cp ../Day2Morning/SingleGenes/$LINE LBscores_below50; done < LBscores_below50.txt
	```

9. Now we need to concatenate these again and run a tree reconstruction of the new supernatrix
	
	```
	cd LBscores_below50
	```
	
	* Copy FASconCAT-G to this folder as well as the _sbatch_Concatenation.sh_ we used yesterday (see yesterday's exercise for this) 
	
	* Modify the _.sh_-file to suite your needs now
	
	```
	sbatch sbatch_Concatenation.sh
	```

10. When it is done, run a tree reconstruction again on the supermatrix as done before

