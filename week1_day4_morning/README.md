# [DATA](data)
The folder "_data_" contains the smaller data files for the morning session of day 4. However, we would advice you to use the data provided on SAGA as it will make the data transfer to your folder much faster.

# [SCRIPTS](scripts)
The folder "_scripts_" contains all the scripts to be run on SAGA. However, we would advice you that you transfer them from SAGA as described in the exercise below.

# [LECTURE](Lecture)
The folder "_Lecture_" contains the lecture from this session.
* [Evolutionary rate & Saturation](https://github.com/ForBioPhylogenomics/tutorials/blob/main/week1_day4_morning/Lecture/Day4_01_EvolutionaryRateSaturation.pdf)

# EXERCISE
1. For getting started, copy all data from the folder "_Day4Morning_" to your folder in the project area in SAGA
	
	```
	cd /cluster/projects/nn9458k/phylogenomics/$YOURNAME
	cp -r ../week1/Day4Morning .
	cd Day4Morning
	```
	
2. For the calculation of the saturation indices based on slope we will need the treefiles and phylip alignment files for each orthologous loci included in the original dataset

	* Copy the _.treefile_ files and relaxed phylip files of the first 100 loci from the excerise of the morning of Day 2 to this folder
	
	```
	sbatch sbatch_TreSpEx_Saturation.sh
	```
	
3. For the c indices will only need the alignment file of the supermatrix with its partitions. Both are already in the folder.
	
	```
	sbatch sbatch_BaCoCa.sh
	```
	
4. Download the following files to your own computer using _scp_
	
	* Correlation_Results/Correlation_Slope_Summary.txt
	* BaCoCa_Results/summarized_frequencies.txt

**On your own computer**<br>

	* Open "summarized_frequencies.txt" in a text editor, delete the first line and save the file

5. Important these two txt-files in R studio using "_Import Database/From text (base)_"; the heading to "_yes_"; RowNames to "_Use first column_"

6. We create now density plots in R to explore the distribution of the data

	* Create a new R script and type in it the following:
	
	```
	x <- density(Correlation_Slope_Summary$Slope)
	plot(x)
	y <- density(Correlation_Slope_Summary$R2)
	plot(y)
	z <- density(log10(summarized_frequencies$c.value))
	plot(z)
	```
	
	* execute the R script
	
	* explore the plots, what could be a reasonable threshold?

**Back on SAGA**<br>

7. We now extract all files, which have a value above your specified threshold and which shall be included; please do the following step for all three values (c value, R2 and slope); one example is given 
	
	```
	awk -F"\t" '{if($26<100)print$1}' < BaCoCa_Results/summarized_frequencies.txt | sed "s/locus/FcC_locus/" | sed "s/$/.phy/" > summarized_frequencies_below100.txt
	mkdir Cvalue_below100
	while read LINE; do cp ../Day2Morning/SingleGenes/$LINE Cvalue_below100; done < summarized_frequencies_below100.txt
	```
	
8. Now we need to concatenate these again and run a tree reconstruction of the new supermatrix
	
	```
	cd Cvalue_above100
	```
	
	* Copy FASconCAT-G to this folder as well as the _sbatch_Concatenation.sh_ we used yesterday (see yesterday's exercise for this)
	 
	* Modify the .sh-file to suite your needs now
	
	```
	sbatch sbatch_Concatenation.sh
	```
	
9. When it is done, run a tree reconstruction again on the supermatrix as done before

