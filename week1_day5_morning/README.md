# [DATA](data)
The folder "_data_" contains the smaller data files for the morning session of day 5. However, we would advice you to use the data provided on SAGA as it will make the data transfer to your folder much faster.

# [SCRIPTS](scripts)
The folder "_scripts_" contains all the scripts to be run on SAGA. However, we would advice you that you transfer them from SAGA as described in the exercise below.

# [LECTURE](Lecture)
The folder "_Lecture_" contains the lecture from this session.
* [Branch length heterogeneity](https://github.com/ForBioPhylogenomics/tutorials/blob/main/week1_day5_morning/Lecture/Day5_01_BranchLengthHeterogeneity.pdf)

# EXERCISE
1. Go to SAGA, copy all data from the folder "_Day5Morning_" to your folder in the project area in SAGA
	
	```
	cd /cluster/projects/nn9458k/phylogenomics/$YOURNAME
	cp -r ../week1/Day5Morning .
	cd Day5Morning
	```
	
2. We will visualize possible long-branch problems on the tree we obtained from the original 100 orthologous genes
	
	```
	sbatch sbatch_AliGroove.sh
	```
	
3. Download the _.svg_ files to computer to look at them

# [RESULTS](Results)
The folder "_Results_" contains the most important results from this session.
