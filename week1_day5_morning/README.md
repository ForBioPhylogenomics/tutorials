# [DATA](data)
The folder "_data_" contains the smaller data files for the afternoon session of day 1. However, we would advice you to use the data provided on SAGA as it will make the data transfer to your folder much faster. SAGA is also the only place, where the larger data files are.

# [SCRIPTS](scripts)
The folder "_scripts_" contains all the scripts to be run on SAGA. However, we would advice you that you transfer them from SAGA as described in the exercise below.

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
