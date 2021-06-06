# [DATA](data)
The folder "_data_" contains the smaller data files for the afternoon session of day 3. However, we would advice you to use the data provided on SAGA as it will make the data transfer to your folder much faster.

# [SCRIPTS](scripts)
The folder "_scripts_" contains all the scripts to be run on SAGA. However, we would advice you that you transfer them from SAGA as described in the exercise below.

# EXERCISE
1. For getting started, copy all data from the folder "_Day3Afternoon_" to your folder in the project area in SAGA
	
	```
	cd /cluster/projects/nn9458k/phylogenomics/$YOURNAME
	cp -r ../week1/Day3Afternoon .
	cd Day3Afternoon
	```
	
2. Run the IQTree script provided in your directory
	
	```
	sbatch IQTree_partition.sh
	```
	
3. Download the two trees and two log files to your local computer.

4. With an additional partitioning analysis, how do these two trees compare to the trees from yesterday afternoon? How do they compare to one another?
