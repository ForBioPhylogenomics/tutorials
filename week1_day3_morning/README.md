# [DATA](data)
The folder "_data_" contains the smaller data files for the morning session of day 3. However, we would advice you to use the data provided on SAGA as it will make the data transfer to your folder much faster.

# [SCRIPTS](scripts)
The folder "_scripts_" contains all the scripts to be run on SAGA. However, we would advice you that you transfer them from SAGA as described in the exercise below.

# EXERCISE
1. Download the data from the Github data directory, or use the _.ufboot_ and _.contree_ files you obtained from the IQTree analysis yesterday afternoon.
	
2. Navigate to the Roguenarok website at: https://rnr.h-its.org/about and click "_submit job_"
	
3. Upload your _.ufboot_ file in the 'bootstrap file' space, and your _.contree_ file in the 'maximum likelihood tree' space, and click _submit_.
		
4. Once the page has loaded, select "_Strict Consensus_" as your threshold to configure your rogue taxon search on the left of the screen. Then select optimise for support, and use the Roguenarok algorithm.　Click '_Do it!_' to start the analysis

5. Repeat the analysis twice, using "_Strict Consensus_" as your threshold and optimising for support, but once using _Leaf Stability Index_ and once using _Taxonomic Instability Index_ as your algorithms. How do the results compare? Are the same sequences selected as rogues each time? Why might that be?

6. Now we will comtinue with the analyses of rogue taxa using PhyUtility.

7. For getting started, copy all data from the folder "_Day3Morning_" to your folder in the project area in SAGA

	```
	cd /cluster/projects/nn9458k/phylogenomics/$YOURNAME
	cp -r ../week1/Day3Morning .
	cd Day3Morning
	```
	
8. For calculating the leaf stability indices using PhyUtility and the BAF for _Argentina_ sp. run the following command

	```
	sbatch sbatch_PhyUtility.sh
	```
	
9. Let's take a look at the leaf stability indices. The lower the value the worse.
	
	```
	tail -40 PhyUtility_results.Matrix_original_supermatrix.fas.ufboot_rooted.txt > LeafStability_Indices_PhyUtility.txt
	sort -k2 -n LeafStability_Indices_PhyUtility.txt
	```
	
10. Download "_Matrix_original_supermatrix.fas.treefile_BAF_Argentina.tre_" to your own computer using _scp_

11. Look at the tree file with FigTree:
	What alternative positions does Argentina have? What is the frequency?
