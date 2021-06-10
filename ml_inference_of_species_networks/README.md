# Maximum-Likelihood Inference of Species Networks

A tutorial on maximum-likelihood phylogenetic inference of species networks<br>
By [Michael Matschiner](https://evoinformatics.group/team.html#michaelmatschiner)

## Summary

Particularly when investigating groups of rapidly diverging species, incomplete lineage sorting (ILS) may not be the only biological process leading to variation among gene trees. Hybridization between species can lead to transfer of genetic material – so-called "introgression" between these, and thus produce gene trees that differ from the species tree. Introgression can therefore add "reticulation edges" to species trees, meaning that these are no longer bi-furcating and become a "network" instead of a tree. The identification of species networks from genomic data is still a challenging task that sees a lot of active method development. The more applicable of the currently available methods infer species networks from sets of alignments or gene trees; while a few SNP-based methods have been developed, these are so far too computationally demanding to be practical. In this tutorial, we will apply one of the earliest methods for the inference of species networks from gene trees, the "[InferNetwork\_ML](https://wiki.rice.edu/confluence/display/PHYLONET/InferNetwork_ML)" method implemented in the program [PhyloNet](https://bioinfocs.rice.edu/phylonet).


## Table of contents

* [Outline](#outline)
* [Dataset](#dataset)
* [Requirements](#requirements)
* [Extraction of alignment blocks](#alignments)
* [Scripted gene tree analyses with BEAST2](#scripted)
* [Inferring species networks with PhyloNet](#phylonet)
	* [Preparing input for PhyloNet](#preparephylonet)
	* [Running the PhyloNet analysis](#runphylonet)
	* [Visualizing species networks with IcyTree](#icytree)
	* [PhyloNet analysis without branch lengths](#nobranches)


<a name="outline"></a>
## Outline

This tutorial demonstrates how phylogenies based on whole-genome sequence data can be used to detect introgression between closely related species. The dataset used in this tutorial is a set of alignments that will be extracted from the whole-genome alignments produced in tutorial [Whole-Genome Alignment](../whole_genome_alignment/README.md). This whole-genome alignment was based on assemblies of five species of *Neolamprologus* cichlid fishes from Lake Tanganyika and the Nile tilapia as the outgroup. Phylogenies will be inferred for each alignment with BEAST2, followed by inference of the species network on the basis of these phylogeny with the "InferNetwork\_ML" method implemented in PhyloNet.


<a name="dataset"></a>
## Dataset

The dataset used in this tutorial is based on the whole-genome alignments produced in tutorial [Whole-Genome Alignment](../whole_genome_alignment/README.md). This whole-genome alignment contains sequences for five species of *Neolamprologus* cichlid fishes from Lake Tanganyika analyzed by [Gante et al. (2016)](https://doi.org/10.1111/mec.13767), as well as for Nile tilapia (*Oreochromis niloticus*), which will serve as an outgroup. 

<center>

| Species ID  | Species name                   | Tribe         |
|-------------|--------------------------------|---------------|
| orenil      | *Oreochromis niloticus*        | Oreochromini  |
| neomar      | *Neolamprologus marunguensis*  | Lamprologini  |
| neogra      | *Neolamprologus gracilis*      | Lamprologini  |
| neobri      | *Neolamprologus brichardi*     | Lamprologini  |
| neooli      | *Neolamprologus olivaceous*    | Lamprologini  |
| neopul      | *Neolamprologus pulcher*       | Lamprologini  |

</center>

Species from this group have previously been suggested to hybridize, based on mitochondrial sequences ([Salzburger et al. 2002](https://doi.org/10.1046/j.0962-1083.2001.01438.x)) and AFLP data ([Sturmbauer et al. 2010](https://doi.org/10.1016/j.ympev.2010.06.018)).

<a name="requirements"></a>
## Requirements

This tutorial requires **FigTree** to be installed. Details about the installation of this tool can be found in tutorial [Bayesian Phylogenetic Inference](../bayesian_phylogeny_inference/README.md).

The following tool is required additionally:

* **PhyloNet:** The [PhyloNet](https://bioinfocs.rice.edu/phylonet) program ([Than et al. 2008](https://doi.org/10.1186/1471-2105-9-322)) implements a range of methods for the inference of species trees and networks in a maximum-likelihood, Bayesian, or parsimony framework. PhyloNet is not available as a module on Saga, but as it is written in Java, it also does not need to be compiled and can thus simply be downloaded and placed in the current directory on Saga. This can be done with the following commands:

		wget https://bioinfocs.rice.edu/sites/g/files/bxs266/f/kcfinder/files/PhyloNet_3.8.2.jar
	
	To run PhyloNet, Java is required. To load Java, use these commands:
	
		module purge
		module load Java/11.0.2

* **babette:** [babette](https://github.com/ropensci/babette) ([Bilderbeek and Etienne (2018)](https://doi.org/10.1111/2041-210X.13032)) is an R package that allows the generation of basic XML input files for BEAST2 through scripting instead of the BEAUti GUI. Installation of this R package is required on Saga. To install it there, use the following commands:

		module load R/4.0.0-foss-2020a
		R
		install.packages("babette", repos='http://cran.us.r-project.org')
		quit(save="no")

* **ape:** [ape](http://ape-package.ird.fr) ([Paradis 2004](https://doi.org/10.1093/bioinformatics/btg412)) is an R package that serves as a multitool for basic phylogenetic analyses and handling of phylogenetic trees. To install it on Saga, use the following commands:

		R
		install.packages("ape", repos='http://cran.us.r-project.org')
		quit(save="no")



<a name="alignments"></a>
## Extraction of alignment blocks

As the whole-genome alignment produced in tutorial [Whole-Genome Alignment](../whole_genome_alignment/README.md) is stored in [HAL format](https://github.com/ComparativeGenomicsToolkit/hal/blob/master/README.md), and this format can not be read by programs for phylogenetic inference, we will need to convert the format of the dataset. At the same time, we are going to reduce the data that will go into the phylogenetic analyses, by extracting only the most suitable alignment blocks from across the whole-genome alignment rather than using the entire whole-genome alignment.

<!--XXX Polish the start, perhaps move the hal2maf conversion to the whole-genome alignment tutorial XXX-->

* As a first step, we will need to convert the whole-genome alignment from HAL format to MAF format, which is human-readable and more accessible for scripts. To do so, write a new Slurm script named `convert_hal_to_maf.slurm` on Saga with the following content:

		#!/bin/bash

		# Job name:
		#SBATCH --job-name=convert_hal_to_maf
		#
		# Wall clock limit:
		#SBATCH --time=5:00:00
		#
		# Processor and memory usage:
		#SBATCH --ntasks=1
		#SBATCH --mem-per-cpu=2G
		#
		# Accounting:
		#SBATCH --account=nn9458k
		#
		# Output:
		#SBATCH --output=convert_hal_to_maf.out

		# Set up job environment.
		set -o errexit  # Exit the script on any error
		set -o nounset  # Treat any unset variables as an error
		module --quiet purge  # Reset the modules to the system default

		# Load the hdf5 module
		module load HDF5/1.10.7-iimpi-2020b

		# Convert from hal to maf format.
		hal2maf --refGenome orenil --onlyOrthologs --noAncestors --maxBlockLen 1000000 --maxRefGap 1000 cichlids_chr5.hal cichlids_chr5.maf #XXX hal2maf is not installed yet! Running hal/bin/hal2maf so far! XXX

<!--Run time: 2:40 hours-->



* Add the script `make_alignments_from_maf.py` to your current directory on Saga, by copying it from `/cluster/projects/nn9458k/phylogenomics/week2/src` or by downloading it from GitHub, using one of the following two commands:

		cp /cluster/projects/nn9458k/phylogenomics/week2/src/make_alignments_from_maf.py .
		
	or

		wget https://raw.githubusercontent.com/ForBioPhylogenomics/tutorials/main/week2_src/make_alignments_from_maf.py

* Have a look at the options of the script `make_alignments_from_maf.py` as explained in its help text:

		module purge
		module load Python/3.8.2-GCCcore-9.3.0
		python make_alignments_from_maf.py -h
		
	You'll see that you need to specify two parameters, namely the name of an input file in MAF format and a prefix for output file name, in this order. In addition, you can specify a maximum number of alignments to be extracted from a MAF file (option `-n`), as well as a minimum length (`-l`), a maximum alignment length (`-k`), a minimum completeness (`-c`), a minimum number of sequences (`-x`), and a minimum number of polymorphic sites (`-m`) that these alignments should have. Finally, a path can be specified for where the output files should be written (`-p`).

* Run this script to extract up to 500 alignments (`-n 500`) with a length of exactly 2,500 bp (`-l 2500 -k 2500`), a completeness of at least 95% (`-c 0.95`), and a minimum of 20 polymorphic sites (`-m 20`), from the whole-genome alignment file `cichlids_chr5.maf`. Importantly, also specify that each extracted alignment should have sequences for all six species (`-x 6`). The output should be written to a directory named `alignments` (`-p alignments`), be in Fasta format `-f fasta`, and have the prefix `cichlids`:

		srun --ntasks=1 --mem-per-cpu=1G --time=00:05:00 --account=nn9458k --pty python make_alignments_from_maf.py cichlids_chr5.maf cichlids -n 500 -l 2500 -k 2500 -m 20 -c 0.95 -x 6 -p alignments -f fasta

	The screen output of this command should indicate that 500 alignment files have been written to directory `alignments`.
	
* Make sure that in fact the directory `alignments` contains 500 alignments in Fasta format:

		ls alignments/*.fasta | wc -l


<a name="scripted"></a>
## Scripted gene tree analyses with BEAST2

As PhyloNet requires "gene" (= in this context, a region of the genome, regardless of whether it's a gene or not) trees rather than alignments as input, we will still need to generate trees for each alignment. While PhyloNet can infer species networks only based on the topologies of gene trees, simulations have shown ([Yu et al. 2014](https://doi.org/10.1073/pnas.1407950111)) that the inference is improved when branch lengths are included. But when using branch lengths in the inference, a requirement of PhyloNet is that the trees also must be ultrametric, meaning that the tips of the phylogeny are all equally distant from the root. Because of this requirement, IQ-TREE is not suited for the inference of the gene trees. While recent versions of IQ-TREE (> 2.0) do come with an option to infer time-calibrated (and thus ultrametric) trees (e.g. the `--date` option) under maximum likelihood, these trees often contain zero-length branches that might not be handled well by PhyloNet. Therefore, it would be preferable to infer the gene trees with BEAST2, which also produces ultrametric trees (as you know), with more realistic branch length distributions than those produced by IQ-TREE. 

However, we obviously don't want to set up XML files with BEAUti for hundreds or thousands of alignments. Instead, it would be far more convenient to generate all XML files automatically with a script. This has been made possible with the release of the R package [babette](https://github.com/ropensci/babette) ([Bilderbeek and Etienne (2018)](https://doi.org/10.1111/2041-210X.13032)), which can produce basic XML files for BEAST only from an alignment in Fasta format and few specifications.

* Write an R script named `make_xml.r` that loads the babette library, reads an input file in Fasta format, and produces an input file for BEAST2 in XML format. The content of this script should be the following, specifying the HKY model of sequence evolution and a chain length of 500,000 iterations:

		# Load the babette library.
		library("babette")

		# Get the command-line arguments.
		args <- commandArgs(trailingOnly = TRUE)
		fasta <- args[1]
		xml <- args[2]

		# Define a site model based on the hky model.
		site_model <- create_hky_site_model()

		# Specify a run length of 500,000 iterations.
		mcmc <-  create_mcmc(chain_length = 500000)

		# Make an xml file with babette.
		create_beast2_input_file(fasta, xml, site_model = site_model, mcmc = mcmc)

* While you're at it, also write a second R script, named `convert_to_newick.r`. We will need this R script to convert the MCC trees produced by BEAST2 for each alignment into Newick format, because this will allow us to combine them into a single input file for PhyloNet more easily. This R script named `convert_to_newick.r` should have the following content:

		# Load the ape library.
		library("ape")
		
		# Get the command-line arguments.
		args <- commandArgs(trailingOnly = TRUE)
		nexus <- args[1]
		nwk <- args[2]
		
		# Read the file with a tree in nexus format.
		tree <- read.nexus(nexus)
		
		# Write a new file with a tree in newick format.
		write.tree(tree, nwk)
		
* To write an XML file for each alignment with script `make_xml.r`, analyze this XML file with BEAST2, summarize the posterior tree distribution in the form of an MCC tree, and then convert that MCC tree to Newick format with script `convert_to_newick.r`, we'll additionally need to write a Slurm script. This Slurm script should be named `prepare_phylonet.slurm` and have the following content:

		#!/bin/bash

		# Job name:
		#SBATCH --job-name=prep_phylonet
		#
		# Wall clock limit:
		#SBATCH --time=2:00:00
		#
		# Processor and memory usage:
		#SBATCH --ntasks=1
		#SBATCH --mem-per-cpu=1G
		#
		# Accounting:
		#SBATCH --account=nn9458k
		#
		# Output:
		#SBATCH --output=prepare_phylonet.out

		# Set up job environment.
		set -o errexit  # Exit the script on any error
		set -o nounset  # Treat any unset variables as an error
		module --quiet purge  # Reset the modules to the system default

		# Load modules
		module load R/4.0.0-foss-2020a
		module load Beast/2.6.4-GCC-9.3.0

		# Set a variable for the directory.
		dir=alignments
		last_char=${1}

		# Make an xml file for beast for every fasta file.
		for fasta in ${dir}/*${last_char}.fasta
		do
			xml=${fasta%.fasta}.xml
			Rscript make_xml.r ${fasta} ${xml}
		done

		# Analyze every xml file with beast.
		cd ${dir} # Go into the alignment directory to run beast there. 
		for xml in *${last_char}.xml
		do
			beast ${xml}
		done
		cd - # Return to previous directory.

		# Make an mcc summary tree for each posterior tree set.
		for trees in ${dir}/*${last_char}.trees
		do
			tre=${trees%.trees}.tre
			treeannotator -burnin 10 -heights mean ${trees} ${tre}
		done

		# Convert each mcc summary tree to newick format.
		for tre in ${dir}/*${last_char}.tre
		do
			nwk=${tre%.tre}.nwk
			Rscript convert_to_newick.r ${tre} ${nwk}
		done

	This script uses a quick-and-dirty way allowing parallelization: The line with `last_char=${1}` stores a variable that can be specified outside of this script, when calling it with `sbatch`. When we then call this script for example with `sbatch prepare_phylonet.slurm 0`, the variable "0" is passed into it and used to process only a subset of files. The line `for fasta in ${dir}/*${last_char}.fasta` for example then loops over all Fasta file in the alignment directory that end in `0.fasta`. We can therefore submit this script ten times with `sbatch prepare_phylonet.slurm 0` to `sbatch prepare_phylonet.slurm 9` to process all alignment files with ten parallel analyses. This will allow us to reduce the run time from over 2 hours to about 20 minutes.

* Submit this new Slurm script ten times with `sbatch` using the following loop:

		for i in {0..9}; do sbatch prepare_phylonet.slurm ${i}; done

	As the script runs, it should first add XML files to the alignment directory, then BEAST2's output files, then the MCC summary trees, and finally files with the MCC trees converted to Newick format. You could occasionally monitor the progress of it by checking of how many files with the endings `.xml`, `.trees`, `.tre`, and `.nwk` have already been written:

		ls alignments/*.xml | wc -l
		ls alignments/*.trees | wc -l
		ls alignments/*.tre | wc -l
		ls alignments/*.nwk | wc -l

* When the Slurm script has finished, make sure that the expected number of files with the ending `.nwk` has been written:

		ls alignments/*.nwk | wc -l


<a name="phylonet"></a>
## Inferring species networks with PhyloNet

<a name="preparephylonet"></a>
### Preparing input for PhyloNet

* Next, we need to write the input file for PhyloNet. Have a look at the examples given on the [website](https://wiki.rice.edu/confluence/display/PHYLONET/InferNetwork_ML) for PhyloNet's `InferNetwork_ML` command, which is the command that we will use, to learn about the format required for it.

	As you'll see, the file format used by PhyloNet is also Nexus format, and that it includes one block with trees as well as one block that is specific for PhyloNet and contains the commands for the PhyloNet analysis.
	
* Still on the website for PhyloNet's `InferNetwork_ML` command, scroll back to the top of the website.

	In the box just below "Usage", you should see that the `InferNetwork_ML` command can be called with just two options, called `geneTreeList` and `numReticulations`. The first of these should be a list of the IDs of the gene trees, in parentheses and separated by commas, and the second should be a number specifying the maximum number of reticulation edges that should be tested. If there are sufficiently strong introgression signals in the dataset, this number of reticulation edges will be placed on the phylogeny by PhyloNet, but if there is not sufficient support for introgression, a lower number or zero reticulation edges will be inferred. It is recommended to first limit the maximum number of reticulations, because the analysis can quickly become very computationally demanding with increasing numbers of reticulation events.
	
	All of the other options can be ignored, except the `-bl` option, with which we can specify that the branch lengths should be used in the inference.

* Write a script named `prepare_phylonet.sh` to collect the tree strings in Newick format from the files with the `.nwk` ending that were written by the `prepare_phylonet.slurm` script, and to write an input file for PhyloNet with that information. The script should have the following content:

		# Set the directory with alignments.
		dir=alignments

		# Set the name of the nexus file.
		nex= phylonet.nex

		# Write a nexus file with all trees.
		echo -e "#NEXUS\n\nBEGIN TREES;\n" > ${nex}
		count=1
		for tre in ${dir}/*.nwk
		do
			echo -n "Tree gt${count} = " >> ${nex}
			cat ${tre} >> ${nex}
			count=$(( ${count} + 1 ))
		done
		echo -e "\nEND;\n\n" >> ${nex}

		# Add a phylonet block to the nexus file.
		echo -e "BEGIN PHYLONET;\n" >> ${nex}
		echo -e "InferNetwork_ML (all) 1 -bl;" >> ${nex}
		echo -e "\nEND;" >> ${nex}

* Execute the new script `prepare_phylonet.sh` with `srun`:

		srun --ntasks=1 --mem-per-cpu=1G --time=00:01:00 --account=nn9458k --pty bash prepare_phylonet.sh
		
* When the script has finished, have a look at the file that it has written, named `phylonet.nex`, for example with `less`:

		less phylonet.nex
		
	On the line with the `InferNetwork_ML` command, you'll see that instead of the list of gene trees that should be provided, the first parameter that is passed to the command is "(all)". This is a special value that is understood by the command as a shorthand for all gene tree IDs.
	
	The second parameter that is passed to the command is the number "1", meaning that we allow the inference of up to 1 reticulation edge in this analysis.
	
	Finally, the `-bl` option is invoked, specifying that we want to use the branch lengths of the gene trees in the inference.


<a name="runphylonet"></a>
### Running the PhyloNet analysis

* To run PhyloNet with the input file `phylonet.nex`, a Slurm script is required as the last script for this analysis. Write this Slurm script with the name `run_phylonet.slurm` and the following content:

		#!/bin/bash

		# Job name:
		#SBATCH --job-name=phylonet
		#
		# Wall clock limit:
		#SBATCH --time=3:00:00
		#
		# Processor and memory usage:
		#SBATCH --ntasks=1
		#SBATCH --mem-per-cpu=10G
		#
		# Accounting:
		#SBATCH --account=nn9458k
		#
		# Output:
		#SBATCH --output=run_phylonet.out

		# Set up job environment.
		set -o errexit  # Exit the script on any error
		set -o nounset  # Treat any unset variables as an error
		module --quiet purge  # Reset the modules to the system default

		# Load the java module.
		module load Java/11.0.2

		# Run phylonet.
		java -jar PhyloNet_3.8.2.jar phylonet.nex

* Then, submit the script `run_phylonet.slurm` with `sbatch`:

		sbatch run_phylonet.slurm

* When the PhyloNet analysis has finished, have a look at the screen output that was sent to the file named `run_phylonet.out`:

		less run_phylonet.out
	
	You should see that the output lists results after 5 runs. This is the default number of searches for the maximum-likelihood species network that PhyloNet performs. In each case, the log-likelihood is reported on the same line as the network, such as in this example:
	
		Results after run #1
		-2.5192768474651035: (((neobri:1.811811478E-4,(neooli:1.322393568E-4,neopul:1.322393568E-4)I6:4.8941791E-5)I4:1.345174406E-4,((neogra:7.89246471E-5)I5#H1:1.5784929419999997E-4::0.9828573145495751,(neomar:1.578492942E-4,I5#H1:7.89246471E-5::0.017142685450424855)I2:7.892464709999999E-5)I3:7.892464710000001E-5)I1:0.007725305613599999,orenil:0.008041004202)I0;

	Unless the search for the maximum likelihood got stuck on local optima of the likelihood surface, all five log-likelihoods and networks should be identical.
	
	Following the results from the five runs, PhyloNet reports the network with the best log-likelihood (which again should be identical to at least one of the first five) on the line after "Inferred Network #1":
	
		Inferred Network #1:
		((neogra:3.724082696E-4,(neomar:3.156985884E-4,(neobri:1.811811478E-4,(neopul:1.322393568E-4,neooli:1.322393568E-4):4.8941791E-5):1.345174406E-4):5.67096812E-5):0.0076685959323999995,orenil:0.008041004202);

	If PhyloNet in fact found support for a reticulation edge (we allowed maximally one), the network should contain two occurrences of the keyword "#H1", indicating the connection points of the reticulation edge. The specification is then not strictly in Newick format anymore, but instead in "extended Newick" format, which is described in Cardona et al. ([2008](https://doi.org/10.1186/1471-2105-9-532)). The format is extremely difficult to read, and FigTree is unable to read it. To visualize the network, other programs must be used. One tool capable of reading and visualizing extended Newick format is Dendroscope ([Huson et al. 2007](https://doi.org/10.1186/1471-2105-8-460)), but we'll here use an easier-to-use alternative, the browser-based tree viewer [IcyTree](https://icytree.org) ([Vaughan 2017](https://doi.org/10.1093/bioinformatics/btx155))



<a name="icytree"></a>
### Visualizing species networks with IcyTree

XXX


* Copy the string with the network (starting with the first opening parenthesis and ending with the semi-colon) and paste it into a new file, named for example `tmp.tre`, on your local computer.

* Open the [IcyTree webpage](https://icytree.org) in a browser and load the file `tmp.tre` that you just wrote.

	**Question 1:** Did PhyloNet infer a network with a reticulation edge or merely a tree without any? [(see answer)](#q1)

* If PhyloNet inferred a reticulation edge, note the log-likelihood of the species network and repeat the inference without allowing any reticulation, and with a maximum of two allowed reticulation edges. To do that, open `phylonet.nex` in a text editor, and replace the line `InferNetwork_ML (all) 1 -bl;` with either `InferNetwork_ML (all) 0 -bl;` or `InferNetwork_ML (all) 2 -bl;`. Then, resubmit the Slurm script `run_phylonet.slurm` and wait for it to finish. The new output should then be appended to the existing file `run_phylonet.out`. Again note the log-likelihoods found for the network without reticulation edges (actually not a network but a tree) and the network with two reticulation edges (if such two edges were supported by the alignments).

* Compare the log-likelihoods for the networks with zero, one, and two reticulation edges. Assuming that the networks with more reticulation edges are identical to those with less reticulation edges except for the added reticulation edge, we can use a likelihood-ratio test to assess whether or not the additional reticulation edge significantly improved the likelihood. The table below may be helpful for this comparison.
 
	| Difference between log-likelihoods | p-value |
	|------------------------------------|---------|
	| < 1.92073                          | > 0.05  |
	| > 1.92073                          | < 0.05  |
	| > 3.31745                          | < 0.01  |
	| > 5.41375                          | < 0.001 |




<!--XXX Why is there no reticulation??? Test with SpeciesNetwork XXX-->

<a name="nobranches"></a>
### PhyloNet analysis without branch lengths

If there is time left, you could repeat the PhyloNet analysis with gene trees without branch lengths. These trees can be inferred under maximum likelihood with IQ-TREE. The comparison of the results of this PhyloNet analysis with the previous one might then shed light on the effect that the types of gene trees and the inclusion of branch lengths have.

* Write a script named `prepare_phylonet2.slurm` to run an IQ-TREE analysis for each Fasta file in the directory `alignments`. The script should have the following content:

		#!/bin/bash

		# Job name:
		#SBATCH --job-name=prep_phylonet2
		#
		# Wall clock limit:
		#SBATCH --time=1:00:00
		#
		# Processor and memory usage:
		#SBATCH --ntasks=1
		#SBATCH --mem-per-cpu=1G
		#
		# Accounting:
		#SBATCH --account=nn9458k
		#
		# Output:
		#SBATCH --output=prepare_phylonet2.out

		# Set up job environment.
		set -o errexit  # Exit the script on any error
		set -o nounset  # Treat any unset variables as an error
		module --quiet purge  # Reset the modules to the system default

		# Load modules
		module load IQ-TREE/2.1.2-foss-2020a

		# Set a variable for the directory.
		dir=alignments

		# Infer a tree for every fasta file.
		for fasta in ${dir}/*.fasta
		do
			iqtree2 -s ${fasta} -o orenil
		done

		# Clean up.
		rm alignments/*.bionj
		rm alignments/*.ckp.gz
		rm alignments/*.iqtree
		rm alignments/*.mldist
		rm alignments/*.model.gz
		rm alignments/*.fasta.log

* Submit this script with `sbatch`:

		sbatch prepare_phylonet2.slurm`
		
	This script should take around 5 minutes to complete.
	
* To prepare the new PhyloNet input file, first copy file `prepare_phylonet.sh` to a new file named `prepare_phylonet2.sh`:

		cp prepare_phylonet.sh prepare_phylonet2.sh

* Edit the new file `prepare_phylonet2.sh` so that it has the following content:

		# Set the directory with alignments.
		dir=alignments

		# Set the name of the nexus file.
		nex=phylonet2.nex

		# Write a nexus file with all trees.
		echo -e "#NEXUS\n\nBEGIN TREES;\n" > ${nex}
		count=1
		for tre in ${dir}/*.treefile
		do
			echo -n "Tree gt${count} = " >> ${nex}
			cat ${tre} >> ${nex}
			count=$(( ${count} + 1 ))
		done
		echo -e "\nEND;\n\n" >> ${nex}

		# Add a phylonet block to the nexus file.
		echo -e "BEGIN PHYLONET;\n" >> ${nex}		echo -e "InferNetwork_ML (all) 1;" >> ${nex}
		echo -e "\nEND;" >> ${nex}

	Note the changes on lines 5, 10, and the second-last line: The output Nexus file is now called `phylonet2.nex` (on line 5), the input tree files on now have the ending `.treefile` (on line 10), and the `InferNetwork_ML` command is now called without the option `-bl` (on the second-last line).
	
* Execute the new script `prepare_phylonet2.sh` with `srun`:

		srun --ntasks=1 --mem-per-cpu=1G --time=00:01:00 --account=nn9458k --pty bash prepare_phylonet2.sh

* Prepare a new Slurm script for this new PhyloNet analysis, by copying file the first Slurm script `run_phylonet.slurm` to a new file named `run_phylonet2.slurm`:

		cp run_phylonet.slurm run_phylonet2.slurm

* Change the job name to `phylonet2` (on line 4), the screen output file name to `run_phylonet2.out` (on line 17), and the input file name for PhyloNet to `phylonet2.nex` (on the last line), so that the Slurm script has the following content:

		#!/bin/bash

		# Job name:
		#SBATCH --job-name=phylonet2
		#
		# Wall clock limit:
		#SBATCH --time=3:00:00
		#
		# Processor and memory usage:
		#SBATCH --ntasks=1
		#SBATCH --mem-per-cpu=10G
		#
		# Accounting:
		#SBATCH --account=nn9458k
		#
		# Output:
		#SBATCH --output=run_phylonet2.out

		# Set up job environment.
		set -o errexit  # Exit the script on any error
		set -o nounset  # Treat any unset variables as an error
		module --quiet purge  # Reset the modules to the system default

		# Load the java module.
		module load Java/11.0.2

		# Run phylonet.
		java -jar PhyloNet_3.8.2.jar phylonet2.nex

* Then, submit this Slurm script with `sbatch`:

		sbatch run_phylonet2.slurm

	This analysis might take longer than the previous one, perhaps over an hour. You could move on to the next tutorial and come back to this one later to check the results of this second PhyloNet analysis.

* When this second PhyloNet analysis has finished, copy once again the string with the network from the screen output file `run_phylonet2.out` to a new file on your local computer, save this file under any name, and load it from the IcyTree website.

	**Question 2:** How is the resulting network affected by the different trees used as input and the fact that branch lengths were not used this time? [(see answer)](#q2)







<br><hr>

<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>

## Answers

<a name="q1"></a>

* **Question 1:** The answer may vary in each case, and depend on the alignment blocks that were extracted from the whole-genome alignment, the gene trees that were inferred by BEAST2, and whether PhyloNet found the maximum likelihood or not. All of these analysis steps included stochastic elements and may not always produce the same result.


<a name="q2"></a>

* **Question 2:** The result might be quite different from the one obtained in the first PhyloNet analysis. In my case, XXX Update with screenshot XXX, as shown in the next screenshot:

	In my view, these differences illustrate how the likelihood surface can be highly complex when species networks are inferred, and that maximum-likelihood approaches can therefore often produce results that are not particularly robust to minor changes in the model, or even to re-analyses with the same model. Species networks obtained with maximum-likelihood approaches should therefore be corroborated with other approaches, such as Bayesian inference (see tutorial [Bayesian Inference of Species Networks](../bayesian_inference_of_species_networks/README.md)) or tests of introgression based on tree-topology comparisons (see tutorial `XXX Update name XXX`) or SNPs (see tutorial [Analysis of Introgression with SNP Data](../analysis_of_introgression_with_snp_data/README.md)), before drawing conclusions about introgression.
