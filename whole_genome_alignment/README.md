# Whole Genome Alignment

A tutorial on multiple whole genome alignment Cactus<br/>
By Ole K. Tørresen

## Summary

Large scale projects such as Earth BioGenome Project ([Lewin et al. (2018)](https://doi.org/10.1073/pnas.1720115115), the Vertebrate Genomes Project [Rhie et al. (2021)](https://doi.org/10.1038/s41586-021-03451-0)), and other similar efforts, aspire to sequence and assemble the genomes of all living species. More and more high quality genome assemblies are therefore being deposited in public databases, and these resources contain great potential and enables rigorous investigation into a multitude of aspects of biology through data - and curiosity driven research. The further generation of high-quality genomes will facilitate great leaps in our understanding of, and ability to answer, core questions within fundamental biology such as in ecology and evolution, and in other fields such as medicine.

However, to actually understand these genomes, we need some way of analysing them. For smaller projects, sequencing and assembling a few species, annotating and analysing each can be done focussing on each in turn. However, for large projects such as Zoonimia (131 new genome assemblies) ([Zoonomia Consortium (2020)](https://doi.org/10.1038/s41586-020-2876-6)) or bird 10k (B10K; 267 new genome assemblies) ([Feng et al. (2020)](https://doi.org/10.1038/s41586-020-2873-9)), using that much time and effort on each genome is not scalable, and therefore they used whole-genome alignments to analyse their datasets. So far, the only whole-genome aligner that has been used successfully on such large datasets is Cactus ([Armstrong et al. 2020](https://doi.org/10.1038/s41586-020-2871-y)).

## Table of contents

* [Outline](#outline)
* [Dataset](#dataset)
* [Requirements](#requirements)
* [Softmasking with Red](#softmask)
* [Creating a "phylogeny" with Mash](#mash)
* [Running Cactus](#cactus)

<a name="outline"></a>
## Outline

In this tutorial we will run six genome assemblies of varying quality through this pipeline. This requires masking repeats in the genomes and a phylogeny (or "phylogeny") of the species involved. Unfortunately, most of these processes take longer than the time we have available in this session, but all intermediate and final files are available on Saga.


<a name="dataset"></a>
## Dataset

The dataset used in this tutorial consists of a subset of species you have used earlier in this course. They are all cichlids, but five of them are quite closely related with Nile tilapia as outgroup.

<center>

| ID      | Species                       | Common name               | Genus            |
|---------|-------------------------------|---------------------------|------------------|
| neobri  | *Neolamprologus brichardi*    | lyretail cichlid          | *Neolamprologus* |
| neomar  | *Neolamprologus marunguensis* | ?                         | *Neolamprologus* |
| neogra  | *Neolamprologus gracilis*     | ?                         | *Neolamprologus* |
| neooli  | *Neolamprologus olivaceous*   | ?                         | *Neolamprologus* |
| neopul  | *Neolamprologus pulcher*      | daffodil cichlid          | *Neolamprologus* |
| orenil  | *Oreochromis niloticus*       | Nile tilapia              | *Oreochromis*    |


</center>

<a name="requirements"></a>
## Requirements

We will use [Red](http://toolsmith.ens.utulsa.edu/) (but via the Python script [redmask](https://github.com/nextgenusfs/redmask)), [Mash](https://github.com/marbl/Mash), [rapidNJ](https://github.com/somme89/rapidNJ) and [Cactus](https://github.com/ComparativeGenomicsToolkit/cactus). Red, Mash and rapidNJ should be available as modules on Saga, while Cactus works best as an container (the installation is non-trivial.)



<a name="softmask"></a>
## Softmasking with Red

A detailed exposition of alignment is outside the scope of this tutorial (and I guess most of you know it quite well). However, to understand this step, we need to review some steps of it. Broadly, to create an alignment two sequences needs to be compared in some manner. The classic algorithms here are [Smith-Waterman](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm) and [Needleman-Wunsch](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm), for respectively local and global alignments. In local alignment the most similar sequence is found between the two input sequences, while a global alignment requires that the whole length of both input sequences are involved. For larger sequences (chromosomes for instance) or for multiple comparisons (between proteins and a genome), the requirements in computation (both memory and CPU) becomes unpractical. To remedy this, different approaches have been developed, such as [BLAST](https://en.wikipedia.org/wiki/BLAST_(biotechnology) where words are connected to the sequences and only sequences with words in common are compared.

Repetitive sequences are sequences that are identical or very similar that occur multiple places in a genome. For instance the short tandem repeat ACT repeated both in tandem (ACTACTACTACT), but also multiple places, or different transposable elements that have copied themselves across the genome. These repetitive sequences will lead to large numbers of comparisons unless they are taken out of the picture in some way. The usual approach is to mask them, either hard mask (replace the repeats with Ns) or soft mask (lower case letters; standard fasta sequence is upper case). If a sequence is hard masked, it is basically gone from the analysis, so the recommendation is to soft mask. Alignments will not be initiated in soft masked regions, but alignments might extend into soft mask regions and therefore lead to less fragmented alignments.

Maybe most common today is to use [RepeatModeler](http://www.repeatmasker.org/RepeatModeler/) [Flynn et al (2020)](https://doi.org/10.1073/pnas.1921046117) to create a species specific repeat library, and use that in [RepeatMasker](https://www.repeatmasker.org/) to mask the genome. That process is a bit circumstantial so we'll take a short cut here (but feel free to use it in your own analysis).

Red [Giris (2015)](https://doi.org/10.1186/s12859-015-0654-5) can do the softmasking in one step instead of two, and is faster. We'll run it via the Python script [redmask](https://github.com/nextgenusfs/redmask) to get readily softmasked genomes.

First, go to your user folder on Saga and copy the needed scripts:

```
cd /cluster/projects/nn9458k/phylogenomics/
mkdir $YOURNAME
cd $YOURNAME
mkdir -p cichlids
cd cichlids
cp ../..//week2/data/cichlids/scripts/* .

```

The redmask command itself is quite straight forward:
```
redmask.py -i genome.fa -o mygenome
```

But we'll wrap in an SBATCH script and run it for all genomes. Unfortunately, redmask requires a Python module not readily available on Saga. To set install it locally (for you as a user), do this:
```
moduleload Biopython/1.72-foss-2018b-Python-3.6.6
pip3 install natsort --user
```
We load the same module that we use in the script to ensure that we install natsort for the correct version of Python.

Then you can do this:

```
sh start_red.sh
```
This little snippet starts one script per genome. Feel free to look around in it.

Be aware that this creates lots of small files which exceeded the quota for the course. run_red.sh now navigates to $USERWORK and runs the job there instead.

When I ran through this, it took 15 to 60 minutes per genome to finish.

If you are unable to get it running for some reason, or want to continue anyhow, you can do this:

```
mkdir -p masked_assemblies
cd masked_assemblies
rsync ../../../week2/data/cichlids/softmasked/* .
cd ..
```
When done, you should see something like this in masked_assemblies
```
-rw-rw-r-- 1 olekto nn9458k   43497196 Jun  9 21:40 neobri.repeats.bed
-rw-rw-r-- 1 olekto nn9458k  275320088 Jun  9 21:40 neobri.repeats.fasta
-rw-rw-r-- 1 olekto nn9458k  865018698 Jun  9 21:40 neobri.softmasked.fa
-rw-rw-r-- 1 olekto nn9458k   45092622 Jun  9 22:16 neogra.repeats.bed
-rw-rw-r-- 1 olekto nn9458k  262196845 Jun  9 22:16 neogra.repeats.fasta
-rw-rw-r-- 1 olekto nn9458k  672846647 Jun  9 22:14 neogra.softmasked.fa
-rw-rw-r-- 1 olekto nn9458k   46494946 Jun  9 22:00 neomar.repeats.bed
-rw-rw-r-- 1 olekto nn9458k  269948250 Jun  9 22:01 neomar.repeats.fasta
-rw-rw-r-- 1 olekto nn9458k  683501275 Jun  9 21:59 neomar.softmasked.fa
-rw-rw-r-- 1 olekto nn9458k   46074867 Jun  9 22:02 neooli.repeats.bed
-rw-rw-r-- 1 olekto nn9458k  267874589 Jun  9 22:03 neooli.repeats.fasta
-rw-rw-r-- 1 olekto nn9458k  681893391 Jun  9 22:01 neooli.softmasked.fa
-rw-rw-r-- 1 olekto nn9458k   45076389 Jun  9 22:03 neopul.repeats.bed
-rw-rw-r-- 1 olekto nn9458k  262152533 Jun  9 22:03 neopul.repeats.fasta
-rw-rw-r-- 1 olekto nn9458k  674817636 Jun  9 22:02 neopul.softmasked.fa
-rw-rw-r-- 1 olekto nn9458k   30337323 Jun  9 21:52 orenil.repeats.bed
-rw-rw-r-- 1 olekto nn9458k  347698289 Jun  9 21:53 orenil.repeats.fasta
-rw-rw-r-- 1 olekto nn9458k 1025835685 Jun  9 21:52 orenil.softmasked.fa
```

<a name="mash"></a>
## Creating a "phylogeny" with Mash

In addition to softmasked genome assemblies, Cactus requires a guide tree. This is because it aligns two and two genomes. For instance, it could start with the two most closely related, aligns those, and then use that alignment to reconstruct an ancestral genome to those two. That reconstructed genome is then used to align against the next. This means that Cactus scales linear with the number of genomes and not quadratically as would have been the case if all genomes were compared against all. However, this approach might miss alignment of sequences where species A and C has it, and B not, but since A and B were closest related, those were aligned first and the reconstructed genome lacked the shared regions between A and C and therefore the resulting alignment didn't show that case.

Creating phylogenies can be a time consuming process, and therefore we will cheat (the "phylogeny" in the title is there for a reason). We will use Mash ([Ondov et al. (2016)](https://doi.org/10.1186/s13059-016-0997-x)), which is a fast genome and metagenome distance estimator. Mash reduces large sequences to compressed sketch representations, enabling much quicker comparisons between sequences compared to older methods.

Since we use the original genome assembly fasta files, and not the softmasked ones, you can run this independent of the softmasking.

If you followed the instructions above, you should have a script called run_mash_triangle.sh in your folder. Submit that:

```
sbatch run_mash_triangle.sh
```

This runs quite quickly, less than 1 minute when I did it. Considering that it is comparing the content of six fasta files containing genomes between 600 and 1000 Mbp, this is quite nice.

The result, in the subfolder mash is this:
```
$cat mash/infile
	6
neobri
neogra	0.00981378
neomar	0.00985499	0.00888862
neooli	0.00904645	0.00948708	0.0090861
neopul	0.00924552	0.00924552	0.00932572	0.00659353
orenil	0.0519203	0.0509537	0.0515304	0.0509537	0.0511448
```

If you recognise a lower-triangle [distance matrix](https://evolution.gs.washington.edu/phylip/doc/distance.html) here, you are quite correct. We need to create a guide tree based on this. However, the program we'll use don't support lower-triangle matrices, but needs a full. So we need to convert it. I rewrite slightly a function I found on the [support forum to Mash](https://github.com/marbl/Mash/issues/9) to do this. Run these commands:

```
module purge
module load Biopython/1.72-foss-2018b-Python-3.6.6
python /cluster/projects/nn9458k/phylogenomics/week2/src/convert_triangle_full.py mash/infile  > full_distance_matrix
```

We need to load that Biopython module so we have the Python libraries that we need (numpy and pandas) available. The file full_distance_matrix should contain this now:
```
6
neobri	0.0	0.00981378	0.00985499	0.00904645	0.00924552	0.0519203
neogra	0.00981378	0.0	0.00888862	0.00948708	0.00924552	0.0509537
neomar	0.00985499	0.00888862	0.0	0.0090861	0.00932572	0.0515304
neooli	0.00904645	0.00948708	0.0090861	0.0	0.00659353	0.0509537
neopul	0.00924552	0.00924552	0.00932572	0.00659353	0.0	0.0511448
orenil	0.0519203	0.0509537	0.0515304	0.0509537	0.0511448	0.0
```

Then we'll run [rapidNJ](https://github.com/somme89/rapidNJ) on the full distance matrix to get a tree. Since this is also a quick process, we'll just run it on the command line and not in a sbatch script:
```
module purge
module load rapidNJ/210609-foss-2020b
rapidnj full_distance_matrix -i pd |sed "s/'//g"  nj_tree
```
We're removing the single quote sign because it would confuse Cactus later on.


The result is this:
```
((neomar:0.0045736,neogra:0.004315):0.00021445,((neopul:0.0033453,neooli:0.0032482):0.00094829,neobri:0.0049009):0.00032849,orenil:0.046583);
```
You can use [IcyTree](https://icytree.org/) (choose 'File' -> 'Enter tree directly') or [Phylogenetic tree (newick) viewer](http://etetoolkit.org/treeview/) (clear the alignment box and paste in the tree in the box on the left) to visualise it online.

When we are done with this part, we have both a set of softmasked files and a guide tree, so we are all set to run Cactus!

<a name="cactus"></a>
## Running Cactus
