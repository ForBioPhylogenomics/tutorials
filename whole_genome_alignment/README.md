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
-rw-rw-r-- 1 olekto nn9244k   43497196 May 16 10:11 neobri.repeats.bed
-rw-rw-r-- 1 olekto nn9244k  275320088 May 16 10:11 neobri.repeats.fasta
-rw-rw-r-- 1 olekto nn9244k  865018698 Jun  5 01:10 neobri.softmasked.fa
-rw-rw-r-- 1 olekto nn9244k   45092622 May 16 13:37 neogra.repeats.bed
-rw-rw-r-- 1 olekto nn9244k  262196845 May 16 13:38 neogra.repeats.fasta
-rw-rw-r-- 1 olekto nn9244k  672846647 May 16 13:35 neogra.softmasked.fa
-rw-rw-r-- 1 olekto nn9244k   46494946 May 16 13:29 neomar.repeats.bed
-rw-rw-r-- 1 olekto nn9244k  269948250 May 16 13:29 neomar.repeats.fasta
-rw-rw-r-- 1 olekto nn9244k  683501275 May 16 13:27 neomar.softmasked.fa
-rw-rw-r-- 1 olekto nn9244k   46074867 May 16 13:50 neooli.repeats.bed
-rw-rw-r-- 1 olekto nn9244k  267874589 May 16 13:50 neooli.repeats.fasta
-rw-rw-r-- 1 olekto nn9244k  681893391 May 16 13:48 neooli.softmasked.fa
-rw-rw-r-- 1 olekto nn9244k   45076389 May 16 13:52 neopul.repeats.bed
-rw-rw-r-- 1 olekto nn9244k  262152533 May 16 13:52 neopul.repeats.fasta
-rw-rw-r-- 1 olekto nn9244k  674817636 May 16 13:49 neopul.softmasked.fa
-rw-rw-r-- 1 olekto nn9244k   30337323 May 16 10:24 orenil.repeats.bed
-rw-rw-r-- 1 olekto nn9244k  347698289 May 16 10:24 orenil.repeats.fasta
-rw-rw-r-- 1 olekto nn9244k 1025835685 Jun  5 01:09 orenil.softmasked.fa
```

<a name="mash"></a>
## Creating a "phylogeny" with Mash

In addition to softmasked genome assemblies, Cactus requires a guide tree. This is because it aligns two and two genomes. For instance, it could start with the two most closely related, aligns those, and then use that alignment to reconstruct an ancestral genome to those two. That reconstructed genome is then used to align against the next. This means that Cactus scales linear with the number of genomes and not quadratically as would have been the case if all genomes were compared against all. However, this approach might miss alignment of sequences where species A and C has it, and B not, but since A and B were closest related, those were aligned first and the reconstructed genome lacked the shared regions between A and C and therefore the resulting alignment didn't show that case.

Creating phylogenies can be a time consuming process, and therefore we will cheat (the "phylogeny" in the title is there for a reason). We will use Mash ([Ondov et al. (2016)](https://doi.org/10.1186/s13059-016-0997-x)), which is a fast genome and metagenome distance estimator.
