# Whole Genome Alignment

A tutorial on multiple whole genome alignment with Cactus<br>
By Ole K. Tørresen

## Summary

Large scale projects such as Earth BioGenome Project ([Lewin et al. (2018)](https://doi.org/10.1073/pnas.1720115115)), the Vertebrate Genomes Project ([Rhie et al. (2021)](https://doi.org/10.1038/s41586-021-03451-0)), and other similar efforts, aspire to sequence and assemble the genomes of all living species. More and more high quality genome assemblies are therefore being deposited in public databases, and these resources contain great potential and enables rigorous investigation into a multitude of aspects of biology through data - and curiosity driven research. The further generation of high-quality genomes will hopefully facilitate great leaps in our understanding of, and ability to answer, core questions within fundamental biology such as in ecology and evolution, and in other fields such as medicine.

However, to actually understand these genomes, we need some way of analysing them. For smaller projects, sequencing and assembling a few species, annotating and analysing each can be done focussing on each in turn. However, for large projects such as Zoonimia (131 new genome assemblies) ([Zoonomia Consortium (2020)](https://doi.org/10.1038/s41586-020-2876-6)) or bird 10k (B10K; 267 new genome assemblies) ([Feng et al. (2020)](https://doi.org/10.1038/s41586-020-2873-9)), using that much time and effort on each genome is not scalable, and therefore they used whole-genome alignments to analyse their datasets. So far, the only whole-genome aligner that has been used successfully on such large datasets is Cactus ([Armstrong et al. 2020](https://doi.org/10.1038/s41586-020-2871-y)).

## Table of contents

* [Outline](#outline)
* [Dataset](#dataset)
* [Requirements](#requirements)
* [Softmasking with Red](#softmask)
* [Creating a guide tree with Mash](#mash)
* [Running Cactus](#cactus)

<a name="outline"></a>
## Outline

In this tutorial we will run six genome assemblies of varying quality through this pipeline. This requires masking repeats in the genomes and a guide tree of the species involved. Unfortunately, most of these processes take longer than the time we have available in this session, but all intermediate and final files are available on Saga.

A short presentation of these concepts are [available here](whole_genome_alignment.pdf).

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

A detailed exposition of alignment is outside the scope of this tutorial (and I guess most of you know it quite well). However, to understand this step, we need to review some steps of it. Broadly, to create an alignment two sequences needs to be compared in some manner. The classic algorithms here are [Smith-Waterman](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm) and [Needleman-Wunsch](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm), for respectively local and global alignments. In local alignment the most similar sequence is found between the two input sequences, while a global alignment requires that the whole length of both input sequences are involved. For larger sequences (chromosomes for instance) or for multiple comparisons (between proteins and a genome), the requirements in computation (both memory and CPU) becomes unpractical. To remedy this, different approaches have been developed, such as [BLAST](https://en.wikipedia.org/wiki/BLAST_(biotechnology)) where words (substrings of the sequences) are connected to the sequences and only sequences with words in common are compared.

Repetitive sequences are sequences that are identical or very similar that occur multiple places in a genome. For instance the short tandem repeat ACT repeated both in tandem (ACTACTACTACT), but also multiple places, or different transposable elements that have copied themselves across the genome. These repetitive sequences will lead to large numbers of comparisons unless they are taken out of the picture in some way. The usual approach is to mask them, either hard mask (replace the repeats with Ns) or soft mask (lower case letters; standard fasta sequence is upper case). If a sequence is hard masked, it is basically gone from the analysis, so the recommendation is to soft mask. Alignments will not be initiated in soft masked regions, but alignments might extend into soft mask regions and therefore lead to less fragmented alignments.

Maybe most common today is to use [RepeatModeler](http://www.repeatmasker.org/RepeatModeler/) ([Flynn et al (2020)](https://doi.org/10.1073/pnas.1921046117)) to create a species specific repeat library, and use that in [RepeatMasker](https://www.repeatmasker.org/) to mask the genome. That process is a bit circumstantial so we'll take a short cut here (but feel free to use it in your own analysis). Running these two on a bird genome (about 1 Gbp) took around 20 hours with 40 CPUs for RepeatModeler and 2 hours with 40 CPUs for RepeatMasker, compared to 15 to 60 minutes with 1 CPU with using Red on the genomes here. RepeatModeler/RepeatMasker might be better, but it is hard to evaluate. If you are interested in the actual repeats themselves (which transposable element occurs where, for instance), running RepeatModeler/RepeatMasker is needed, but for just masking the repeats, Red is likely sufficient.

Red ([Giris (2015)](https://doi.org/10.1186/s12859-015-0654-5)) can do the softmasking in one step instead of two, and is faster. We'll run it via the Python script [redmask](https://github.com/nextgenusfs/redmask) to get readily softmasked genomes.

First, go to your user folder on Saga and copy the needed scripts (you might have done part of this before):

```
cd /cluster/projects/nn9458k/phylogenomics/
mkdir $YOURNAME
cd $YOURNAME
mkdir -p cichlids
cd cichlids
cp ../../week2/data/cichlids/scripts/* .

```

You can also look at the scripts in the folder [scripts](scripts), but it is easiest to get them on Saga.

The redmask command itself is quite straight forward:
```
redmask.py -i genome.fa -o mygenome
```

But we'll wrap in an SBATCH script and run it for all genomes. Unfortunately, redmask requires a Python module not readily available on Saga. To set install it locally (for you as a user), do this:
```
module purge
module load Biopython/1.72-foss-2018b-Python-3.6.6
pip3 install natsort --user
```
We load the same module that we use in the script to ensure that we install natsort for the correct version of Python.

Then you can do this:

```
sh start_red.sh
```
This little snippet starts one script per genome. Feel free to look around in it.

Be aware that this creates lots of small files which exceeded the quota for the course. run_red.sh now navigates to `$USERWORK` and runs the job there instead.

When I ran through this, it took 15 to 60 minutes per genome to finish.

If you are unable to get it running for some reason, or want to continue anyhow, you can do this:

```
mkdir -p masked_assemblies
cd masked_assemblies
rsync ../../../week2/data/cichlids/softmasked/* .
cd ..
```
When done, you should see something like this in masked_assemblies:
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
## Creating a guide tree with Mash

In addition to softmasked genome assemblies, Cactus requires a guide tree. This is because it aligns a handful of genomes at a time (2-5). For instance, it could start with the two most closely related, aligns those, and then use that alignment to reconstruct an ancestral genome to those two. That reconstructed genome is then used to align against the next. This means that Cactus scales linear with the number of genomes and not quadratically as would have been the case if all genomes were compared against all.

Creating phylogenies can be a time consuming process, and therefore we will cheat. We will use Mash ([Ondov et al. (2016)](https://doi.org/10.1186/s13059-016-0997-x)), which is a fast genome and metagenome distance estimator. Mash reduces large sequences to compressed sketch representations, enabling much quicker comparisons between sequences compared to older methods.

Since we use the original genome assembly fasta files, and not the softmasked ones, you can run this independent of the softmasking above.

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

If you recognise a lower-triangle [distance matrix](https://evolution.gs.washington.edu/phylip/doc/distance.html) here, you are quite correct. We need to create a guide tree based on this. However, the program we'll use don't support lower-triangle matrices, but needs a full. So we need to convert it. I rewrote slightly a function I found on the [support forum to Mash](https://github.com/marbl/Mash/issues/9) to do this. Run these commands:

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

Multiple genome alignment programs align multiple genomes toward each other, unlike pairwise genome alignment programs (such as Mummer and Minimap2) which only aligns two and two. They have existed for more than a decade, but several of the older were reference based, meaning that one genome was the reference and all others were aligned towards that. Any regions shared by the non-reference genomes would not be represented in the final alignment, a quite large limitation. [Cactus](https://github.com/ComparativeGenomicsToolkit/cactus), or Progressive Cactus as the authors also call it, is one of the few modern reference-free multiple genome aligners. The only other I know of is [SibeliaZ](https://github.com/medvedevgroup/SibeliaZ), but in my experience hasn't performed as well as Cactus.

Unfortunately, Cactus is not an easy program to use. There are likely multiple reasons for this, some might be that it parts are quite old (it builds directly on software older than a decade), which can lead to it being a bit clunky and not so easy to update those old parts. It is tightly connected to a system called Toil, which handles running jobs, and if that has some issues (we'll come back to that later), it can be hard to get working. For these reasons it can be difficult to install properly. Luckily they provide a Docker container, and Saga has support for that via Singularity, so we will use that.

Earlier, and when working with this tutorial, I spent quite some time getting Cactus to run successfully. I think I know how now, but we might have some issues which I'll try to explain how to avoid. This might be a bit technical, but I think it is good to know.

The first limitation is that [Toil](https://github.com/DataBiosphere/toil) uses a lot of hard links via [os.link](https://docs.python.org/3/library/os.html#os.link) (and not soft links) in Python to link files between folders, especially jobStore and workDir (overview of the jobs running and a temporary directory, respectively). Unfortunately, the filesystem on Saga, BeeGFS, [does not support hard linking between folders](https://www.beegfs.io/wiki/FAQ#hardlinks). That means that using --workDir as an option for Cactus will just not work, because files will be hard linked. You can also run into issues regarding temporary directories with this limitation. Some of this is discussed as [an issue on the Toil github page](https://github.com/DataBiosphere/toil/issues/2232).

The second is that Cactus uses /tmp quite a bit (unless specifying --workDir, but as explained above, doesn't work on Saga). The normal nodes on Saga seem to only have about 19 GB /tmp partitions, so larger jobs will quickly fill up these and then fail. Specifying a different temporary folder will not work, because that will be affected with the hard linking issue again. The bigmem nodes (specify #SBATCH --partition=bigmem in the sbatch script) have 303 GB /tmp partitions and will not as easily fill up those (but it only looks like the largest and not the medium sized ones have this).

(It might be possible to use $LOCALSCRATCH, see here: https://documentation.sigma2.no/files_storage/clusters.html#job-scratch-area-on-local-disk. I have tried it with Comparative Annotation Toolkit, and couldn't get it working, but the code is different between CAT and Cactus, so it might work. I haven't tested it yet.)

In the future it is likely that Toil will change their code to not hard link or to handle those cases better at least, or if you run on a different system than Saga, you might not have these limitations. However, I think it is good to point out that your and our analyses might stop or be delayed by such choices as which filesystem is used on our computation cluster (and that choice is usually without us and it would be hard to imagine all different cases such as this).

As you might understand, Cactus is a beast and therefore is quite computing hungry. When I ran the dataset, it used 30 hours using 32 CPUs on a bigmem node. We have been given a reservation for nodes on Saga, 128 CPUs from 13:15:00 and for 26 hours. All of the participants cannot run the full dataset, nothing will finish in time then. We'll have 2 people running the full dataset, with the rest running a reduced one, using only chromosome 5 from Nile tilapia and the corresponding sequences from the other species (see [Getting only chr5](#chr5) below if you'll like to see how I did that).

All cannot run this in the course folder, that is, in /cluster/projects/nn9458k/phylogenomics/$USERNAME, because that might exceed the limitations of what the course have been given. So we need to run this in $USERWORK. I have set this up in the scripts, and hopefully it will work fine.

To start, we need to grab the container. You don't have to do this, I put it into the work scripts, but this is basically the command:
```
singularity pull --name cactus-v1.3.0.sif docker://quay.io/comparative-genomics-toolkit/cactus:v1.3.0
```

If you copied the scrips in the first part, it is basically to just do:

```
sbatch run_cactus_chr5.sh
```

And everything should work fine. In the end you'll get cichlids_chr5.hal and cichlids_chr5.maf files in your /cluster/projects/nn9458k/phylogenomics/$USERNAME folder, in addition to halValidation.chr5.txt and halStats.chr5.txt, a validation and a statistics file of the HAL file. The whole process took about 5 hours when I tried it.

halValidation.chr5.txt should basically just contain:
```
File valid
```

While halStats.chr5.txt should contain this:
```
hal v2.1
((neomar:0.0045736,neogra:0.004315)Anc1:0.00021445,((neopul:0.0033453,neooli:0.0032482)Anc3:0.00094829,neobri:0.0049009)Anc2:0.00032849,orenil:0.046583)Anc0;

GenomeName, NumChildren, Length, NumSequences, NumTopSegments, NumBottomSegments
Anc0, 3, 28764338, 471, 0, 562986
Anc1, 2, 22457867, 1162, 482173, 155250
neomar, 0, 22119859, 1278, 152061, 0
neogra, 0, 21364358, 1287, 148253, 0
Anc2, 2, 29214047, 356, 555446, 254028
Anc3, 2, 22168327, 1182, 133108, 140058
neopul, 0, 20898136, 1358, 133718, 0
neooli, 0, 21525845, 1293, 136975, 0
neobri, 0, 43197405, 18, 315265, 0
orenil, 0, 39714817, 1, 692510, 0
```
Running the full alignment is quite similar to the one above, but we need some extra information for the sbatch command since we have gotten two bigmem nodes reserved.

```
sbatch --reservation=nn9458k run_cactus_all.sh
```
If that is successful, you'll see cichlids_all.hal, halValidation.all.txt and halStats.all.txt with the content of the stats file this:

```
hal v2.1
((neomar:0.0045736,neogra:0.004315)Anc1:0.00021445,((neopul:0.0033453,neooli:0.0032482)Anc3:0.00094829,neobri:0.0049009)Anc2:0.00032849,orenil:0.046583)Anc0;

GenomeName, NumChildren, Length, NumSequences, NumTopSegments, NumBottomSegments
Anc0, 3, 681329676, 32958, 0, 16672331
Anc1, 2, 659155609, 75358, 15727138, 6780370
neomar, 0, 668470958, 89823, 6854395, 0
neogra, 0, 657998073, 91312, 6746422, 0
Anc2, 2, 662754018, 18173, 15927409, 5320178
Anc3, 2, 660323875, 78298, 4933681, 6370077
neopul, 0, 659878940, 94163, 6380869, 0
neooli, 0, 666859030, 91778, 6457234, 0
neobri, 0, 847910432, 9099, 6072207, 0
orenil, 0, 1005681550, 2460, 22880110, 0
```

You'll work further with the cichlids_chr5.maf file tomorrow.

One of the main advantages of having a multiple whole genome alignment is that you can pinpoint conserved sequences across the species. In many cases that will be exons, and a HAL file can be use for [comparative annotation](https://github.com/ComparativeGenomicsToolkit/Comparative-Annotation-Toolkit). This is outside the scope of this tutorial.

You have come to the end and are done with the tutorial. Congratulations!

<a name="chr5"></a>
## Getting only chr5
This is for information purposes only, and to get this kinda complete.

To create a reduced dataset we'll only use chromosome 5 from Nile tilapia and the corresponding sequences from the other species. For each species I did this:
```
module load StdEnv minimap2/2.17-GCC-8.3.0

mkdir -p paf_alignments

cd paf_alignments

if [ ${1} != 'orenil' ]; then
	minimap2 -t 10  -cx asm20 ../../data/orenil.fasta ../../data/${1}.fasta > $1_orenil.asm20.paf
fi
```
That is, mapped all against Nile tilapia. Chromosome 5 is called NC_031970.2 in Nile tilapia, so we extract all sequences that map to it at a quality higher than 50 and with an alignment length longer than 5000 bp:
```
module load StdEnv seqtk/1.3-foss-2018b

for i in $(ls ../data/*fasta); do
	j=${i%.fasta}
	l=${j##*/}
	echo $l
	#>NC_031970.2 is LG5/chromosome 5
	#quality of more than 50 and longer than 5000 bp alignmentx
	grep NC_031970.2 paf_alignments/${l}_orenil.asm20.paf |awk '$12 > 50' |awk '$10 > 5000'| cut -f 1 | sort -u > paf_alignments/${l}_maps_to_chr5
	seqtk subseq masked_assemblies/${l}.softmasked.fa paf_alignments/${l}_maps_to_chr5 > masked_assemblies/${l}.softmasked.chr5.fa
done

samtools faidx  masked_assemblies/orenil.softmasked.fa NC_031970.2 > masked_assemblies/orenil.softmasked.chr5.fa
```
In the end the *softmasked.chr5.fa files were created.

If you are reading this and thinking "Didn't you basically do what we did above in a fraction of the time?", then yes, you are correct. There might be some more alignments done inside Cactus, but it is clearly not the quickest whole genome aligner out there. Minimap2 is much faster, but it is only pairwise, and there are presently no way of using Minimap2 instead of LastZ (which is in Cactus) now. However, [Cactus does have a pangenome pipeline implemented](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md), which actually uses Minimap2 but assumes genomes from one species and they need to do some tricks (splitting up into chromosomes) to get it to work.
