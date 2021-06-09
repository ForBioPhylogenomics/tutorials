# Whole Genome Alignment

A tutorial on multiple whole genome alignment Cactus<br/>
By Ole K. Tørresen

## Summary

Large scale projects such as Earth BioGenome Project [Lewin et al. (2018)](https://doi.org/10.1073/pnas.1720115115), the Vertebrate Genomes Project [Rhie et al. (2021)](https://doi.org/10.1038/s41586-021-03451-0), and other similar efforts, aspire to sequence and assemble the genomes of all living species. More and more high quality genome assemblies are therefore being deposited in public databases, and these resources contain great potential and enables rigorous investigation into a multitude of aspects of biology through data - and curiosity driven research. The further generation of high-quality genomes will facilitate great leaps in our understanding of, and ability to answer, core questions within fundamental biology such as in ecology and evolution, and in other fields such as medicine.

However, to actually understand these genomes, we need some way of analysing them. For smaller projects, sequencing and assembling a few species, annotating and analysing each can be done focussing on each in turn. However, for large projects such as Zoonimia (131 new genome assemblies) [Zoonomia Consortium. (2020)](https://doi.org/10.1038/s41586-020-2876-6) or bird 10k (B10K; 267 new genome assemblies) [Feng et al. (2020)](https://doi.org/10.1038/s41586-020-2873-9), using that much time and effort on each genome is not scalable, and therefore they used whole-genome alignments to analyse their datasets. So far, the only whole-genome aligner that has been used successfully on such large datasets is Cactus [Armstrong et al. 2020](https://doi.org/10.1038/s41586-020-2871-y).

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

We will use RED, Mash, RapidNJ and Cactus. RED, Mash and RapidNJ should be available as modules on Saga, while Cactus works best as an container (the installation is non-trivial.)




<a name="softmask"></a>
## Softmasking with Red
