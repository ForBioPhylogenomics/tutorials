# Requirements

Most analyses of the tutorials of the ForBio Phylogenomics course will be conducted using the [Saga server](https://documentation.sigma2.no/hpc_machines/saga.html) provided by [Sigma2](https://www.sigma2.no), where the required tools are already be installed and available through modules. There are a few exceptions, however, where the use of Graphical User Interface (GUI) programs will be required or beneficial; these should therefore be installed locally.

## Graphical User Interface programs

<!--* **AliView:** To visualize sequence alignments, the software [AliView](http://www.ormbunkar.se/aliview/) ([Larsson 2014](https://academic.oup.com/bioinformatics/article/30/22/3276/2391211)) is recommended. The installation of AliView is described at [http://www.ormbunkar.se/aliview/](http://www.ormbunkar.se/aliview/) and should be possible on all operating systems.
-->
* **PAUP\*:** A command-line version of [PAUP\*](http://paup.phylosolutions.com) will be installed on Saga; however, it may be easier to get to know the relevant functions of PAUP\* by using the GUI version of it. Unfortunately, this GUI version does not run on MacOS 10.15 (Catalina) or newer, but on other systems, the program can be installed using the instructions and precompiled versions available on [http://phylosolutions.com/paup-test/](http://phylosolutions.com/paup-test/).

* **FigTree:** The program [FigTree](http://tree.bio.ed.ac.uk/software/figtree/) by Andrew Rambaut is a very intuitive and useful tool for the visualization and (to a limited extent) manipulation of phylogenies encoded in [Newick](http://evolution.genetics.washington.edu/phylip/newicktree.html) format. Executables for Mac OS X, Linux, and Windows are provided on [https://github.com/rambaut/figtree/releases](https://github.com/rambaut/figtree/releases).

* **BEAST2:** The BEAST2 package, including BEAUti, BEAST2 itself, TreeAnnotator, and other tools can be downloaded from the BEAST2 website [https://www.beast2.org](https://www.beast2.org). As all these programs are written in Java, compilation is not required, and all programs should work on Mac OS X, Linux, and Windows. I recommend downloading, where possible, the program versions that include the Java Runtime Environment, which may prevent conflicts with Java versions that may already be installed on your machine.<br>

* **Tracer:** Just like BEAST2, Tracer is written in Java and should work on your system without problems. The program can be downloaded for Mac OS X, Linux, or Windows from [https://github.com/beast-dev/tracer/releases](https://github.com/beast-dev/tracer/releases). The file with the extension `.dmg` is for Mac OS X, the one with the extension `.tgz` is for Linux, and the Windows version is the file ending in `.zip`.