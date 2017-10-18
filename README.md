# Bash Exome Pipeline #

This repository collects scripts and documents useful to get up and running with NGS Target Region Sequencing.

## What is this repository for ##

The main script implements the following steps:

* Verification of input integrity, quality checks, read trimming and primer contamination removal;
* BWA alignment;
* BAM conversion, sorting and indexing;
* Duplicates removal, as they result as PCR amplification bias;
* A local realignment around known IN-DELs position, carried on to delete the other artifacts;
* Quality score recalibration to refine some oddness caused by sequencing and mapping on quality scores;
* Variants (SNP and DIP) calling from the filtered mapping data obtained from the previous steps;

* Version Beta
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

### How do I get set up? ###

* Dependencies:

  * GATK
  * Picard tools
  * Samtools
  * BWA

* Ubuntu 16.04 LTS setup
  * GATK
    * dowload latest version from Broad site
    * unzip into your "tools" directory e.g. ~/bin in my case
    * create a symlink: ln -s GenomeAnalysisTK-3.x-y gatk
  * Picard tools
    * dowload latest version from Broad site
    * unzip into your "tools" directory e.g. ~/bin in my case
    * create a symlink: ln -s Picard_x.y picard
  * Samtools

    ln -s samtools-1.6 samtools
    cd samtools


    A few libs may be missing in your system. ./configure will complain about that, so add them:

    sudo apt-get install libncurses5 libncurses5-dev libncursesw5 libncursesw5-dev libbz2-1.0 libbz2-dev liblzma5 liblzma-dev

    ./configure

  * BWA
    * download latest version from: [Sourceforge](http://bio-bwa.sourceforge.net)

    * ln -s bwa_0.7.12/bwa.kit/ bwa
* Configuration


* Deployment instructions

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact