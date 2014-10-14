# Bash Exome Pipeline #

This repository collects scripts and documents useful to get up and running with NGS Target Region Sequencing.

### What is this repository for? ###

The main script implements the following steps:

*Verification of input integrity, quality checks, read trimming and primer contamination removal;
*BWA alignment;
*BAM conversion, sorting and indexing;
*Duplicates removal, as they result as PCR amplification bias;
*A local realignment around known IN-DELs position, carried on to delete the other artifacts;
*Quality score recalibration to refine some oddness caused by sequencing and mapping on quality scores;
*Variants (SNP and DIP) calling from the filtered mapping data obtained from the previous steps;

* Version Beta
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

### How do I get set up? ###

* Summary of set up
* Configuration
* Dependencies:
**GATK
**Samtools
* Deployment instructions

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact