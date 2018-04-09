# AccNetPhylA

[This project is in development. Documentation is fairly light. You are welcomed to use this software, but please expect it to change in non-trivial ways. Feel free to return your impressions.]

## Contents 

* [Introduction](https://github.com/meb-team/AccNetPhylA/blob/master/README.md#introduction)
* [Installation](https://github.com/meb-team/AccNetPhylA/blob/master/README.md#installation)
* [How to use ?](https://github.com/meb-team/AccNetPhylA/blob/master/README.md#how-to-use)
* [Bugs](https://github.com/meb-team/AccNetPhylA/blob/master/README.md#bugs)
* [Citation](https://github.com/meb-team/AccNetPhylA/blob/master/README.md#citation)

## Introduction

This Perl script AccNetPhylA (Accessory Constellation Network with Phylogenetic Analysis) takes up most of the analysis of AccNET (Accessory Constellation Network) (Val F. Lanza et al., 2016) by replacing the clustering tool kClust (Hauser et al., 2013) with CD-HIT (Weizhong Li et al., 2001) for analysis with nucleotide and protein sequences. 

"AccNET is a comparative genomic tool for accessory genome analysis using bipartite networks. The software has been designed to be compatible with most of the Network Analysis software (i.e. Cytoscape, Gephi or R)."

## Installation 

### Prerequisites

 * [Perl](https://www.perl.org/) - scripting language
 * [CD-HIT](http://weizhongli-lab.org/cd-hit/) - Clustering program
 * [MUSCLE](https://www.drive5.com/muscle/) - Tool for MUltiple Sequence Comparison by Log-Expectation
 * [trimAl](http://trimal.cgenomics.org/) - Tool for removal of spurious sequences or poorly aligned regions
 * [PHYLIP](http://evolution.genetics.washington.edu/phylip.html) - PHYlogeny Inference Programs
 * [PhyML](http://www.atgc-montpellier.fr/phyml/) - Phylogeny software

### GitHub

Choose somewhere to put it, for example in your home directory (no root access required):

```bash
cd $HOME
```

Clone the latest version of the repository:

```bash
git clone https://github.com/meb-team/AccNetPhylA.git
```

## How to use ?

```
    Usage : perl AccNetPhylA.pl -in *.fasta -seq aa
    
            perl AccNetPhylA.pl [-in file.fasta] [-seq aa] [-cp '-c 0.9 -n 5 -g 1 -aS 0.9 -T 6 -M 5000'] 
	    [-out file1csv] [-tblout file2.csv] [-threshold 1] [-fast no] [-clean yes]


    Global options:

    -Help|help|h, produces this help file. 
	
    -in, Files to analyze as input.	
	
    -seq, Proteic or nucleic sequence (nt|aa).

    -cp, CD-HIT parameters. /!\ Output filename no required. Default('-c 0.9 -n 5 -g 1 -aS 0.9 -T 6 -M 5000')
	
    -out, Network filename. Default(Network.csv)
	
    -tblout, Table filename. Default(Table.csv)
	
    -threshold,	Percent of genomes to consider coregenome (values > 1 includes coregenome). Default(1)
	
    -fast, Skip the phylogenetic distance determination. Default(no)
	
    -clean, Remove and class files. Default(yes)
	
    -dir, Directory's name. Default(Analyse)
  
```


## Bugs

* Submit problems or requests here: https://github.com/meb-team/AccNetPhylA/issues


## Citation

### Authors
* Siguret Clea - [cleasiguret](https://github.com/cleasiguret)
