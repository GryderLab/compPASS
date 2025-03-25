#!/bin/bash

# reference paths
DATASOURCE="../input"
SCRIPTSOURCE="../scripts"
OUTPUTDIR="../outputs"
GENEREF="../references/gene_coord_hg38.tsv"
CONDAPATH="/mnt/vstor/SOM_GENE_BEG33/users/yxa181/mambaforge/envs/"
GCT="../references/gene_categories.txt" # adds length, class, genomic, and more(!) for characterizing genes

# regions of interest for Pol2 pileup
pro_start=-800
pro_end=-30
tssr_start=-30
tssr_end=300
gene_start=300
gene_end=0
tesr_start=0
tesr_end=4000

# cutoff for considering group of genes as a class
COUNTCUTOFF=20
# system config for downstream intensive processes (generating bigWig, matrix, profilePlots)
PROCESSORS= 1 # more preferred

# profile plot variables
