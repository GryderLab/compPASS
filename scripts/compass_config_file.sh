#!/bin/bash

# folder location - should be in scripts folder of compPASS main folder
LOCATION="/mnt/vstor/SOM_GENE_BEG33/software/compPASS/"

# reference paths
DATASOURCE="/mnt/vstor/SOM_GENE_BEG33/ChIP_seq/hg38/DATA/"
SCRIPTSOURCE=$LOCATION"scripts/"
REFSOURCE=$LOCATION"references/"
OUTPUTDIR="/mnt/vstor/SOM_GENE_BEG33/ChIP_seq/hg38/projects/Pol2_ratios/compPASS_outs/"
GENEREF=$LOCATION"references/gene_2K_coord_hg38.tsv"
GCT=$LOCATION"references/gene_categories_binary.txt" # adds length, class, genomic, and more(!) for characterizing genes

# conda path
CONDAPATH="/mnt/vstor/SOM_GENE_BEG33/mamba/miniforge3/envs/" # < -- SET THIS VALUE
CONDAVAL=$CONDAPATH"compPASS" # < -- SET THIS VALUE

# regions of interest for Pol2 pileup
pro_start=-800
pro_end=-30
tssr_start=-30
tssr_end=300
gene_start=300 # added to actual_gene_end - actual_gene_start if sense
gene_end=0
tesr_start=0 # starts actual_gene_end
tesr_end=4000

# exclude chromosome Y? set to T if model is female-derived
EXCLUDE_CHR_Y=T

# use RNA-seq

# system config for downstream intensive processes (generating bigWig, matrix, profilePlots)
THREADCOUNT=8 # more preferred

# profile plot variables
# note: genes are plotted such that output always reads sense, left to right
leftgenebp="3000"  # extends past promoter region
genebodybp="5000"  # gene body regions scale, used for relative length
rightgenebp="4000" # matches tesr region
