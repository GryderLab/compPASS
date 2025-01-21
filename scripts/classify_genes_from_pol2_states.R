# plot abundance, classify and chart out types of genes
library(dplyr)
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)

options(scipen=999)

Args <- commandArgs(trailingOnly=T)

project_dir <- Args[1]
prefiltered.tb <- Args[2]
path_to_genecat <- Args[3]
.GlobalEnv$gene_start <- Args[4]

#gene_cat_table <- read.table(path_to_genecat, sep="\t",header=T)
gene_cat_table <- read.table(path_to_genecat, sep="\t",header=T)

# filter table for genes of interest
filterAndLabelGenes <- function(pol2ratios, gene_cat_table) {
  # filter by genes included in categories table from reference
  ## Note: this autofilters out genes with length less than the TSSR_extension+50bp
  
  # filter 1: region length
  pol2ratios <- pol2ratios[pol2ratios$RegionLength > .GlobalEnv$gene_start+50,]
  
  colnames(gene_cat_table)[4] <- "GeneID"
  # filter 2: filter by genes included in gene_categories table
  in.genes <- gene_cat_table$GeneID
  pol2ratios <- pol2ratios[pol2ratios$GeneID %in% in.genes,]
  pol2ratios <- left_join(pol2ratios, gene_cat_table, by="GeneID")
  
  # filter 3: filter by amount of Pol2 in Promoter and TSSR regions of both
  .GlobalEnv$orig.ratios <- pol2ratios
  total_pol2.prom <- apply(pol2ratios[,c("Prom1","Prom2")],1, sum)/2
  total_pol2.tssr <- apply(pol2ratios[,c("TSSR1","TSSR2")],1, sum)/2
  quant_cut.prom <- quantile(total_pol2.prom, probs=seq(0,1,0.1))[3]
  quant_cut.tssr <- quantile(total_pol2.tssr, probs=seq(0,1,0.1))[3]
  p.filter <- total_pol2.prom > quant_cut.prom & total_pol2.tssr > quant_cut.tssr
  pol2ratios <- pol2ratios[p.filter,]
  
  .GlobalEnv$orig.ratios <- orig.ratios[!orig.ratios$GeneID %in% pol2ratios$GeneID,]
  return(pol2ratios)
}

filt.tb <- filterAndLabelGenes(input.synth, gene_cat_table) # reminder: clean gene table reference, duplicates not in refseqhg38
filt.tb <- filterAndLabelGenes(prefiltered.tb, gene_categories)

# get histograms of regions of interest across both samples

# get hexbin contour plots of promoter and tssr


threeTestCheck <- function(pol2ratios, cols, isHi){
  outlist <- apply(pol2ratios[, c(cols)], 1, function(x){
    test.a <- ifelse(isHi, x[1] > 0, x[1] < 0)
    test.b <- ifelse(isHi, x[2] < x[1], x[2] > x[1])
    test.c <- ifelse(isHi, isHiOutlier(x[3]), isLoOutlier(x[3]))
    return(test.a & test.b & test.c)
  })
  return(outlist)
}

isStalling <- function(pol2ratios) {
  threeTestCheck(pol2ratios, c("L2FC_Prom","L2FC_TSSR","L2FC_StR"), T)
}
isPausing <- function(pol2ratios) {
  threeTestCheck(pol2ratios, c("L2FC_TSSR","L2FC_GB","L2FC_PaR"), T)
}
isOverloading <- function(pol2ratios) {
  threeTestCheck(pol2ratios, c("L2FC_GB","L2FC_TESR","L2FC_UnR"), T)
}

isEntering <- function(pol2ratios) {
  threeTestCheck(pol2ratios, c("L2FC_Prom","L2FC_TSSR","L2FC_StR"), F)
}
isReleasing <- function(pol2ratios) {
  threeTestCheck(pol2ratios, c("L2FC_TSSR","L2FC_GB","L2FC_PaR"), F)
}
isUnloading <- function(pol2ratios) {
  threeTestCheck(pol2ratios, c("L2FC_GB","L2FC_TESR","L2FC_UnR"), F)
}