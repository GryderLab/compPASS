# Main R script for evaluating profile of chromatin reads from Pol2 ChIP/HiChIP data
# Yaw Asante | yxa181@case.edu | January 15th, 2025 | Gryder Lab, CWRU

# source functions
source("funct_compas.R")

### Global Variables
NUMREGIONS=4 # number of ROIs for read capture 

# reads column
SAMPLEA.COL <- 7
SAMPLEB.COL <- 8

# feature column
STRANDLOC <- 4
GENAMELOC <- 5

# header
hdr <- c("Prom","TSSR","GB","TESR","StR","PaR","UnR")

options(scipen = 999)

### Input Handling
Args <- commandArgs(trailingOnly=T)
project_dir <- Args[1]
source_file_name <- Args[2]
prescan_scalefactors <- Args[3]

# global references for region lengths
.GlobalEnv$pro_start=Args[4]
.GlobalEnv$tssr_start=Args[5]
.GlobalEnv$gene_start=Args[6]
.GlobalEnv$tesr_end=Args[7]

prescan_scalefactors <- unlist(strsplit(prescan_scalefactors, ","))

scalefactors <- c(1,1)
if(length(prescan_scalefactors) < 2){
	print("Only one scalefactor found, using no scalefactor ...")
} else if (sum(is.numeric(prescan_scalefactors)) < 2){
	print("Non-numeric scale factors found, using no scalefactor ...")
} else {
	scalefactors <- prescan_scalefactors
}

# divide counts by length of bin and scale by factor
# sample A as column 7, sample B as column 8
# capture each line of array by jumping 4 at a time

input.tb <- read.table(paste0(project_dir,"/",source_file_name), sep="\t",header=F)

# sanity check, resort by chr, gene-name and coordinates
input.tb <- input.tb[order(input.tb[,1], input.tb[,5], input.tb[,2]),]

if (dim(input.tb)[1] %% NUMREGIONS != 0){
	print("Error, number of segment regions is not a multiple of number of regions of interest, exiting!"
	stop()
}


k <- 1
results.tb <- NULL
while(k <= dim(input.tb)[1]){
	cur.rows <- input.tb[k:(k+(NUMREGIONS-1)),]
	strand <- cur.rows[k,NUMREGIONS]
	
	# gene name
	gene_id <- cur.rows[1,GENAMELOC]
	gene_len <- -1
	
	row.contents <- c()
	if (strand == "+") {
		gene_len <- (cur.rows[3,3] - cur.rows[3,2])+gene_start
		row.contents <- calculateVals(cur.rows[1:NUMREGIONS,SAMPLEA.COL]*scalefactors[1], cur.rows[1:NUMREGIONS,SAMPLEB.COL]*scalefactors[2], gene_len)
	} else {
		gene_len <- cur.rows[2,3] - cur.rows[2,2])+gene_start
		row.contents <- calculateVals(cur.rows[NUMREGIONS:1,SAMPLEA.COL]*scalefactors[1], cur.rows[NUMREGIONS:1,SAMPLEB.COL]*scalefactors[2], gene_len)
	}
	
	results.tb <- rbind(results.tb, c(gene_id, gene_len, row.contents))
	k <- k + NUMREGIONS
}
results.tb <- as.data.frame(results.tb)
colnames(results.tb) <- c(paste0(hdr, rep(1,NUMREGIONS)), paste0(hdr, rep(2,NUMREGIONS)), paste0(rep("L2FC_",NUMREGIONS), hdr))

# write out full table to output
write.table(results.tb, file=paste0(project_dir,"/comparison_tb.tsv", row.names=F, col.names=T, quote=F, sep="\t")


