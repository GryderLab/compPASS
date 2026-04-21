# Main R script for evaluating profile of chromatin reads
## from Pol2 ChIP/HiChIP data

# Yaw Asante | yaw.asante@case.edu
## January 15th, 2025 | Gryder Lab, CWRU


### Global Variables
NUMREGIONS <- 4 # number of ROIs for read capture
NUMCOLS <- 7 # number of total data columns for each sample

# reads column
SAMPLEA.COL <- 8
SAMPLEB.COL <- 9

# feature column
STRANDLOC <- 4
GENAMELOC <- 5

# header for 4 regions and 3 ratios
hdr <- c("Prom", "TSSR", "GB", "TESR", "StR", "PaR", "UnR")

# prevent sci-notation truncation
options(scipen = 999)

### Input Handling
args <- commandArgs(trailingOnly = TRUE)
project_dir <- args[1]
source_file_name <- args[2]
prescan_readcounts <-  args[3]

# global references for region lengths - defaults
.GlobalEnv$pro_start <- -800
.GlobalEnv$tssr_start <- -30
.GlobalEnv$gene_start <- 300
.GlobalEnv$tesr_end <- 4000

.GlobalEnv$pro_start <- as.numeric(args[4])
.GlobalEnv$tssr_start <- as.numeric(args[5])
.GlobalEnv$gene_start <- as.numeric(args[6])
.GlobalEnv$tesr_end <-  as.numeric(args[7])

script_source <- args[8]
exclude_chrY <- args[9]
output_name <- args[10]

if (exclude_chrY != "T" || exclude_chrY != "F"){
  exclude_chrY <- "F"
}
exclude_chrY <- as.logical(exclude_chrY)

# source functions
source(paste0(script_source, "/r_functions_compPASS.R"))

prescan_readcounts <- unlist(strsplit(prescan_readcounts, ","))

readcounts <- c(1, 1)

if (length(prescan_readcounts) < 2) {
  print("Only one scalefactor found, using no scalefactor ...")
  scalefactors <- c(1, 1)
} else {
  readcounts <- as.numeric(prescan_readcounts)
  scalefactors <- round(1000000 / readcounts, 3)
}


# divide counts by length of bin and scale by factor
# sample A as column 7, sample B as column 8
# capture each line of array by jumping 4 at a time

input.tb <- read.table(paste0(project_dir, "/", source_file_name),
                       sep = "\t", header = FALSE)

# sanity check, resort by chr, gene-name and coordinates
input.tb <- input.tb[order(input.tb[, 1], input.tb[, 5], input.tb[, 2]), ]

if (dim(input.tb)[1] %% NUMREGIONS != 0) {
print("Error, number of segment regions is not a multiple
of number of regions of interest, exiting!")
  stop()
}

# for every set of four regions, capture ratios and region changes
k <- 1
results.tb <- NULL
while(k <= dim(input.tb)[1]) {
  cur.rows <- input.tb[k:(k + (NUMREGIONS - 1)), ]
  cur_chr <- cur.rows[1, 1]
  if (exclude_chrY && cur_chr == "chrY") {
    k <- k + NUMREGIONS
    next;
  }

  strand <- cur.rows[1, STRANDLOC]
  
  # gene name
  gene_id <- cur.rows[1, GENAMELOC]
  gene_len <- -1
  
  row.contents <- c()
  if (strand == "+") {
    gene_len <- as.numeric(cur.rows[3, 3]) -
                as.numeric(cur.rows[3, 2]) +
                .GlobalEnv$gene_start

    row.contents.pre <- calculateVals(
                    cur.rows[1:NUMREGIONS, SAMPLEA.COL] * scalefactors[1],
                    cur.rows[1:NUMREGIONS, SAMPLEB.COL] * scalefactors[2],
                    gene_len
                    )
    row.contents <- as.numeric(row.contents.pre)
  } else {
    gene_len <- (as.numeric(cur.rows[2,3]) - as.numeric(cur.rows[2,2])) + .GlobalEnv$gene_start
    row.contents.pre <- calculateVals(
                      cur.rows[NUMREGIONS:1, SAMPLEA.COL] * scalefactors[1],
                      cur.rows[NUMREGIONS:1, SAMPLEB.COL] * scalefactors[2],
                      gene_len
                    )
    row.contents <- as.numeric(row.contents.pre)
  }
  
  results.tb <- rbind(results.tb,
                      c(gene_id, gene_len, as.numeric(row.contents)))
  k <- k + NUMREGIONS
}
results.tb <- as.data.frame(results.tb)
colnames(results.tb) <- c("GeneID", "GeneLength",
                          paste0(hdr, rep(1, NUMCOLS)), # sample 1 cols
						  paste0(hdr, rep(2, NUMCOLS)), # sample 2 cols
						  paste0(rep("L2FC_", NUMCOLS), hdr) # l2fc cols
						  )

# write out full table to output
write.table(x = results.tb, file = output_name,
            row.names = FALSE, col.names = TRUE,
			quote = FALSE, sep = "\t")

cat(" - Pol2 table made.\n")