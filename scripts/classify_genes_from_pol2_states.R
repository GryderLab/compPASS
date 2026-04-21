# plot abundance, classify and chart out types of genes

# functions

# cutoffs for significant change defined in Asante et al
## from triplicate Pol2 ChIP-seq in PDAC cells
PROMOTER_L2FC_CUTOFF <- 0.799
TSSR_L2FC_CUTOFF <- 0.799
GENEBODY_L2FC_CUTOFF <- 0.399
TESR_L2FC_CUTOFF <- 0.799

# L2FC for ratios determined to be 1
## 95% of PDAC genes had L2FC < 1
RATIO_L2FC_CUTOFF <- 0.999

# Fold difference expected of starting
RATIO_THRESHOLD <- 1.999

ratioGainCheck <- function(pol2ratios, cols, isHi, cutoff){
  outlist <- apply(pol2ratios[, c(cols)], 1, function(x){
    numGain <- x[1] > cutoff
    denomLoss <- x[2] < 0
    ratioGain <- x[3] > RATIO_L2FC_CUTOFF
    postRatio <- x[4] >= RATIO_THRESHOLD

    result <- numGain & denomLoss & ratioGain & postRatio
    return(result)
  })
  return(outlist)
}

ratioLossCheck <- function(pol2ratios, cols, isHi, cutoff){
  outlist <- apply(pol2ratios[, c(cols)], 1, function(x){
    numLoss <- x[1] < 0
    denomGain <- x[2] > cutoff
    ratioLoss <- x[3] < -(RATIO_L2FC_CUTOFF)
    preRatio <-  x[4] < RATIO_THRESHOLD

    result <- numLoss & denomGain & ratioLoss & preRatio
    return(result)
  })
  return(outlist)
}

isStalling <- function(pol2ratios) {
  ratioGainCheck(pol2ratios, c("L2FC_Prom", "L2FC_TSSR", "L2FC_StR", "StR2"),
                 TRUE, PROMOTER_L2FC_CUTOFF)
}

isPausing <- function(pol2ratios) {
  ratioGainCheck(pol2ratios, c("L2FC_TSSR", "L2FC_GB", "L2FC_PaR", "PaR2"),
                 TRUE, TSSR_L2FC_CUTOFF)
}

isClogging <- function(pol2ratios) {
  ratioGainCheck(pol2ratios, c("L2FC_GB", "L2FC_TESR", "L2FC_UnR", "UnR2"),
                 TRUE, GENEBODY_L2FC_CUTOFF)
}

isEntering <- function(pol2ratios) {
  ratioLossCheck(pol2ratios, c("L2FC_Prom", "L2FC_TSSR", "L2FC_StR", "StR1"),
                 FALSE, PROMOTER_L2FC_CUTOFF)
}

isReleasing <- function(pol2ratios) {
  ratioLossCheck(pol2ratios, c("L2FC_TSSR", "L2FC_GB", "L2FC_PaR", "PaR1"),
                 FALSE, TSSR_L2FC_CUTOFF)
}

isUnloading <- function(pol2ratios) {
  ratioLossCheck(pol2ratios, c("L2FC_GB", "L2FC_TESR", "L2FC_UnR", "UnR1"),
                 FALSE, GENEBODY_L2FC_CUTOFF)
}



args <- commandArgs(trailingOnly = TRUE)

project_dir <- args[1]
ratio_table <- args[2]
script_source <- args[3]
path_to_genecat <- args[4]
ctrl_name <- args[5]
case_name <- args[6]
output_name <- args[7]

path_to_ratio_table <- paste0(project_dir, "/", ratio_table)

# parse inputs
source(paste0(script_source,"/r_functions_compPASS.R"))
prefiltered.tb <- read.table(path_to_ratio_table, sep = "\t", header = TRUE)
gene_cat_table <- read.table(path_to_genecat, sep = "\t", header = TRUE)
curname <- paste0(ctrl_name, "_vs_", case_name)

# filter table for genes of interest
filterAndLabelGenes <- function(pol2ratios, gene_cat_table) {

  # filter by amount of Pol2 in Promoter and TSSR regions of both
  .GlobalEnv$orig.ratios <- pol2ratios
  total_pol2.prom <- apply(pol2ratios[, c("Prom1","Prom2")], 1, sum) / 2
  total_pol2.tssr <- apply(pol2ratios[, c("TSSR1","TSSR2")], 1, sum) / 2
  quant_cut.prom <- quantile(total_pol2.prom, probs=seq(0, 1, 0.1))[4]
  quant_cut.tssr <- quantile(total_pol2.tssr, probs=seq(0, 1, 0.1))[4]

  p.filter <- total_pol2.prom > quant_cut.prom &
              total_pol2.tssr > quant_cut.tssr 
  pol2ratios <- pol2ratios[p.filter, ]
  
  #.GlobalEnv$orig.ratios <- orig.ratios[!orig.ratios$GeneID %in% pol2ratios$GeneID,]
  return(pol2ratios)
}

# reminder: clean gene table reference, duplicates not in refseqhg38
pol2ratios <- filterAndLabelGenes(prefiltered.tb, gene_cat_table)

# additional variables
#pol2ratios$DProm <- pol2ratios$Prom2 - pol2ratios$Prom1
#pol2ratios$DTSSR <- pol2ratios$TSSR2 - pol2ratios$TSSR1
#pol2ratios$DGB <- pol2ratios$GB2 - pol2ratios$GB1
#pol2ratios$DTESR <- pol2ratios$TESR2 - pol2ratios$TESR1
samp1.main_cols <- c("Prom1", "TSSR1", "GB1", "TESR1")

pol2ratios$TotalPol2.samp1 <- apply(pol2ratios[, samp1.main_cols], 1, sum)

samp2.main_cols <- c("Prom2", "TSSR2", "GB2", "TESR2")
pol2ratios$TotalPol2.samp2 <- apply(pol2ratios[, samp2.main_cols], 1, sum)



# get histograms of regions of interest across both samples
density_plot.name <- paste0(project_dir, "/", curname, "_density_plot.pdf")
out.plot <- makeDensityDistributionPlots(pol2ratios)
ggsave(filename= density_plot.name, out.plot, width = 12, height = 6)

# get hexbin contour plots of promoter and tssr
q1.x <- round(quantile(log2(prefiltered.tb$Prom1 + prefiltered.tb$TSSR1 + 0.1),
                       probs = seq(0, 1, 0.05))[20], 3)
q1.y <- round(quantile(log2(prefiltered.tb$Prom2 + prefiltered.tb$TSSR2 + 0.1),
                       probs = seq(0, 1, 0.05))[20], 3)

q2.x <- round(quantile(log2(pol2ratios$Prom1+pol2ratios$TSSR1 + 0.1),
                       probs = seq(0, 1, 0.05))[20], 3)
q2.y <- round(quantile(log2(pol2ratios$Prom2+pol2ratios$TSSR2 + 0.1),
                       probs = seq(0, 1, 0.05))[20], 3)

max.x <- max(q1.x, q2.x)
max.y <- max(q1.y, q2.y)
max.o <- ceiling(max(max.x, max.y))
min_value_hex <- -4
TEST_MAX_BINVAL <- 500

pre.filt.hex <- ggplot(prefiltered.tb) +
  geom_hex(aes(x=log2(Prom1 + TSSR1 + 0.1), y = log2(Prom2+TSSR2+0.1)), bins = 50) +
  scale_x_continuous(paste0("Log2 ", ctrl_name," Prom + TSSR RPM"), limits=c(min_value_hex, max.o))+
  scale_y_continuous(paste0("Log2 ", case_name," Prom + TSSR RPM"), limits=c(min_value_hex, max.o))+
  scale_fill_gradient(limits = c(0, TEST_MAX_BINVAL),
                      low = "black", high = "cyan",
                      #oob = scales::oob_squish,
                      na.value = "white" ) +
  ggtitle("Pre-Filtering Distribution")+
  theme_light()
  
postfilt.hex <- ggplot(pol2ratios)+
  geom_hex(aes(x=log2(Prom1+TSSR1+0.1), y=log2(Prom2+TSSR2+0.1)), bins=50)+
  scale_x_continuous(paste0("Log2 ", ctrl_name," Prom + TSSR RPM"), limits=c(min_value_hex, max.o))+
  scale_y_continuous(paste0("Log2 ", case_name," Prom + TSSR RPM"), limits=c(min_value_hex, max.o))+
  scale_fill_gradient(limits = c(0, TEST_MAX_BINVAL),
                      low = "black", high = "cyan",
                      #oob = scales::oob_squish,
                      na.value = "white" ) +
  ggtitle("Post-Filtering Distribution")+
  theme_light()

out_hex <- plot_grid(pre.filt.hex, postfilt.hex, nrow=1, ncol=2)
ggsave(filename=paste0(project_dir,"/",curname,"_hex_plots.pdf"), out_hex, width=12, height=6)

# plot distribution of reads as boxplots, by regions
#curr.regions <- c("Prom","TSSR","GB","TESR")
#boxplts_list <- list()
#i <- 1
#for (val in curr.regions){
#	subtb <- pol2ratios[,c("GeneID",paste0(val,"1"),paste0(val,"2"))]
#	colnames(subtb) = c("GeneID", ctrl_name, case_name)
#	dist.pivot <- pivot_longer(subtb, names_to="Sample", values_to="RPM", !"GeneID")
#	
#	label_order <- factor(ctrl_name, case_name)
#	outplot <- ggplot(dist.pivot,aes(x=factor(Sample, levels=label_order),
#	            y=log2(RPM+0.1),fill=Sample))+
#				geom_boxplot()+
#				theme_light()+
#				scale_x_discrete("")+
#				scale_y_continuous("Log2 RPM Density")+
#				ggtitle(paste0(val," Distribution"))+
#				theme(legend.position="none")
#	boxplts_list[[i]] <- outplot
#	i <- i +1
#}
#out_dists = plot_grid(plotlist=boxplts_list, ncol=4, nrow=1)
#ggsave(filename=paste0(project_dir,"/",curname,"_box_dist_plots.pdf"), out_dists, width=8, height=4)


# assign status based on ratios
depleted_genes <- prefiltered.tb$GeneID[! prefiltered.tb$GeneID %in% pol2ratios$GeneID]
  
# genes where every segment has  loss
loss.positions <- (log2(pol2ratios$Prom2 + pol2ratios$TSSR2 +
                        pol2ratios$GB2 + pol2ratios$TESR2) -
                   log2(pol2ratios$Prom1 + pol2ratios$TSSR1 +
                        pol2ratios$GB1 + pol2ratios$TESR1)) < (-1)
total_loss_genes <- pol2ratios$GeneID[loss.positions]
  
# genes where every segment has a major gain
gain.positions <- (log2(pol2ratios$Prom2 + pol2ratios$TSSR2 +
                        pol2ratios$GB2 + pol2ratios$TESR2) -
                   log2(pol2ratios$Prom1 + pol2ratios$TSSR1 +
                        pol2ratios$GB1 + pol2ratios$TESR1)) > 1

total_gain_genes <- pol2ratios$GeneID[gain.positions]

main.classes <- c("Stall", "Entry",
                  "Pause", "Release",
                  "Clog", "Unload",
                  "Gain", "Loss")

pol2ratios$stall <- ifelse(isStalling(pol2ratios), "Stall","")
pol2ratios$pause <- ifelse(isPausing(pol2ratios)  , "Pause","")
pol2ratios$clog <-  ifelse(isClogging(pol2ratios) , "Clog","")
  
pol2ratios$entry <-   ifelse(isEntering(pol2ratios), "Entry","")
pol2ratios$release <- ifelse(isReleasing(pol2ratios) , "Release","")
pol2ratios$unload <-  ifelse(isUnloading(pol2ratios) , "Unload","")
premergedclass <- pol2ratios |> dplyr::select(stall,entry,pause,release,clog,unload) |> 
  apply(1, paste,collapse="")
  
pol2ratios$gain <- ifelse(premergedclass == "" & pol2ratios$GeneID %in% total_gain_genes,
                          "Gain","")
pol2ratios$loss <- ifelse(premergedclass == "" & pol2ratios$GeneID %in% total_loss_genes,
                          "Loss","")
pol2ratios$mergedclass <- pol2ratios |>
                          dplyr::select(stall,entry,
                                        pause,release,
                                        clog,unload,
                                        gain,loss) |> 
                          apply(1, paste,collapse="")

pol2ratios$mergedclass[pol2ratios$mergedclass == ""] <- "None"

# make ratio plots
se.plot <- makeRatioPlot(pol2ratios, "stall","entry","StR","L2FC_Prom","L2FC_TSSR",
                         "L2FC Promoter Reads","L2FC TSSR Reads",
                         "Stalling","Entering", curname)
pr.plot <- makeRatioPlot(pol2ratios, "pause","release","PaR","L2FC_TSSR","L2FC_GB",
                         "L2FC TSSR Reads","L2FC Gene Body Reads",
                         "Pausing","Releasing", curname)
cu.plot <- makeRatioPlot(pol2ratios, "clog","unload","UnR", "L2FC_GB","L2FC_TESR",
                         "L2FC Gene Body Reads","L2FC TESR Reads",
                         "Clogging","Unloading", curname)
  
gout <- plot_grid(se.plot, pr.plot, cu.plot, ncol = 3, nrow = 1)
ggsave(filename = paste0(project_dir, "/", curname, "_ratio_plots.pdf"),
       gout, width = 12, height = 5)

# make MA plots for gain and loss

# don't rescale since values not scaled by density

#pol2ratios$TotalPol2.samp1 <- apply(pol2ratios[,2:6], 1,  function(x){
#  gene_len <- x[2]
#  
#  prom <- x[3] * (.GlobalEnv$region_tssr_start - .GlobalEnv$region_pro_start + 1)
#  tssr <- x[4] * (.GlobalEnv$region_gene_start - .GlobalEnv$region_tssr_start + 1)
#  gb   <- x[5]*gene_len
#  tesr <- x[6] * (.GlobalEnv$region_tesr_end + 1)
#  totpol2 <- prom + tssr + gb + tesr
#  return(totpol2/ (gene_len+.GlobalEnv$region_tesr_end - .GlobalEnv$region_pro_start))
#})
#pol2ratios$TotalPol2.samp2 <- apply(pol2ratios[,c(2,10:13)], 1,   function(x){
#  gene_len <- x[2]
#  
#  prom <- x[3] * (.GlobalEnv$region_tssr_start - .GlobalEnv$region_pro_start + 1)
#  tssr <- x[4] * (.GlobalEnv$region_gene_start - .GlobalEnv$region_tssr_start + 1)
#  gb   <- x[5]*gene_len
#  tesr <- x[6] * (.GlobalEnv$region_tesr_end + 1)
#  totpol2 <- prom + tssr + gb + tesr
#  return(totpol2/ (gene_len+.GlobalEnv$region_tesr_end - .GlobalEnv$region_pro_start))
#})


pol2ratios$Delta_TotPol2 <- round(
                            pol2ratios$TotalPol2.samp2 -
                            pol2ratios$TotalPol2.samp1, 
                            3)

pol2ratios$L2FC_TotPol2 <- round(
                           log2(pol2ratios$TotalPol2.samp2 + 0.1) -
                           log2(pol2ratios$TotalPol2.samp1 + 0.1),
                           3)

gain_ratios <- subset(pol2ratios, pol2ratios$mergedclass == "Gain")
loss_ratios <- subset(pol2ratios, pol2ratios$mergedclass == "Loss")
ma.plot <- ggplot() +
          geom_point(data=pol2ratios,
                     aes(x = Delta_TotPol2, y = L2FC_TotPol2),
                     color = "grey") +
            geom_point(data = gain_ratios,
                      aes(x = Delta_TotPol2, y = L2FC_TotPol2),
                      color="red") +
            geom_point(data = loss_ratios,
                      aes(x = Delta_TotPol2, y = L2FC_TotPol2),
                      color = "blue") +
            scale_x_continuous("Delta Total Gene Pol2 (RPM)") +
            scale_y_continuous(paste0("Log2 ",case_name, " / ",
                               ctrl_name, " Total Gene Pol2")) +
            ggtitle(paste0("Gain (red) and Loss (blue) in Total Pol2
                            \nOverall Gain (n = ",
                    nrow(gain_ratios),
            ") and Overall Loss (n = ", nrow(loss_ratios),")")) +
            geom_hline(yintercept=-0.99, linetype="dashed") +
            geom_hline(yintercept=0.99, linetype="dashed") +
            theme_light()
ggsave(filename = paste0(project_dir, "/", curname,"_gain_loss_plots.pdf"),
       ma.plot, width = 5, height = 4)


# replot but standardized to fit window
ma.max_y <- max(abs(summary(pol2ratios$L2FC_TotPol2)))
ma.max_y <- ceiling(ma.max_y / 0.5) * 0.5 # to nearest

ma.max_x <- ceiling(log2(max(abs(pol2ratios$Delta_TotPol2))))
ma.plot2 <- ggplot() +
            geom_point(data = pol2ratios,
                      aes(x = log2(abs(Delta_TotPol2)), y = L2FC_TotPol2),
                      color="grey") +
            geom_point(data = gain_ratios,
                      aes(x = log2(abs(Delta_TotPol2)), y = L2FC_TotPol2),
                      color="red") +
            geom_point(data = loss_ratios,
                       aes(x = log2(abs(Delta_TotPol2)), y = L2FC_TotPol2),
                       color="blue") +
            scale_x_continuous("Log2 Abs Diff Total Gene Pol2 (RPM)",
            limits=c(0, ma.max_x)) +
            scale_y_continuous(paste0("Log2 ", case_name, " / ",
                               ctrl_name, " Total Gene Pol2"),
                               limits=c(-ma.max_y, ma.max_y)) +
            ggtitle(paste0("Gain (red) and Loss (blue) in Total Pol2
                            \nOverall Gain (n = ",
                    nrow(gain_ratios),
                    ") and Overall Loss (n = ",
                    nrow(loss_ratios),")")) +
            theme_light()
ggsave(filename = paste0(project_dir, "/", curname, "_gain_loss_plots_fitted.pdf"),
       ma.plot2, width=5, height=4)
# do not plot mixed modes
#upset_plot.name <- paste0(project_dir, "/", curname,"_upset_plot.pdf")
#out_upset <- makeUpsetPlot(pol2ratios, 0)

#pdf(file=upset_plot.name, width=8, height=6)
#print(out_upset)
#invisible(dev.off())
#ggsave(filename=upset_plot.name, out_upset, width=8, height=6)

### Create final ratios table
# assign final classes based on greatest ratio
final.classes <- c()
index <- 1
while(index <= nrow(pol2ratios)) {
  row.pol2ratio <- pol2ratios[index,]
  if (!row.pol2ratio$mergedclass %in% c(main.classes, "None")) {
    str.val <- row.pol2ratio$L2FC_StR
    par.val <- row.pol2ratio$L2FC_PaR
    unr.val <- row.pol2ratio$L2FC_UnR
    val.vector <- c(str.val, par.val, unr.val)

    get.col <- which(abs(val.vector) == max(abs(val.vector)))
    select.val <- val.vector[get.col]
    if(select.val == str.val) {
      final.classes <- append(final.classes,
                     ifelse (select.val > 0, "Stall", "Entry"))
    } else if(select.val == par.val) {
      final.classes <- append(final.classes,
                     ifelse (select.val > 0, "Pause", "Release"))
    } else {
      final.classes <- append(final.classes,
                     ifelse (select.val > 0, "Clog", "Unload"))
    }
  } else {
    final.classes <- append(final.classes, row.pol2ratio$mergedclass)
  }
  index <- index + 1
}
pol2ratios$ratioclass <- final.classes

# output table
final.ratios <- pol2ratios[, c("GeneID","GeneLength",
                                "TotalPol2.samp1", "TotalPol2.samp2",
                                "L2FC_StR", "L2FC_PaR", "L2FC_UnR",
                                "Delta_TotPol2", "L2FC_TotPol2",
                                "ratioclass")]

## ---- robust write checks ----
# ensure project dir and gene-lists dir exist
if (!dir.exists(project_dir)) {
  dir.create(project_dir, recursive = TRUE)
}
outdir2 <- paste0(project_dir, "/", curname, "_genelists")
if (!dir.exists(outdir2)) {
  dir.create(outdir2, recursive = TRUE)
}

# ensure ratioclass length matches rows (catch earlier bugs)
if (exists("final.classes")) {
  if (length(final.classes) != nrow(pol2ratios)) {
    stop("Length of final.classes (", length(final.classes),
         ") does not match nrow(pol2ratios) (", nrow(pol2ratios), ").")
  }
}

# coerce final table to plain data.frame and show structure
final.ratios <- as.data.frame(final.ratios)
#message("About to write final.ratios with ", nrow(final.ratios), " rows and ",
#        ncol(final.ratios), " cols. Column classes:")
#print(sapply(final.ratios, function(x) class(x)[1]))

# write with full path (safer)
out_final_path <- file.path(output_name)
tryCatch({
  write.table(x = final.ratios, file = out_final_path, sep="\t",
              row.names = FALSE, col.names = TRUE, quote = FALSE)
  message(" - Wrote final ratios to: ", out_final_path)
}, error = function(e) {
  stop("Failed to write final.ratios: ", conditionMessage(e))
})

write.table(file = output_name, x=final.ratios,
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


columns_of_interest <- c("Prom","TSSR","GB","TESR","StR","PaR","UnR")
col_select <- c(paste0(columns_of_interest, "1"), 
				paste0(columns_of_interest, "2"), 
				paste0("L2FC_", columns_of_interest),
				"ratioclass")

write.table(file=paste0(project_dir, "/full_table_",
                        curname, "_pol2_ratios.tsv"),
			      x=pol2ratios,#[, c("GeneID", "GeneLength", col_select)],
            sep = "\t", row.names = FALSE,
            col.names = TRUE, quote = FALSE)

cat(" - Gene classifications made.\n")


### plot bar plot of final ratios
pol2ratios.count <- table(pol2ratios$mergedclass[pol2ratios$mergedclass %in% main.classes])
df.count <- data.frame(class = names(pol2ratios.count),
                       counts = as.numeric(pol2ratios.count))

bar.out <- ggplot(df.count)+
  geom_col(aes(x=factor(class, levels = main.classes),
           y = counts), fill = "blue")+
  scale_y_continuous("Gene Counts") +
  scale_x_discrete("Ratio Class") +
  theme_light()
ggsave(filename = paste0(project_dir, "/", curname, "_barplot.pdf"),
                       bar.out, width = 5, height = 6)

#bar.max_y <- ceiling(log10(max(df.count$counts)))
bar.out2 <- ggplot(df.count)+
  geom_col(aes(x = factor(class, levels = main.classes),
               y = counts), fill = "blue")+
  scale_y_log10("Gene Counts") +
  scale_x_discrete("Ratio Class") +
  theme_light()
ggsave(filename = paste0(project_dir, "/", curname, "_barplot_logscaled.pdf"),
                       bar.out2, width = 5, height = 6)


total.genelist <- (read.table(paste0(script_source, "/gene_2K_coord_hg38.tsv"),
                    sep = "\t", header = TRUE))$GENE


### Donut Plot
gene.donut.tb <- data.frame(GeneID = total.genelist, start.class = "Depleted")

gene.donut.tb$start.class[gene.donut.tb$GeneID %in% pol2ratios$GeneID] <- "No Major Shift"

perturbed.list <- pol2ratios$GeneID[pol2ratios$ratioclass != "None"]
gene.donut.tb$start.class[gene.donut.tb$GeneID %in% perturbed.list] <- "Shifted"

temp.tb <- table(gene.donut.tb$start.class)

donut.vals <- data.frame(class = names(temp.tb),
                       perc = round(as.numeric(temp.tb) /
                              length(total.genelist), 2))

donut.vals$ymin <- c(0, donut.vals$perc[1],
                     donut.vals$perc[1] + donut.vals$perc[2])
donut.vals$ymax <- c(donut.vals$perc[1],
                     donut.vals$perc[1] + donut.vals$perc[2], 1)

# Set colors
donut.colors <- c("Depleted" = "black",
                  "No Major Shift" = "gray",
                  "Shifted" = "blue")

# Create donut plot
donut.plot <- ggplot(donut.vals) +
  geom_rect(aes(ymax = ymax, ymin = ymin,
                xmax = 4, xmin = 2.5,
                fill = class), color = "white") +
  geom_text(aes(x = 3.25, y = (ymin + ymax)/2,
                label = paste0(round(perc * 100, 1), "%")),
            color = "white", size=5, fontface = "bold") +
  scale_fill_manual(values = donut.colors) +
  coord_polar(theta = "y") +
  xlim(c(1, 4)) +  # donut width
  theme_void()

ggsave(filename=paste0(project_dir, "/", curname, "_donutplot.pdf"),
       donut.plot, width = 5, height = 4)


### output genelists
# uses cutoff of 30 as built-in
outdir2 <-  paste0(project_dir, "/", curname, "_genelists")
for(cur.class in pol2ratios$ratioclass) {
    # skip unclassified
    if (cur.class == "None") {
      next
    }
    cur.df <- pol2ratios[pol2ratios$ratioclass == cur.class, ]
    if (nrow(cur.df) > 30) {
      write.table(x = unique(cur.df$GeneID), 
                  file = paste0(outdir2, "/", curname, "_",
                                cur.class, "_genes.txt"),
                  row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
}

# output summary by class
outfile <- paste0(project_dir, "/", curname, "_gene_report.txt")

cat(
  paste0("Initial considered genes: ", nrow(prefiltered.tb), "\n"),
  paste0("Total depleted genes: ", length(depleted_genes), "\n"),
  paste0("Total considered genes: ", nrow(pol2ratios), "\n"),
  paste0("Genes experiencing substantial change: ",
          sum(pol2ratios$ratioclass != "None"), "\n"),
  paste0("Stalled Genes: ", sum(pol2ratios$ratioclass == "Stall"), "\n"),
  paste0("Entering Genes: ",    sum(pol2ratios$ratioclass == "Entry"), "\n"),
  paste0("Pausing Genes: ",  sum(pol2ratios$ratioclass == "Pause"), "\n"),
  paste0("Releasing Genes: ",sum(pol2ratios$ratioclass == "Release"), "\n"),
  paste0("Clogging Genes: ", sum(pol2ratios$ratioclass == "Clog"), "\n"),
  paste0("Unloading Genes: ",sum(pol2ratios$ratioclass == "Unload"), "\n"),
  paste0("Overall Gain Genes: ", length(total_gain_genes), "\n"),
  paste0("Overall Loss Genes: ", length(total_loss_genes), "\n"),
  file = outfile
)



# Cumulative distribution Plots
str.cdf.plot <- makeRatioCDPlot(pol2ratios, "StR", "Stalling Ratio",
                            ctrl_name, case_name, 
                            plot.title=paste0(curname, " Stall Ratio Comparison"))
ggsave(filename=paste0(project_dir,"/",curname,"_stalling_cdf.pdf"),
                        str.cdf.plot, height = 5, width = 7)

par.cdf.plot <- makeRatioCDPlot(pol2ratios, "PaR", "Pausing Ratio",
                            ctrl_name, case_name, 
                            plot.title=paste0(curname, " Pause Ratio Comparison"))
ggsave(filename=paste0(project_dir,"/",curname,"_pausing_cdf.pdf"),
                        par.cdf.plot, height = 5, width = 7)

unr.cdf.plot <- makeRatioCDPlot(pol2ratios, "UnR", "Unloading Ratio",
                            ctrl_name, case_name, 
                            plot.title=paste0(curname, " Unloading Ratio Comparison"))
ggsave(filename=paste0(project_dir,"/",curname,"_unloading_cdf.pdf"),
                        unr.cdf.plot, height = 5, width = 7)

### enrichment plots
dotplot_name <- paste0(project_dir, "/",curname, "_enrichment.pdf")
dotdata_name <- paste0(project_dir, "/",curname, "_enrichment_table.txt")
dot.data <- makeEnrichmentDotPlot(pol2ratios, gene_cat_table, main.classes)
dotplot <- dot.data[[1]]
ggsave(filename=dotplot_name, dotplot, height = 6, width = 6)

write.table(x=as.data.frame(dot.data[[2]]), file=dotdata_name,
            sep="\t", row.names = FALSE,
            col.names = TRUE, quote = FALSE)
cat(" - Enrichment plots made.\n")

### Radar plot
pdf(file=paste0(project_dir, "/",
               curname, "_radar_plot_100perc.pdf"),
               width=6,height=6)
print(getRadar(pol2ratios, 100))
graphics.off()