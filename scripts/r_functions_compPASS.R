# libraries
# process Pol2 compare table
options(warn = -1)

myPaths <- .libPaths()
conda_R_path <- myPaths[grepl("/envs/compPASS/", myPaths)]

suppressPackageStartupMessages(library(pheatmap,
                               lib.loc = conda_R_path))
suppressPackageStartupMessages(library(dplyr,
                               lib.loc = conda_R_path))
suppressPackageStartupMessages(library(tidyr,
                               lib.loc = conda_R_path))
suppressPackageStartupMessages(library(cowplot,
                               lib.loc = conda_R_path))
suppressPackageStartupMessages(library(ggplot2,
                               lib.loc = conda_R_path))
suppressPackageStartupMessages(library(ggrepel,
                               lib.loc = conda_R_path))
suppressPackageStartupMessages(library(UpSetR,
                               lib.loc = conda_R_path))

if (!require(fmsb)) {
  install.packages("fmsb", repo = "cran.case.edu")
  suppressPackageStartupMessages(library(fmsb,
                                 version="0.7.6"))
}


options(scipen = 999)

# Functions for calculating and plotting R data
CHR_EXCLUSION_LIST <- c("chrY")

REGION_COUNT <- 4

.GlobalEnv$region_pro_start <- -800
.GlobalEnv$region_tssr_start <- -30
.GlobalEnv$region_gene_start <- 300
.GlobalEnv$region_tesr_end <- 4000

# Fold difference expected of starting
RATIO_THRESHOLD <- 1.999

# writes out a tab-delim table with unix line-endings
writeTable_unix <- function(x, y, cols = FALSE) {
  output.file <- file(y, "wb")
  write.table(file = output.file,
              x, quote = FALSE, row.names = FALSE,
              col.names = cols, sep = "\t")
  close(output.file)
}

# functions for generate_comparison_table_from_multicov.R
calculateVals <- function(vct.a, vct.b, gene_len) {
  
  #scaled_vct.a <- c(vct.a[1] / ((.GlobalEnv$region_tssr_start - .GlobalEnv$region_pro_start)+1),
  #        vct.a[2] / ((.GlobalEnv$region_gene_start - .GlobalEnv$region_tssr_start)+1),
  #        vct.a[3] / (gene_len+1),
  #        vct.a[4] / ((.GlobalEnv$region_tesr_end)+1))
          
  #scaled_vct.b <- c(vct.b[1] / ((.GlobalEnv$region_tssr_start - .GlobalEnv$region_pro_start)+1),
  #        vct.b[2] / ((.GlobalEnv$region_gene_start - .GlobalEnv$region_tssr_start)+1),
  #        vct.b[3] / (gene_len+1),
  #        vct.b[4] / ((.GlobalEnv$region_tesr_end)+1))
          
  #scaled_vct.a <- round(scaled_vct.a,6)
  #scaled_vct.b <- round(scaled_vct.b,6)
  # get sample_a solo values
  out.a <- c(vct.a, get_stall_ratio(vct.a),
             get_pause_ratio(vct.a, gene_len),
             get_unloading_ratio(vct.a, gene_len))[1:7]

  # get sample_b solo values
  out.b <- c(vct.b, get_stall_ratio(vct.b),
             get_pause_ratio(vct.b, gene_len),
             get_unloading_ratio(vct.b, gene_len))[1:7]
  # get comparison l2fc values
  j <- 1
  comp.vector <- c()
  while(j <= length(out.a)) {
      cur.l2fc <- round(log2(as.numeric(out.b[j]) + 0.1) -
                        log2(as.numeric(out.a[j]) + 0.1), 3)
      comp.vector <- append(comp.vector, cur.l2fc)
      j <- j + 1
  }
  return(c(out.a, out.b, comp.vector))
}

get_stall_ratio <- function(x) {
  ret_val <- 0
  prom_len <- .GlobalEnv$region_tssr_start - .GlobalEnv$region_pro_start + 1
  tssr_len <- .GlobalEnv$region_gene_start - .GlobalEnv$region_tssr_start + 1
  x <- as.numeric(x)
  if (x[2] != 0) {
    ret_val <- ((x[1] + 0.01) / prom_len) / ((x[2] + 0.01) / tssr_len)
  }
  return(round(ret_val, 3))
}

get_pause_ratio <- function(x, glen) {
  ret_val <- 0
  tssr_len <- .GlobalEnv$region_gene_start - .GlobalEnv$region_tssr_start + 1
  x <- as.numeric(x)
  if (x[3] != 0) {
    ret_val <- ((x[2] + 0.01) / tssr_len) / ((x[3] + 0.01) / (glen + 1))
  }
  return(round(ret_val, 3))
}

get_unloading_ratio <- function(x, glen) {
  ret_val <- 0
  tesr_len <- .GlobalEnv$region_tesr_end + 1
  x <- as.numeric(x)
  if (x[4] != 0) {
    ret_val <- ((x[3] + 0.01)/(glen + 1)) / ((x[4] + 0.01) / tesr_len)
  }
  return(round(ret_val, 3))
}
# note change

# functions for classify_genes_from_pol2_states.R
isHiOutlier <- function(x) {
  # original program
  cur.iqr <- IQR(x)
  low.bd <- summary(x)[2]
  hi.bd <- summary(x)[5]

  return(x >= (1.5 * cur.iqr + hi.bd))
}

isLoOutlier <- function(x) {
  # original program
  cur.iqr <- IQR(x)
  low.bd <- summary(x)[2]
  hi.bd <- summary(x)[5]

  return(x <= (low.bd - (1.5 * cur.iqr)))
}


getTukeyRangeNonNeg <- function(n.vect){
  iqr <- IQR(n.vect)
  low.val <- summary(n.vect)[2] - (2 * iqr)
  if(low.val < 0){
    nonz <- n.vect[n.vect != 0]
    low.val = min(nonz)
    }
  
  hi.val <- summary(n.vect)[5] + (2 * iqr)
  return(c(low.val, hi.val))
  
}

getTukeyRangePick <- function(x1, x2, isNonNeg = FALSE){
  rng1 <- range(x1)[2] - range(x1)[1]
  rng2 <- range(x2)[2] - range(x2)[1]
  
  if(rng2 > rng1){
    return(getTukeyRange(x2, isNonNeg))
  } else {
    return(getTukeyRange(x1, isNonNeg))
  }
}

getTukeyRange <- function(n.vect, isNonNeg = FALSE){
  if (isNonNeg){
    return(getTukeyRangeNonNeg(n.vect))
  }
  iqr <- IQR(n.vect)
  low.val <- summary(n.vect)[2] - (2 * iqr)
  
  hi.val <- summary(n.vect)[5] + (2 * iqr)
  return(c(low.val, hi.val))
}

# generate Pol2 plots for each ROI for each version
generateHistogramOneSample <- function(input.df, target.col, plottitle,
                               x.label, y.label, bincount, mincutoff = 1,
                               determine_scale_x = TRUE, scale_x = c(0, 1),
                               select.fill = "cyan", select.col = "black") {
  input.df2 <- input.df[input.df[, target.col] > mincutoff, ]
  if (determine_scale_x){
    scale_x <- getTukeyRangeNonNeg((input.df2[, target.col]))
  }

  p <- ggplot(input.df2) +
    geom_histogram(aes(x = log2(.data[[target.col]])), bins = bincount,
                   fill = select.fill, color = select.col) +
    ggtitle(plottitle) +
    scale_x_continuous(x.label, limits = log2(scale_x)) +
    theme_light()
  return(p)
}

# for each region in case and control, plot RPM density
makeDensityDistributionPlots <- function(pol2ratios) {
  regions.all <- c("Prom", "TSSR", "GB", "TESR")
  region.names.all <- c("Promoter", "TSSR", "Gene Body", "TESR")

  i <- 1
  plots.out <- list()
  while(i <= length(regions.all)) {
    plots.out[[i]] <- makeRPMDistributionPlot(
      pol2ratios,
      region.names.all[i], regions.all[i])
    i <- i + 1
  }

  out.plot <- plot_grid(plotlist = plots.out, ncol=REGION_COUNT, nrow = 1)
  return(out.plot)
}

makeRatioPlot <- function(pol2ratios, class1, class2,
                          sizer, x.val, y.val, x.scale,
                          y.scale, red.group, blue.group,
                          curname) {
    # set size and shape based on
  pol2ratios$shapeClass <- "middle"
  pol2ratios$shapeClass[pol2ratios[[paste0(sizer,"2")]] >= RATIO_THRESHOLD] <- "highEnd" # triangle
  pol2ratios$shapeClass[pol2ratios[[paste0(sizer,"1")]] >= RATIO_THRESHOLD] <- "highStart" # square
  
  red.df <- subset(pol2ratios, pol2ratios[, class1] != "")
  blue.df <- subset(pol2ratios, pol2ratios[, class2] != "")

  x.max.val <- round(max(pol2ratios[[x.val]]) / 0.5, 0) * 0.5
  x.min.val <- round(min(pol2ratios[[x.val]]) / 0.5, 0) * 0.5
 
  y.max.val <- round(max(pol2ratios[[y.val]]) / 0.5, 0) * 0.5
  y.min.val <- round(min(pol2ratios[[y.val]]) / 0.5, 0) * 0.5

  #x.max.val <- x.vals[1]
  #x.min.val <- x.vals[2]
  #
  #y.max.val <- y.vals[1]
  #y.min.val <- y.vals[2]
  named.vector.shapes <- c("highStart" = 22, "middle" = 1, "highEnd" = 25)

  out <- ggplot(pol2ratios, aes(x = .data[[x.val]],
        #y=.data[[y.val]], shape=.data[["shapeClass"]]))+
        y = .data[[y.val]])) +
    geom_point(alpha = 0.1) +
    geom_point(data = red.df, color="red") +
    geom_point(data = blue.df, color="blue") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_y_continuous(y.scale, limits=c(y.min.val, y.max.val)) +
    scale_x_continuous(x.scale, limits=c(x.min.val, x.max.val)) +
    scale_shape_manual(values = named.vector.shapes) +
    ggtitle(paste0(red.group, " (red) vs ", blue.group, " (blue)\n",
            red.group, " (n = ",dim(red.df)[1], "), ",
            blue.group, " , (n = ",dim(blue.df)[1], ")")) +
    guides(shape = "none") +
    theme_bw()

  return(out)
}

# for a given region, make RPM density
makeRPMDistributionPlot <- function(pol2ratios, region.name, region.short) {
  target.value1 <- paste0(region.short, "1")
  target.value2 <- paste0(region.short, "2")
  
  ctrl.limits <- getTukeyRange(pol2ratios[[target.value1]], isNonNeg = TRUE)
  case.limits <- getTukeyRange(pol2ratios[[target.value2]], isNonNeg = TRUE)

  overall.max <- max(c(ctrl.limits,case.limits))
  limits <- c(0, overall.max)
  
  ctrl.plot <- ggplot(pol2ratios)+
    geom_histogram(aes(x = .data[[target.value1]]), bins = 500, color = "gray") +
    scale_x_continuous(paste0("Distribution of ", region.name, " RPM"),
                       limits = limits) +
    scale_y_continuous("Count") +
    theme_light()
  
  case.plot <- ggplot(pol2ratios)+
    geom_histogram(aes(x=.data[[target.value2]]), bins = 500, color = "cyan") +
    scale_x_continuous(paste0("Distribution of ", region.name, " RPM"),
                       limits = limits) +
    scale_y_continuous("Count") +
    theme_light()

   # mix.plot <- ggplot(pol2ratios)+
   # geom_histogram(aes(x=.data[[target.value1]]), bins=500, color="gray", alpha=0.9)+
   # geom_histogram(aes(x=.data[[target.value2]]), bins=500, color="cyan", alpha=0.3)+
   # scale_x_continuous(paste0("Distribution of ", region.name, " RPM Density"),
   #                    limits=limits)+
   # scale_y_continuous("Count")+
   # ggtitle("Control (gray), Case (blue)")+
   # theme_light()
  
  out.plot <- plot_grid(plotlist = list(ctrl.plot, case.plot), ncol = 1, nrow = 2)
  return(out.plot)
}

# make plot of mixed categories from ratios, for all classes > min.count
makeUpsetPlot <- function(pol2ratios, min.count) {
  total.intersects <- unique(pol2ratios$mergedclass)
  total.intersects <- total.intersects[total.intersects != ""]
  main.classes <- c("Stall", "Entry",
                    "Pause", "Release",
                    "Clog", "Unload",
                    "Gain", "Loss")
  
  ct.tb <- table(pol2ratios$mergedclass)
  ct.tb.in <- ct.tb[ct.tb >= min.count]
  ct.tb.in <- ct.tb.in[names(ct.tb.in) != ""]
  mix.classes <- names(ct.tb.in)
  mix.classes <- mix.classes[!mix.classes %in% c("", main.classes)]
  
  total.intersects.n <- length(mix.classes)
  set.order <- c(main.classes, mix.classes)
  
  # upset plot
  upset.tb <- list(Stall=which(pol2ratios$stall == "Stall"),
                   Pause=which(pol2ratios$pause == "Pause"),
                   Clog=which(pol2ratios$clog == "Clog"),
                   Entry=which(pol2ratios$entry == "Entry"),
                   Release=which(pol2ratios$release == "Release"),
                   Unload=which(pol2ratios$Unload == "Unload"),
                   Gain=which(premergedclass == "" & pol2ratios$GeneID %in% total_gain_genes ),
                   Loss=which(premergedclass == "" & pol2ratios$GeneID %in% total_loss_genes ))
  upset.plot <- upset(fromList(upset.tb), sets=rev(main.classes), keep.order=T, group.by = "degree",
                      text.scale=1.5)
  return(upset.plot)
}

# cumulative distribution plot (will be appended with 1 or 2 for either sample)
makeRatioCDPlot <- function(pol2ratios, ratio.short, value.name,
                            grp1.name, grp2.name, plot.title = "") {
  vct.nm1 <- paste0(ratio.short, "1")
  vct.nm2 <- paste0(ratio.short, "2")
  
  cdf.fct1 <- ecdf(log2(pol2ratios[[vct.nm1]] + 0.1))
  cdf.fct2 <- ecdf(log2(pol2ratios[[vct.nm2]] + 0.1))
  
  cumdf.vct1 <- cdf.fct1(log2(pol2ratios[[vct.nm1]] + 0.1))
  cumdf.vct2 <- cdf.fct2(log2(pol2ratios[[vct.nm2]] + 0.1))
  
  tb1 <- data.frame(GeneID=pol2ratios$GeneID,
                    group=rep(grp1.name, length(cumdf.vct1)),
                    log_val=log2(pol2ratios[[vct.nm1]] + 0.1),
                    cumdist=cumdf.vct1)
  tb2 <- data.frame(GeneID=pol2ratios$GeneID,
                    group=rep(grp2.name, length(cumdf.vct2)),
                    log_val=log2(pol2ratios[[vct.nm2]] + 0.1),
                    cumdist=cumdf.vct2)

  sig.val <- (wilcox.test(tb1$log_val, tb2$log_val))$p.value
  sig.statement <- ifelse(sig.val < 0.05, "Distributions significantly differ",
                          "No significant difference")
  sig.val.out <- formatC(sig.val, digits = 3)

  plot.tb <- rbind(tb1, tb2)

  min.val <- floor(min(plot.tb$log_val)/2) * 2
  max.val <- ceiling(max(plot.tb$log_val)/2) * 2

  plot.title.out <- ifelse(plot.title == "", 
                           paste0(value.name, " Comparison"),
                           plot.title)

  colors.out <- c("black", "red")
  cplot <- ggplot(plot.tb) +
    geom_line(aes(x=log_val, y=cumdist*100, 
                  color=factor(group, levels = c(grp1.name, grp2.name)))) +
    scale_color_manual("Groups", values = colors.out) +
    scale_x_continuous(paste0("Log2 ", value.name),
                       breaks=seq(min.val, max.val, 2),
                       limits=c(min.val, max.val)) +
    scale_y_continuous("Cumulative Distribution", 
                        limits=c(0,100), breaks = seq(0, 100, 20)) +
    ggtitle(paste0(plot.title.out, "\n",
                   sig.statement, " (p.val = ", sig.val.out,
                   ", n=", dim(pol2ratios)[1]," )")) +
    theme_light()
  
  return(cplot)
}

# make dot plot with gene category enrichment
makeEnrichmentDotPlot <- function(pol2ratios, gene_cat, focus.grps) {
  # left join gene_cat with decided class
  
  # rows from column names - first two are ID and length
  gene.cat.names <- colnames(gene_cat)[3:ncol(gene_cat)]
  
  out.tb <- NULL
  for(cur.cat in gene.cat.names) {
    pvals.list <- c()
    for (grp in focus.grps) {
      # get hypergeometric enrichment score for each group
      select_set <- pol2ratios$GeneID[pol2ratios$GeneID %in% gene_cat$GeneID[gene_cat[[cur.cat]]]]
      pvals.list <- append(pvals.list,
                           getCategoryPvalue(pol2ratios,
                                             select_set,
                                             grp))
    }
    out.tb <- rbind(out.tb, c(cur.cat, pvals.list))
  }
  out.tb <- as.data.frame(out.tb)
  colnames(out.tb) <- c("Gene Class", focus.grps)
  out.tb[,c(2:ncol(out.tb))] <- apply(out.tb[,c(2:ncol(out.tb))], 2, as.numeric)
  
  out.tb.pv <- pivot_longer(out.tb, cols=c(-"Gene Class"),
                            names_to = "outcome",values_to = "pLog_hyper")
  y.axis.order <- c("Stall", "Entry",
                    "Pause", "Release",
                    "Clog", "Unload",
                    "Gain", "Loss")
  #x.axis.order
  out.tb.pv$isSig <- ifelse(out.tb.pv$pLog_hyper > 1.3, TRUE, FALSE)
  out.plot <- ggplot(out.tb.pv)+
    geom_point(aes(x = `Gene Class`, y = factor(outcome, levels=focus.grps),
                   size = pLog_hyper, color = isSig)) +
    scale_y_discrete("") +
    scale_x_discrete("") +
    scale_color_manual(values = c("grey", "orange")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))

  # melt table
  return(list(out.plot, out.tb))
  
  # bubble plot

}

getCategoryPvalue <- function(pol2ratios, select_set, target_var) {
  gene_subset <- pol2ratios[pol2ratios$GeneID %in% select_set, ]

  # Total number of genes
  N <- nrow(pol2ratios)
  
  # Total number of "successes" in population (target_var == TRUE)
  K <- length(select_set)# i.e. genes where category is T
  
  # Size of subset
  n <- sum(pol2ratios$mergedclass == target_var)
  
  # Number of "successes" in subset
  k <- sum(gene_subset$mergedclass == target_var)
  
  p_val <- (phyper(k - 1, K, N - K, n, lower.tail = FALSE))
  return(round(-log10(p_val),3))
}

#compareRatiosToGEX <- function(pol2ratios, tpm.file) {
  #
#}

# creates a radar plot showing percentages from 0% to a set %
## of perturbed genes occupying one of the eight classes of
## Pol2 activity shift
getRadar <- function(results.case, set_max = 100) {
  MAX.VAL <- set_max
  main.classes <- c("Stall", "Entry",
                    "Pause", "Release",
                    "Clog", "Unload",
                    "Gain", "Loss"
                    )
  temp.rtb <- table(results.case$ratioclass[results.case$ratioclass != "None"])
  tot.perturbed <- sum(results.case$ratioclass != "None")
  results.case.tb <- data.frame(RatioClass=names(temp.rtb),
                                Percentage=round(as.numeric(temp.rtb) /
                                                 tot.perturbed * 100,2))#,
  #Count=as.numeric(temp.rtb))
  results.case.tb2 <- pivot_wider(results.case.tb,
                                  names_from = "RatioClass", 
                                  values_from = "Percentage")
  for (cur.class in main.classes) {
    if (!cur.class %in% colnames(results.case.tb2)) {
      results.case.tb2[[cur.class]] <- 0
    }
  }
  
  radar_order <- c("Gain", "Stall",
                   "Pause", "Clog",
                   "Loss", "Entry",
                   "Release", "Unload")
  results.case.tb2 <- results.case.tb2[, radar_order]
  num.cat <- length(main.classes)
  
  results.case.tb_out <- rbind(rep(MAX.VAL, num.cat),
                               rep(0, num.cat),
                               results.case.tb2)
  
  outplot <- radarchart(results.case.tb_out,
             axistype = 0,
             seg = 4, # change as a multiple of max value
             cglty = 1,
             axislabcol = "black",
             pfcol = rgb(0, 0.1, 0.9, 0.2), # red, green, blue ratio, then alpha
             vlcex = 0.9)
  return(outplot)
}