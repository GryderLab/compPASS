# Make RNA comparison to Pol2 ratios
# untested

# Load libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# Arguments: two TPM files, one ratio_table file, output PDF
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: script.R <tpm1.tsv> <tpm2.tsv> <ratio_table.tsv> <output.pdf>")
}

tpm1_file <- args[1]
tpm2_file <- args[2]
ratio_table_file <- args[3]
output_pdf <- args[4]

# Read input (TPM and ratio table)
tpm1 <- read.table(tpm1_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
tpm2 <- read.table(tpm2_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
ratio_table <- read.table(ratio_table_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Ensure GeneID column exists
stopifnot("GeneID" %in% colnames(tpm1))
stopifnot("GeneID" %in% colnames(tpm2))
stopifnot("GeneID" %in% colnames(ratio_table))

# Filter TPM files down to only genes in ratio_table
tpm1 <- tpm1[tpm1$GeneID %in% ratio_table$GeneID, ]
tpm2 <- tpm2[tpm2$GeneID %in% ratio_table$GeneID, ]

# Merge TPM values into ratio_table
merged <- ratio_table %>%
  left_join(tpm1 %>% dplyr::select(GeneID, TPM1 = TPM), by = "GeneID") %>%
  left_join(tpm2 %>% dplyr::select(GeneID, TPM2 = TPM), by = "GeneID")

# Reshape L2FC values for plotting
plot_df <- merged %>%
  dplyr::select(GeneID, ratioclass, dplyr::starts_with("L2FC_")) %>%
  pivot_longer(cols = starts_with("L2FC_"),
               names_to = "Comparison",
               values_to = "L2FC")

# Plot
p <- ggplot(plot_df, aes(x = ratioclass, y = L2FC, fill = ratioclass)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~Comparison, scales = "free_y") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Ratio Class", y = "Log2 Fold Change",
       title = "L2FC by Ratio Class")

# Save to PDF
ggsave(output_pdf, p, width = 10, height = 6)
