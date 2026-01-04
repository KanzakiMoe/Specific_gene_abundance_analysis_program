#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(stringr)
library(ggplot2)

infile  <- "results/anammox_coverm.tsv"
outfile_table <- "results/anammox_abundance_summary.tsv"
outfile_plot  <- "results/anammox_abundance_barplot.pdf"

df <- read_tsv(infile, show_col_types = FALSE)

# Extract Mean / Covered Fraction columns
mean_cols <- grep(" Mean$", colnames(df), value = TRUE)
cov_cols  <- grep(" Covered Fraction$", colnames(df), value = TRUE)

if (length(mean_cols) == 0 || length(mean_cols) != length(cov_cols)) {
  stop("Mean and Covered Fraction columns do not match")
}

# resolve sample names
samples <- str_replace(mean_cols, " Mean$", "")

# calculate abundance
result <- data.frame(
  sample = samples,
  abundance = sapply(seq_along(samples), function(i) {
    sum(df[[mean_cols[i]]] * df[[cov_cols[i]]])
  })
)

# export table
write_tsv(result, outfile_table)

# plot barplot
p <- ggplot(result, aes(x = sample, y = abundance)) +
  geom_col(width = 0.6) +
  xlab("Sample") +
  ylab("Anammox gene abundance\n(Mean × Covered Fraction)") +
  theme_bw(base_size = 14)

ggsave(outfile_plot, p, width = 5, height = 4)
