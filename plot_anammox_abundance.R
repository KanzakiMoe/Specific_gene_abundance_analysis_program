#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(stringr)

# read CoverM output
df <- read_tsv("results/coverm_all.tsv", show_col_types = FALSE)

# recognize column names
mean_cols <- grep(" Mean$", colnames(df), value = TRUE)
cov_cols  <- grep(" Covered Fraction$", colnames(df), value = TRUE)

if (length(mean_cols) != length(cov_cols)) {
  stop("Mean and Covered Fraction columns do not match")
}

# resolve sample names
samples_mean <- str_remove(mean_cols, " Mean$")
samples_cov  <- str_remove(cov_cols, " Covered Fraction$")

if (!all(samples_mean == samples_cov)) {
  stop("Sample names in Mean and Covered Fraction columns do not match")
}

samples <- samples_mean

# change to long format
mean_long <- df %>%
  select(Contig, all_of(mean_cols)) %>%
  pivot_longer(
    cols = -Contig,
    names_to = "sample",
    values_to = "mean_coverage"
  ) %>%
  mutate(sample = str_remove(sample, " Mean$"))

cov_long <- df %>%
  select(Contig, all_of(cov_cols)) %>%
  pivot_longer(
    cols = -Contig,
    names_to = "sample",
    values_to = "covered_fraction"
  ) %>%
  mutate(sample = str_remove(sample, " Covered Fraction$"))

df_long <- left_join(mean_long, cov_long,
                     by = c("Contig", "sample"))

# summarise total mean coverage and mean covered fraction per sample
summary_df <- df_long %>%
  group_by(sample) %>%
  summarise(
    total_mean_coverage = sum(mean_coverage, na.rm = TRUE),
    mean_covered_fraction = mean(covered_fraction, na.rm = TRUE),
    .groups = "drop"
  )

# export table
write_tsv(summary_df, "results/coverm_summary.tsv")

# plot barplot
p <- ggplot(summary_df,
            aes(x = sample, y = total_mean_coverage)) +
  geom_col(width = 0.6) +
  theme_bw() +
  labs(
    title = "ANAMMOX contig abundance",
    x = "Sample",
    y = "Total Mean Coverage"
  )

ggsave("results/coverm_abundance_barplot.png",
       p, width = 6, height = 4)

# ===== Pie chart based on BAM mapping (Option B) =====
bam_df <- read_tsv("results/bam_mapped_summary.tsv", show_col_types = FALSE)

pie_df <- bam_df %>%
  transmute(
    sample,
    Anammox = mapped,
    `Non-anammox` = total - mapped
  ) %>%
  pivot_longer(
    cols = c(Anammox, `Non-anammox`),
    names_to = "type",
    values_to = "reads"
  ) %>%
  group_by(sample) %>%
  mutate(
    fraction = reads / sum(reads),
    label = sprintf("%.1f%%", fraction * 100)
  )

p_pie <- ggplot(pie_df,
                aes(x = "", y = fraction, fill = type)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 4
  ) +
  facet_wrap(~ sample) +
  theme_void() +
  labs(title = "Fraction of ANAMMOX-mapped reads")

ggsave(
  "results/coverm_abundance_pie.png",
  p_pie,
  width = 6,
  height = 4
)

message("All plots generated successfully")
