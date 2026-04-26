#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggsci)
})

in_path <- "/scratch/p312334/project/10--HGT_isolates/2--Phenotype_association/results/infalmmation/with_specific_meds/HGT__Lab1Calprotectin/df_sp.csv"
out_dir <- "/scratch/p312334/project/10--HGT_isolates/2--Phenotype_association/results/infalmmation/Plot"
out_pdf <- file.path(out_dir, "Lab1Calprotectin_cut200_vs_phylogeny_loess.pdf")
out_binned_csv <- file.path(out_dir, "Lab1Calprotectin_cut200_vs_phylogeny_binned.csv")
bin_count <- 30
loess_span <- 1.2

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

df <- read_csv(in_path, show_col_types = FALSE) %>%
  transmute(
    phylogeny = as.numeric(phylogeny),
    n_total_pairs = as.numeric(n_total_pairs),
    n_hgt_pairs = as.numeric(n_hgt_pairs),
    Lab1Calprotectin = as.numeric(Phenotype)
  ) %>%
  filter(
    !is.na(phylogeny),
    !is.na(n_total_pairs),
    !is.na(n_hgt_pairs),
    !is.na(Lab1Calprotectin),
    n_total_pairs > 0
  ) %>%
  mutate(
    hgt_rate = n_hgt_pairs / n_total_pairs,
    calpro_group = if_else(Lab1Calprotectin > 200, "> 200", "<= 200")
  ) %>%
  filter(phylogeny <= 2) %>%
  mutate(calpro_group = factor(calpro_group, levels = c("<= 200", "> 200")))

if (nrow(df) < 10) {
  stop("Too few rows after filtering.")
}

rng <- range(df$phylogeny, na.rm = TRUE)
if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) {
  stop("Invalid phylogeny range.")
}

bins <- seq(rng[1], rng[2], length.out = bin_count + 1)
df_binned <- df %>%
  mutate(bin = cut(phylogeny, breaks = bins, include.lowest = TRUE)) %>%
  group_by(calpro_group, bin) %>%
  summarise(
    phylogeny_mean = mean(phylogeny, na.rm = TRUE),
    hgt_rate_mean = mean(hgt_rate, na.rm = TRUE),
    n_rows = n(),
    .groups = "drop"
  ) %>%
  group_by(calpro_group) %>%
  mutate(n_bins = n()) %>%
  ungroup() %>%
  filter(n_bins >= 6)

write_csv(df_binned, out_binned_csv)

palette <- c("<= 200" = "#1f77b4", "> 200" = "#d62728")

n_label <- df %>%
  count(calpro_group) %>%
  mutate(txt = paste0(as.character(calpro_group), ": n=", n)) %>%
  summarise(lbl = paste(txt, collapse = "\n")) %>%
  pull(lbl)

p_label <- paste(
  n_label,
  "LRT p (unadjusted) = 1.452234e-05",
  "LRT p (medication-adjusted) = 1.751148e-04",
  sep = "\n"
)

p <- ggplot(df_binned, aes(x = phylogeny_mean, y = hgt_rate_mean, color = calpro_group, fill = calpro_group)) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1.1, span = loess_span, alpha = 0.18) +
  scale_color_manual(values = palette, drop = FALSE) +
  scale_fill_manual(values = palette, drop = FALSE) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(
    title = "Lab1Calprotectin (>200 vs <=200)",
    x = "Phylogenetic distance",
    y = "Mean HGT rate (binned)",
    color = "Lab1Calprotectin"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  ) +
  annotate("text", x = -Inf, y = Inf, label = p_label, hjust = -0.02, vjust = 1.1, size = 3.2, color = "black")

ggsave(out_pdf, p, width = 6.5, height = 4.8, device = cairo_pdf)
message("Saved: ", out_pdf)
message("Saved: ", out_binned_csv)
