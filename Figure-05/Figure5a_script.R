#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(forcats)
  library(scales)
})

in_path <- "/scratch/p312334/project/10--HGT_isolates/4--Function_IBD_association/COG_function_secondary_cluster/COG24_IBD_summary_annotated.csv"
out_pdf <- "/scratch/p312334/project/10--HGT_isolates/4--Function_IBD_association/COG_function_secondary_cluster/COG24_IBD_secondary_cluster_downsampling_lollipop_fully100.pdf"
out_csv <- "/scratch/p312334/project/10--HGT_isolates/4--Function_IBD_association/COG_function_secondary_cluster/COG24_IBD_secondary_cluster_downsampling_fully100.csv"

cat_colors <- c(
  "C" = "#D73027",
  "B" = "#4575B4",
  "H" = "#1A9850",
  "R" = "#984EA3",
  "K" = "#FF7F00",
  "U" = "#A65628",
  "L" = "#4DAF4A"
)

df <- read_csv(in_path, show_col_types = FALSE) %>%
  filter(prop_beta_pos %in% c(0, 1), down_prop_sig == 1) %>%
  mutate(
    label = paste0(category_name, " (", letter, ")"),
    down_label = percent(down_prop_sig, accuracy = 1),
    label = fct_reorder(label, median_or)
  )

write_csv(df, out_csv)

p <- ggplot(df, aes(x = median_or, y = label, color = letter)) +
  geom_segment(aes(x = 1, xend = median_or, y = label, yend = label), linewidth = 0.7, alpha = 0.85) +
  geom_point(aes(size = down_prop_sig), alpha = 0.95) +
  geom_text(aes(label = down_label), hjust = -0.2, size = 3.2, color = "#444444") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#888888", linewidth = 0.5) +
  scale_color_manual(
    values = cat_colors
  ) +
  scale_size_continuous(
    name = "Downsampling\nsig. proportion",
    labels = percent_format(accuracy = 1),
    range = c(2.5, 8),
    limits = c(0, 1)
  ) +
  coord_cartesian(xlim = c(min(0.75, min(df$median_or, na.rm = TRUE)), max(df$median_or, na.rm = TRUE) + 1.2)) +
  theme_classic(base_size = 13) +
  theme(
    axis.title.y = element_blank(),
    legend.title = element_text(size = 11),
    legend.position = "right",
    plot.title = element_text(face = "bold")
  ) +
  labs(
    title = "COG24 Categories: Secondary-Cluster Downsampling Robustness",
    subtitle = "Only categories with 100% sign consistency and 100% significant downsampling runs",
    x = "Median odds ratio (IBD / non-IBD) across downsampling runs",
    color = "COG24\ncategory"
  )

ggsave(out_pdf, p, width = 10.5, height = 7.8)
