#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
})

result_dir <- "/scratch/p312334/project/10--HGT_isolates/2--Phenotype_association/results/all_phenotypes_aggregate_model"
input_path <- file.path(result_dir, "stable_phenotypes_or_plot_data_final5.csv")
filtered_path <- file.path(result_dir, "stable_phenotypes_or_plot_data_final5_down100_noIBD.csv")
pdf_path <- file.path(result_dir, "stable_phenotypes_or_lrt_lollipop_final5_R.pdf")
png_path <- file.path(result_dir, "stable_phenotypes_or_lrt_lollipop_final5_R.png")

plot_df <- read_csv(input_path, show_col_types = FALSE) %>%
  filter(
    phenotype != "IBD",
    round(down_prop_LRT_p_lt_0_05, 12) == 1
  ) %>%
  mutate(
    phenotype = factor(
      phenotype,
      levels = rev(c(
        "Anemia",
        "PetsChild0to15Y.Any",
        "HeartRateComplains",
        "Urbanicity"
      ))
    ),
    direction = factor(direction, levels = c("OR < 1", "OR > 1"))
  )

write_csv(plot_df, filtered_path)

p <- ggplot(plot_df, aes(y = phenotype, x = OR_used, colour = direction)) +
  geom_segment(aes(x = 1, xend = OR_used, yend = phenotype), linewidth = 1.4, alpha = 0.8) +
  geom_point(aes(size = neglog10_LRT_p), alpha = 0.95) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey45", linewidth = 0.7) +
  scale_colour_manual(
    values = c("OR < 1" = "#2f76c7", "OR > 1" = "#c83b2d"),
    name = "Direction"
  ) +
  scale_size_continuous(
    name = expression(-log[10]("main LRT P")),
    breaks = c(25, 50, 75),
    range = c(4, 14)
  ) +
  scale_x_continuous(
    limits = c(0.8, 2.85),
    breaks = c(1.0, 1.5, 2.0, 2.5)
  ) +
  labs(
    title = "Stable phenotype associations",
    x = "Odds ratio",
    y = "Phenotype",
    caption = "Final phenotypes pass: main LRT P < 0.05, leave-top5 P < 0.05\nDownsampling P < 0.05 in 100% iterations, direction consistency = 1"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 18, hjust = 0, margin = margin(b = 8)),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 11, colour = "black"),
    axis.line = element_line(linewidth = 0.7, colour = "black"),
    axis.ticks = element_line(linewidth = 0.7, colour = "black"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 11),
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = unit(0.35, "cm"),
    plot.caption = element_text(size = 8.8, hjust = 0, colour = "grey30", lineheight = 1.05),
    plot.margin = margin(t = 8, r = 12, b = 6, l = 8)
  ) +
  guides(
    colour = guide_legend(order = 1, override.aes = list(size = 4, linewidth = 1.4)),
    size = guide_legend(order = 2)
  )

ggsave(pdf_path, p, width = 7.4, height = 4.6, device = cairo_pdf)
ggsave(png_path, p, width = 7.4, height = 4.6, dpi = 300)
