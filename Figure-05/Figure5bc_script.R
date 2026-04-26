#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(forcats)
  library(stringr)
  library(scales)
})

base_dir <- "/scratch/p312334/project/10--HGT_isolates/4--Function_IBD_association/COG_function_secondary_cluster/robust_category_focus"
out_pdf <- file.path(base_dir, "BCH_cog_direction_consistency_fdr05_nspeciesge2_dom80_combined.pdf")

cat_meta <- tibble::tribble(
  ~category, ~category_name, ~cat_color,
  "B", "Chromatin structure and dynamics", "#4575B4",
  "C", "Energy production and conversion", "#D73027",
  "H", "Coenzyme transport and metabolism", "#1A9850"
)

df <- bind_rows(lapply(cat_meta$category, function(cat) {
  p <- file.path(base_dir, cat, "cog_direction_consistency_fdr05_nspeciesge2_dom80.csv")
  if (!file.exists(p)) {
    return(NULL)
  }
  read_csv(p, show_col_types = FALSE) %>%
    mutate(category = cat)
})) %>%
  left_join(cat_meta, by = "category") %>%
  mutate(
    short_name = str_replace(COG24_FUNCTION_NAME, " \\(PDB:.*$", ""),
    short_name = str_replace(short_name, " \\(PUBMED:.*$", ""),
    label = paste0(COG24_FUNCTION_ID, "\n", str_trunc(short_name, 70)),
    dominant_direction = factor(dominant_direction, levels = c("nonIBD_higher", "IBD_higher")),
    dominant_label = percent(dominant_prop, accuracy = 1),
    neglog10_q = -log10(pmax(min_q, 1e-300)),
    # cap extreme significance so a single point does not distort the figure
    neglog10_q_cap = pmin(neglog10_q, 40),
    facet_title = paste0(category, ": ", category_name)
  ) %>%
  group_by(category) %>%
  mutate(label = fct_reorder(label, n_species_pairs)) %>%
  ungroup()

fill_values <- setNames(cat_meta$cat_color, cat_meta$category)

p <- ggplot(df, aes(x = n_species_pairs, y = label)) +
  geom_segment(aes(x = 0, xend = n_species_pairs, yend = label), linewidth = 0.7, color = "#bdbdbd") +
  geom_point(aes(size = neglog10_q_cap, fill = category), shape = 21, color = "black", stroke = 0.35, alpha = 0.95) +
  geom_text(aes(label = dominant_label), hjust = -0.15, size = 3.2, color = "#444444") +
  facet_grid(facet_title ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(values = fill_values, guide = "none") +
  scale_size_continuous(
    name = expression("Capped " * -log[10] * "(min q)"),
    range = c(3, 12)
  ) +
  coord_cartesian(xlim = c(0, max(df$n_species_pairs, na.rm = TRUE) + 1.2)) +
  theme_classic(base_size = 13) +
  theme(
    axis.title.y = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 0, face = "bold"),
    legend.position = "right",
    plot.title = element_text(face = "bold")
  ) +
  labs(
    title = "Consistent Enriched COGs in Categories B, C, and H",
    subtitle = "Filtered to FDR < 0.05, at least 2 species pairs, and >80% direction consistency",
    x = "Number of significant species pairs"
  )

ggsave(out_pdf, p, width = 11, height = 12)

