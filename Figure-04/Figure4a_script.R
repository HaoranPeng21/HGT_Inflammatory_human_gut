library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  idx <- match(flag, args)
  if (!is.na(idx) && idx < length(args)) args[[idx + 1]] else default
}

base_dir <- "/scratch/p312334/project/10--HGT_isolates/2--Phenotype_association/comparision_plot"
input_path <- get_arg("--input", file.path(base_dir, "IBD_nonIBD_shared_species_pair_HGT_plasmidHGT_nonplasmidHGT_rates_gp10_comparison.csv"))
meta_path <- "/scratch/p312334/project/10--HGT_isolates/data/__binstoget_isolates_summary_v4.1_skani999.csv"
output_pdf <- get_arg("--output_pdf", file.path(base_dir, "IBD_nonIBD_phylogeny_nonplasmidHGT_plasmidHGT_comparison_top5.pdf"))
output_png <- get_arg("--output_png", file.path(base_dir, "IBD_nonIBD_phylogeny_nonplasmidHGT_plasmidHGT_comparison_top5.png"))
subtitle_text <- get_arg("--subtitle", "Non-plasmid HGT and Plasmid HGT only; Top species pairs labeled")
top_n <- as.integer(get_arg("--top_n", "5"))
top_n_low_phy <- as.integer(get_arg("--top_n_low_phy", as.character(top_n)))
top_n_high_phy <- as.integer(get_arg("--top_n_high_phy", as.character(top_n)))
phy_threshold <- as.numeric(get_arg("--phy_threshold", "1"))
ncol_layout <- as.integer(get_arg("--ncol", "2"))

species_map <- read_csv(meta_path, show_col_types = FALSE) %>%
  select(secondary_cluster, Jay_Species) %>%
  filter(!is.na(secondary_cluster), secondary_cluster != "", !is.na(Jay_Species), Jay_Species != "") %>%
  distinct() %>%
  count(secondary_cluster, Jay_Species, sort = TRUE) %>%
  group_by(secondary_cluster) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  transmute(secondary_cluster, species_name = gsub("_", " ", Jay_Species))

plot_df <- read_csv(input_path, show_col_types = FALSE) %>%
  select(
    species_pair,
    phylogeny,
    nonplasmid_HGT_rate_IBD,
    nonplasmid_HGT_rate_non_IBD,
    plasmid_HGT_rate_IBD,
    plasmid_HGT_rate_non_IBD
  ) %>%
  pivot_longer(
    cols = -c(species_pair, phylogeny),
    names_to = c("metric", "cohort"),
    names_pattern = "^(.*)_rate_(IBD|non_IBD)$",
    values_to = "rate"
  ) %>%
  mutate(
    cohort = factor(cohort, levels = c("non_IBD", "IBD"), labels = c("non-IBD", "IBD")),
    metric = factor(
      metric,
      levels = c("nonplasmid_HGT", "plasmid_HGT"),
      labels = c("Non-plasmid HGT", "Plasmid HGT")
    )
  )

line_df <- plot_df %>%
  select(species_pair, phylogeny, metric, cohort, rate) %>%
  pivot_wider(names_from = cohort, values_from = rate) %>%
  transmute(
    species_pair,
    phylogeny,
    metric,
    x = phylogeny,
    xend = phylogeny,
    y = `non-IBD`,
    yend = IBD,
    rate_change = IBD - `non-IBD`
  )

label_df <- bind_rows(
  line_df %>%
    filter(phylogeny < phy_threshold) %>%
    group_by(metric) %>%
    arrange(desc(rate_change), .by_group = TRUE) %>%
    slice_head(n = top_n_low_phy) %>%
    ungroup(),
  line_df %>%
    filter(phylogeny > phy_threshold) %>%
    group_by(metric) %>%
    arrange(desc(rate_change), .by_group = TRUE) %>%
    slice_head(n = top_n_high_phy) %>%
    ungroup()
) %>%
  mutate(
    cluster_1 = sub("-.*$", "", species_pair),
    cluster_2 = sub("^.*-", "", species_pair)
  ) %>%
  left_join(species_map, by = c("cluster_1" = "secondary_cluster")) %>%
  rename(species_1 = species_name) %>%
  left_join(species_map, by = c("cluster_2" = "secondary_cluster")) %>%
  rename(species_2 = species_name) %>%
  mutate(label = paste(coalesce(species_1, cluster_1), coalesce(species_2, cluster_2), sep = " | "))

line_limit <- max(abs(line_df$rate_change), na.rm = TRUE)
if (!is.finite(line_limit) || line_limit == 0) {
  line_limit <- 1
}

p <- ggplot(plot_df, aes(x = phylogeny, y = rate)) +
  geom_segment(
    data = line_df,
    aes(x = x, xend = xend, y = y, yend = yend, color = rate_change),
    inherit.aes = FALSE,
    linewidth = 0.6,
    alpha = 0.9
  ) +
  geom_point(data = filter(plot_df, cohort == "non-IBD"), color = "#4DBBD5", size = 1.8, alpha = 0.9) +
  geom_point(data = filter(plot_df, cohort == "IBD"), color = "#E64B35", size = 1.8, alpha = 0.9) +
  geom_text(
    data = label_df,
    aes(x = phylogeny, y = yend, label = label),
    inherit.aes = FALSE,
    nudge_x = 0.03,
    nudge_y = 0.015,
    size = 3,
    color = "black",
    check_overlap = TRUE
  ) +
  facet_wrap(~ metric, ncol = ncol_layout) +
  scale_colour_gradient2(
    low = "#4DBBD5",
    mid = "#BFBFBF",
    high = "#E64B35",
    midpoint = 0,
    limits = c(-line_limit, line_limit),
    guide = "none"
  ) +
  theme_classic(base_size = 13) +
  labs(
    title = "HGT rates by phylogeny distance",
    subtitle = subtitle_text,
    x = "Phylogeny distance",
    y = "Genome-pair HGT rate"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )

plot_width <- ifelse(ncol_layout == 1, 9, 11)
plot_height <- ifelse(ncol_layout == 1, 9, 5.5)
ggsave(output_pdf, p, width = plot_width, height = plot_height)
ggsave(output_png, p, width = plot_width, height = plot_height, dpi = 300)

message("Saved plot: ", output_pdf)
message("Saved plot: ", output_png)
