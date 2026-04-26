#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggsci)
  library(stringr)
})

pair_path <- "/scratch/p312334/project/10--HGT_isolates/data/all_within_individual_cross_species_genome_pairs_2000bp_hgt_flag.csv"
trait_path <- "/scratch/p312334/project/10--HGT_isolates/data/species_pair_traits_nonzero_gt10.tsv"
dist_path <- "/scratch/p312334/project/10--HGT_isolates/data/Species_pair_distances_ribosomal42_tree.csv"
trait_list_path <- "/scratch/p312334/project/10--HGT_isolates/2--Phenotype_association/results/trait_batch_analysis/trait_summary_both_cohorts_main_p_lt_0.05_downsampling_prop_ge_0.95_direction.csv"
trait_p_path <- "/scratch/p312334/project/10--HGT_isolates/2--Phenotype_association/results/trait_batch_analysis/trait_summary_both_cohorts_main_p_lt_0.05.csv"
out_dir <- "/scratch/p312334/project/10--HGT_isolates/2--Phenotype_association/results/trait_batch_analysis/Plot1"
bin_count <- 30
loess_span <- 1.2

target_trait_ids <- c(
  "gram_positive",
  "growth_cellobiose",
  "enzyme_activity_alkaline_phosphatase_ec3_1_3_1",
  "presence_of_motility",
  "obligate_anaerobic"
)

trait_palettes <- list(
  "gram_positive" = setNames(pal_npg("nrc")(3), c("0", "1", "2")),
  "growth_cellobiose" = setNames(pal_jama()(3), c("0", "1", "2")),
  "enzyme_activity_alkaline_phosphatase_ec3_1_3_1" = setNames(pal_jama()(6)[4:6], c("0", "1", "2")),
  "presence_of_motility" = setNames(pal_simpsons()(3), c("0", "1", "2")),
  "obligate_anaerobic" = setNames(pal_futurama()(3), c("0", "1", "2"))
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pair_df <- read_csv(pair_path, show_col_types = FALSE) %>%
  transmute(
    species_pair = as.character(species_pair),
    HGT = as.integer(HGT)
  ) %>%
  filter(!is.na(species_pair), !is.na(HGT)) %>%
  group_by(species_pair) %>%
  summarise(
    n_total_pairs = n(),
    n_hgt_pairs = sum(HGT == 1L, na.rm = TRUE),
    hgt_rate = n_hgt_pairs / n_total_pairs,
    .groups = "drop"
  )

dist_df <- read_csv(dist_path, show_col_types = FALSE) %>%
  transmute(
    species_pair = as.character(Species_Pair),
    phylo_dist = as.numeric(Average_Distance)
  )

trait_df <- read_tsv(trait_path, show_col_types = FALSE) %>%
  transmute(
    species_pair = as.character(Species_Pair),
    !!!rlang::syms(setdiff(names(.), c("Group1", "Group2", "Species_Pair", "Average_Distance", "species_1", "species_2", "tax_name_1", "tax_name_2")))
  )

trait_list <- read_csv(trait_list_path, show_col_types = FALSE) %>%
  filter(trait_id %in% target_trait_ids)

trait_p <- read_csv(trait_p_path, show_col_types = FALSE) %>%
  transmute(
    trait_id,
    US_main_p = as.numeric(US_main_p),
    NL_main_p = as.numeric(NL_main_p)
  )

trait_list <- trait_list %>%
  left_join(trait_p, by = "trait_id")

for (i in seq_len(nrow(trait_list))) {
  trait_label <- trait_list$trait[[i]]
  trait_id <- trait_list$trait_id[[i]]
  if (!trait_label %in% names(trait_df)) next

  palette <- trait_palettes[[trait_id]]
  if (is.null(palette)) next

  df <- pair_df %>%
    left_join(dist_df, by = "species_pair") %>%
    left_join(trait_df[, c("species_pair", trait_label)], by = "species_pair") %>%
    rename(trait_value = !!trait_label) %>%
    mutate(trait_value = as.character(trait_value)) %>%
    filter(trait_value %in% c("0", "1", "2")) %>%
    mutate(trait_value = factor(trait_value, levels = c("0", "1", "2")))

  if (nrow(df) < 30) next

  df <- df %>% filter(phylo_dist <= 2)
  if (nrow(df) < 10) next

  rng <- range(df$phylo_dist, na.rm = TRUE)
  if (!all(is.finite(rng)) || rng[1] == rng[2]) next

  bins <- seq(rng[1], rng[2], length.out = bin_count + 1)
  df_binned <- df %>%
    mutate(bin = cut(phylo_dist, breaks = bins, include.lowest = TRUE)) %>%
    group_by(trait_value, bin) %>%
    summarise(
      phylo_dist_mean = mean(phylo_dist, na.rm = TRUE),
      hgt_rate_mean = mean(hgt_rate, na.rm = TRUE),
      n_pairs = n(),
      .groups = "drop"
    ) %>%
    group_by(trait_value) %>%
    mutate(n_bins = n()) %>%
    ungroup() %>%
    filter(n_bins >= 6)

  if (nrow(df_binned) < 5) next

  p_us <- trait_list$US_main_p[[i]]
  p_nl <- trait_list$NL_main_p[[i]]
  p_label <- sprintf(
    "US LRT p = %s\nNL LRT p = %s",
    ifelse(is.finite(p_us), format(p_us, digits = 3, scientific = TRUE), "NA"),
    ifelse(is.finite(p_nl), format(p_nl, digits = 3, scientific = TRUE), "NA")
  )

  p <- ggplot(df_binned, aes(x = phylo_dist_mean, y = hgt_rate_mean, color = trait_value, fill = trait_value)) +
    geom_smooth(method = "loess", se = TRUE, linewidth = 1.1, span = loess_span, alpha = 0.18) +
    scale_color_manual(values = palette, drop = FALSE) +
    scale_fill_manual(values = palette, drop = FALSE) +
    coord_cartesian(ylim = c(0, NA)) +
    labs(
      title = trait_label,
      x = "Phylogenetic distance",
      y = "Mean HGT rate (binned)",
      color = "Trait (pair)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      legend.key = element_rect(fill = "white", color = NA)
    ) +
    annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = p_label,
      hjust = -0.02,
      vjust = 1.1,
      size = 3.2,
      color = "black"
    )

  out_pdf <- file.path(out_dir, paste0(str_replace_all(trait_label, "[^A-Za-z0-9]+", "_"), ".pdf"))
  ggsave(out_pdf, p, width = 6.5, height = 4.8, device = grDevices::pdf)
}
