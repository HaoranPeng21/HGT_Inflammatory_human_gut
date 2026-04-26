#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(lme4)
  library(lmtest)
})

pair_path <- "/scratch/p312334/project/10--HGT_isolates/data/all_within_individual_cross_species_genome_pairs_2000bp_hgt_flag.csv"
pheno_path <- "/scratch/p312334/project/10--HGT_isolates/data/full_phen_RISEUP_20260329.csv"
out_dir <- "/scratch/p312334/project/10--HGT_isolates/2--Phenotype_association/results/all_phenotypes_aggregate_model/shared_species_pair_hgt_rate_plots_gp10"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

MIN_GENOME_PAIRS <- 10
ibd_color <- "#E64B35"
healthy_color <- "#4DBBD5"
diff_mid <- "#D9D9D9"

normalize_id <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x)] <- ""
  ifelse(grepl("^p\\d+$", x), sprintf("P%02d", as.integer(sub("^p", "", x))), x)
}

build_pheno_lookup <- function(pheno_df) {
  bind_rows(
    pheno_df %>% mutate(join_id = normalize_id(OldID)),
    pheno_df %>% mutate(join_id = normalize_id(sample_id))
  ) %>%
    filter(join_id != "") %>%
    group_by(join_id) %>%
    slice(1) %>%
    ungroup()
}

format_p <- function(x) {
  if (is.na(x)) return("NA")
  if (x == 0) return("< 1e-300")
  if (x < 1e-4) return(format(x, scientific = TRUE, digits = 3))
  format(x, digits = 4)
}

pair_df <- read_csv(pair_path, show_col_types = FALSE)
pheno_df <- read_csv(pheno_path, show_col_types = FALSE)
pheno_lookup <- build_pheno_lookup(pheno_df)

annotated_df <- pair_df %>%
  mutate(join_id = normalize_id(MGS_ID)) %>%
  left_join(pheno_lookup %>% select(join_id, IBD), by = "join_id") %>%
  mutate(
    IBD = suppressWarnings(as.numeric(IBD)),
    phenotype_group = case_when(
      IBD == 1 ~ "IBD",
      IBD == 0 ~ "non-IBD",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(phenotype_group), !is.na(species_pair), !is.na(phylogeny), !is.na(HGT)) %>%
  mutate(phenotype_group = factor(phenotype_group, levels = c("IBD", "non-IBD")))

agg_df <- annotated_df %>%
  transmute(
    species_pair = as.character(species_pair),
    phylogeny = as.numeric(phylogeny),
    phenotype_group,
    flag = as.integer(HGT)
  ) %>%
  group_by(species_pair, phenotype_group) %>%
  summarise(
    total_genome_pair = n(),
    positive_genome_pair = sum(flag == 1L, na.rm = TRUE),
    transfer_rate = positive_genome_pair / total_genome_pair,
    phylogeny = dplyr::first(phylogeny),
    .groups = "drop"
  )

shared_species_pairs <- agg_df %>%
  filter(total_genome_pair > MIN_GENOME_PAIRS) %>%
  count(species_pair, name = "n_groups") %>%
  filter(n_groups == 2) %>%
  pull(species_pair)

plot_df <- agg_df %>%
  filter(species_pair %in% shared_species_pairs) %>%
  mutate(species_pair = factor(species_pair))

wide_df <- plot_df %>%
  select(species_pair, phenotype_group, transfer_rate) %>%
  pivot_wider(names_from = phenotype_group, values_from = transfer_rate) %>%
  mutate(diff = IBD - `non-IBD`)

plot_long <- plot_df %>%
  left_join(wide_df %>% select(species_pair, diff), by = "species_pair")

model_long <- plot_df %>%
  transmute(
    species_pair = factor(species_pair),
    phylogeny = as.numeric(phylogeny),
    group = factor(phenotype_group, levels = c("IBD", "non-IBD")),
    total = total_genome_pair,
    positive = positive_genome_pair,
    negative = total_genome_pair - positive_genome_pair
  ) %>%
  filter(total > 0, positive >= 0, negative >= 0) %>%
  droplevels()

fit_null <- glmer(
  cbind(positive, negative) ~ phylogeny + (1 | species_pair),
  family = binomial(),
  data = model_long,
  nAGQ = 0,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
)

fit_full <- glmer(
  cbind(positive, negative) ~ group + phylogeny + (1 | species_pair),
  family = binomial(),
  data = model_long,
  nAGQ = 0,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
)

lrt <- lmtest::lrtest(fit_null, fit_full)
lrt_df <- as.data.frame(lrt)
p_col <- grep("Pr\\(", colnames(lrt_df), value = TRUE)
lrt_p <- if (length(p_col)) as.numeric(lrt_df[2, p_col[1]]) else NA_real_
coef_table <- coef(summary(fit_full))
group_term <- "groupnon-IBD"
beta <- coef_table[group_term, "Estimate"]
se <- coef_table[group_term, "Std. Error"]
or <- exp(beta)
or_low <- exp(beta - 1.96 * se)
or_high <- exp(beta + 1.96 * se)

write_csv(model_long, file.path(out_dir, "HGT_IBD_vs_nonIBD_shared_species_pair_gp10_glmm_input.csv"))
write_csv(lrt_df, file.path(out_dir, "HGT_IBD_vs_nonIBD_shared_species_pair_gp10_glmm_lrt.csv"))
write_csv(
  tibble::tibble(
    metric = "HGT",
    term = group_term,
    beta = beta,
    SE = se,
    OR_nonIBD_vs_IBD = or,
    OR_low_95Wald = or_low,
    OR_high_95Wald = or_high,
    LRT_p = lrt_p,
    n_species_pairs = nlevels(model_long$species_pair),
    n_rows = nrow(model_long)
  ),
  file.path(out_dir, "HGT_IBD_vs_nonIBD_shared_species_pair_gp10_glmm_summary.csv")
)
write_csv(plot_long %>% arrange(species_pair, phenotype_group), file.path(out_dir, "HGT_IBD_vs_nonIBD_shared_species_pair_gp10_rates.csv"))

stats_label <- paste0(
  "GLMM LRT p = ", format_p(lrt_p),
  "\nOR non-IBD/IBD = ", format(or, digits = 3),
  " [", format(or_low, digits = 3), ", ", format(or_high, digits = 3), "]"
)

p <- ggplot(plot_long, aes(x = phenotype_group, y = transfer_rate, group = species_pair, color = diff)) +
  geom_line(aes(alpha = abs(diff)), linewidth = 0.7) +
  geom_point(aes(alpha = abs(diff)), size = 2) +
  scale_color_gradient2(low = healthy_color, mid = diff_mid, high = ibd_color, midpoint = 0, guide = "none") +
  scale_alpha(range = c(0.3, 1), guide = "none") +
  scale_x_discrete(labels = c("IBD", "non-IBD")) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(color = c(ibd_color, healthy_color), face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.margin = margin(10, 15, 10, 10)
  ) +
  labs(
    title = "HGT: IBD (left) vs non-IBD (right)",
    subtitle = paste0("Shared species pairs only; both groups genome pair > ", MIN_GENOME_PAIRS, "; n = ", n_distinct(plot_long$species_pair)),
    x = NULL,
    y = "Genome-pair HGT rate"
  ) +
  annotate("text", x = 1.5, y = max(plot_long$transfer_rate, na.rm = TRUE) * 1.05, label = stats_label, size = 3.5, lineheight = 0.9)

ggsave(file.path(out_dir, "HGT_IBD_vs_nonIBD_shared_species_pair_gp10_rate.pdf"), p, width = 6, height = 4.5)
