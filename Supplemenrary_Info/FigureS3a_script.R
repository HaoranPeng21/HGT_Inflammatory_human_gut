#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
  library(lme4)
  library(lmtest)
})

in_path <- "/scratch/p312334/project/10--HGT_isolates/data/secondary_cluster_species_pair_genome_pairs_by_MGS_ID.csv"
phylogeny_path <- "/scratch/p312334/project/10--HGT_isolates/data/Species_pair_distances_ribosomal42_tree.csv"
out_dir <- "/scratch/p312334/project/10--HGT_isolates/1--General_infomation/Between_Within"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

diff_low <- "#7396BA"
diff_mid <- "#D9D9D9"
diff_high <- "#80A481"

df <- read_csv(in_path, show_col_types = FALSE) %>%
  mutate(
    within_HGT_rate = as.numeric(within_HGT_rate),
    between_HGT_rate = as.numeric(between_HGT_rate),
    within_individual_genome_pairs = as.numeric(within_individual_genome_pairs),
    between_individual_genome_pairs = as.numeric(between_individual_genome_pairs)
  )

phylogeny_df <- read_csv(phylogeny_path, show_col_types = FALSE) %>%
  transmute(
    species_pair = as.character(Species_Pair),
    phylogeny = as.numeric(Average_Distance)
  ) %>%
  distinct(species_pair, .keep_all = TRUE)

model_df <- df %>%
  filter(
    within_individual_genome_pairs > 20,
    between_individual_genome_pairs > 20
  ) %>%
  select(
    species_pair,
    within_individual_genome_pairs,
    between_individual_genome_pairs,
    within_HGT,
    between_HGT,
    within_HGT_rate,
    between_HGT_rate
  ) %>%
  left_join(phylogeny_df, by = "species_pair") %>%
  filter(!is.na(phylogeny)) %>%
  transmute(
    species_pair,
    phylogeny,
    within_total = within_individual_genome_pairs,
    between_total = between_individual_genome_pairs,
    within_HGT = as.numeric(within_HGT),
    between_HGT = as.numeric(between_HGT),
    within = within_HGT_rate,
    between = between_HGT_rate,
    diff = between - within
  )

model_long <- model_df %>%
  select(species_pair, phylogeny, within_total, between_total, within_HGT, between_HGT) %>%
  pivot_longer(
    cols = c(within_total, between_total, within_HGT, between_HGT),
    names_to = c("group", ".value"),
    names_pattern = "(within|between)_(total|HGT)"
  ) %>%
  mutate(
    group = factor(group, levels = c("within", "between")),
    non_HGT = total - HGT,
    species_pair = factor(species_pair)
  ) %>%
  filter(total > 0, HGT >= 0, non_HGT >= 0) %>%
  droplevels()

fit_null <- glmer(
  cbind(HGT, non_HGT) ~ phylogeny + (1 | species_pair),
  family = binomial(),
  data = model_long,
  nAGQ = 0,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
)

fit_full <- glmer(
  cbind(HGT, non_HGT) ~ group + phylogeny + (1 | species_pair),
  family = binomial(),
  data = model_long,
  nAGQ = 0,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
)

lrt <- lmtest::lrtest(fit_null, fit_full)
p_col <- grep("Pr\\(", colnames(as.data.frame(lrt)), value = TRUE)
glmm_p <- if (length(p_col)) as.numeric(as.data.frame(lrt)[2, p_col[1]]) else NA_real_
coef_table <- coef(summary(fit_full))
group_term <- "groupbetween"
glmm_beta <- coef_table[group_term, "Estimate"]
glmm_se <- coef_table[group_term, "Std. Error"]
glmm_or <- exp(glmm_beta)
glmm_or_low <- exp(glmm_beta - 1.96 * glmm_se)
glmm_or_high <- exp(glmm_beta + 1.96 * glmm_se)

write_csv(model_long, file.path(out_dir, "secondary_cluster_MGS_ID_within_between_HGT_rate_glmm_input.csv"))
write_csv(as.data.frame(lrt), file.path(out_dir, "secondary_cluster_MGS_ID_within_between_HGT_rate_glmm_lrt.csv"))
write_csv(
  tibble(
    term = group_term,
    beta = glmm_beta,
    SE = glmm_se,
    OR = glmm_or,
    OR_low_95Wald = glmm_or_low,
    OR_high_95Wald = glmm_or_high,
    LRT_p = glmm_p,
    n_species_pairs = nlevels(model_long$species_pair),
    n_rows = nrow(model_long)
  ),
  file.path(out_dir, "secondary_cluster_MGS_ID_within_between_HGT_rate_glmm_summary.csv")
)

plot_df <- df %>%
  filter(
    within_individual_genome_pairs > 20,
    between_individual_genome_pairs > 20
  ) %>%
  semi_join(model_df %>% select(species_pair), by = "species_pair") %>%
  transmute(
    species_pair,
    within = within_HGT_rate,
    between = between_HGT_rate,
    diff = between - within
  ) %>%
  pivot_longer(
    cols = c(within, between),
    names_to = "mode_label",
    values_to = "HGT_rate"
  ) %>%
  mutate(mode_label = factor(mode_label, levels = c("within", "between")))

title <- "Isolates: within (left) vs between (right)"
subtitle <- paste0(
  "MGS_ID; genome pairs > 20; n = ",
  n_distinct(plot_df$species_pair)
)
stats_label <- paste0(
  "GLMM LRT p ", if (isTRUE(glmm_p == 0)) "< 1e-300" else paste0("= ", format(glmm_p, digits = 4)),
  "\nOR between/within = ", format(glmm_or, digits = 3),
  " [", format(glmm_or_low, digits = 3), ", ", format(glmm_or_high, digits = 3), "]"
)

p <- ggplot(plot_df, aes(x = mode_label, y = HGT_rate, group = species_pair, color = diff)) +
  geom_line(aes(alpha = abs(diff)), linewidth = 0.7) +
  geom_point(aes(alpha = abs(diff)), size = 2) +
  scale_color_gradient2(
    low = diff_low,
    mid = diff_mid,
    high = diff_high,
    midpoint = 0,
    limits = c(-1, 1),
    guide = "none"
  ) +
  scale_alpha(range = c(0.3, 1), guide = "none") +
  theme_classic(base_size = 13) +
  labs(title = title, subtitle = subtitle, x = NULL, y = "HGT rate") +
  annotate(
    "text",
    x = 1.5,
    y = max(plot_df$HGT_rate, na.rm = TRUE) * 1.05,
    label = stats_label,
    size = 3.5,
    lineheight = 0.9
  )

png_path <- file.path(out_dir, "secondary_cluster_MGS_ID_within_between_HGT_rate.png")
pdf_path <- file.path(out_dir, "secondary_cluster_MGS_ID_within_between_HGT_rate.pdf")
ggsave(png_path, p, width = 6, height = 4.5, dpi = 300)
ggsave(pdf_path, p, width = 6, height = 4.5)

message("input=", in_path)
message("png=", png_path)
message("pdf=", pdf_path)
message("n_species_pairs=", n_distinct(plot_df$species_pair))
message("glmm_lrt_p=", glmm_p)
message("glmm_or_between_vs_within=", glmm_or)
