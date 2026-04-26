#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggsci)
})

base_dir <- "/scratch/p312334/project/10--HGT_isolates/3.2--Pangenome"
summary_file <- file.path(base_dir, "species_fdr05_gene_hgt_summary.tsv")
fisher_file <- file.path(base_dir, "species_fdr05_hgt_fisher.tsv")
out_pdf <- file.path(base_dir, "species_fdr05_hgt_barplot_log10.pdf")
out_png <- file.path(base_dir, "species_fdr05_hgt_barplot_log10.png")
out_data <- file.path(base_dir, "species_fdr05_hgt_barplot_log10_data.tsv")

dat <- read.delim(summary_file, check.names = FALSE, stringsAsFactors = FALSE)
fisher <- read.delim(fisher_file, check.names = FALSE, stringsAsFactors = FALSE)

dat <- dat[dat$n_fdr05_genes > 10, ]
dat <- merge(
  dat,
  fisher[, c("species", "fisher_bh_fdr", "bh_fdr_lt_0.05")],
  by = "species",
  all.x = TRUE,
  sort = FALSE
)

dat$species_label <- sub("__([0-9]+)_([0-9]+)$", " (\\1_\\2)", dat$species)
dat$species_label <- gsub("_", " ", dat$species_label)
dat$species_label <- gsub(" \\(([0-9]+) ([0-9]+)\\)$", " (\\1_\\2)", dat$species_label)
dat$species_label <- ifelse(
  dat$bh_fdr_lt_0.05 == "TRUE",
  paste0(dat$species_label, " *"),
  dat$species_label
)

order_species <- dat$species_label[order(dat$n_fdr05_genes, decreasing = TRUE)]
dat$species_label <- factor(dat$species_label, levels = rev(order_species))

plot_dat <- rbind(
  data.frame(
    species = dat$species,
    species_label = dat$species_label,
    direction = "Health",
    class = "All significant",
    count = dat$n_health_sig,
    total_count = dat$n_health_sig,
    hgt_count = dat$n_health_sig_in_hgt,
    stringsAsFactors = FALSE
  ),
  data.frame(
    species = dat$species,
    species_label = dat$species_label,
    direction = "IBD",
    class = "All significant",
    count = dat$n_ibd_sig,
    total_count = dat$n_ibd_sig,
    hgt_count = dat$n_ibd_sig_in_hgt,
    stringsAsFactors = FALSE
  ),
  data.frame(
    species = dat$species,
    species_label = dat$species_label,
    direction = "Health",
    class = "HGT subset",
    count = dat$n_health_sig_in_hgt,
    total_count = dat$n_health_sig,
    hgt_count = dat$n_health_sig_in_hgt,
    stringsAsFactors = FALSE
  ),
  data.frame(
    species = dat$species,
    species_label = dat$species_label,
    direction = "IBD",
    class = "HGT subset",
    count = dat$n_ibd_sig_in_hgt,
    total_count = dat$n_ibd_sig,
    hgt_count = dat$n_ibd_sig_in_hgt,
    stringsAsFactors = FALSE
  )
)

plot_dat$direction <- factor(plot_dat$direction, levels = c("Health", "IBD"))
plot_dat$class <- factor(plot_dat$class, levels = c("All significant", "HGT subset"))
plot_dat$log10_count_plus1 <- log10(plot_dat$count + 1)
plot_dat$bar_type <- paste(plot_dat$class, plot_dat$direction)
plot_dat$bar_type <- factor(
  plot_dat$bar_type,
  levels = c("All significant Health", "HGT subset Health", "All significant IBD", "HGT subset IBD"),
  labels = c("Sig Health", "HGT Health", "Sig IBD", "HGT IBD")
)
plot_dat$hgt_percent <- ifelse(plot_dat$total_count > 0, 100 * plot_dat$hgt_count / plot_dat$total_count, NA_real_)
plot_dat$percent_label <- ifelse(
  is.na(plot_dat$hgt_percent),
  "NA",
  sprintf("%d/%d (%.1f%%)", plot_dat$hgt_count, plot_dat$total_count, plot_dat$hgt_percent)
)

write.table(plot_dat, out_data, sep = "\t", quote = FALSE, row.names = FALSE)

cols <- ggsci::pal_npg()(10)
bar_cols <- c(
  "Sig Health" = "#C9ECF3",
  "HGT Health" = cols[2],
  "Sig IBD" = "#F6CAC4",
  "HGT IBD" = cols[1]
)

total_dat <- plot_dat[plot_dat$class == "All significant", ]
hgt_dat <- plot_dat[plot_dat$class == "HGT subset", ]

p <- ggplot() +
  geom_col(
    data = total_dat,
    aes(x = species_label, y = log10_count_plus1, fill = bar_type, group = direction),
    position = position_dodge2(width = 0.78, preserve = "single"),
    width = 0.70,
    alpha = 1,
    color = NA
  ) +
  geom_col(
    data = hgt_dat,
    aes(x = species_label, y = log10_count_plus1, fill = bar_type, group = direction),
    position = position_dodge2(width = 0.78, preserve = "single"),
    width = 0.70,
    alpha = 1,
    color = "white",
    linewidth = 0.18
  ) +
  geom_text(
    data = total_dat,
    aes(
      x = species_label,
      y = log10_count_plus1,
      label = percent_label,
      group = direction,
      color = direction
    ),
    position = position_dodge2(width = 0.78, preserve = "single"),
    hjust = -0.08,
    size = 2.8,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = bar_cols, breaks = names(bar_cols), name = NULL) +
  scale_color_manual(values = c(Health = cols[2], IBD = cols[1]), guide = "none") +
  scale_y_continuous(
    name = expression("Gene count  [" * log[10]("(count + 1) scale]")),
    breaks = log10(c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000) + 1),
    labels = c("1", "2", "5", "10", "20", "50", "100", "200", "500", "1000", "2000"),
    expand = expansion(mult = c(0, 0.15))
  ) +
  coord_flip(clip = "off") +
  labs(
    x = NULL,
    title = "FDR < 0.05 pangenome genes and HGT-overlapping subset",
    subtitle = "Pale = significant genes; solid overlay = HGT subset; labels = HGT/Significant (%); * Fisher BH-FDR < 0.05"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, color = "grey30"),
    axis.text.y = element_text(size = 9, color = "grey10"),
    axis.text.x = element_text(size = 9, color = "grey20"),
    axis.title.x = element_text(size = 10, margin = margin(t = 8)),
    legend.position = "top",
    legend.justification = "left",
    legend.text = element_text(size = 9),
    panel.grid.major.x = element_line(color = "grey88", linewidth = 0.25),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(8, 46, 8, 8)
  )

ggsave(out_pdf, p, width = 10.8, height = 5.8, device = cairo_pdf)
ggsave(out_png, p, width = 10.8, height = 5.8, dpi = 320)

message(out_pdf)
message(out_png)
