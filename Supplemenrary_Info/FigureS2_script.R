library(ggplot2)
library(scales)

folder <- "/scratch/p312334/project/10--HGT_isolates/1--General_infomation/Fixation_100/directional_prevalence_v4.1_skani9999"
infile <- file.path(folder, "HGT_table_2000bp_isolates_within_plasmid.gene_individual_prevalence.directional.csv")
hist_out <- file.path(folder, "HGT_individual_gene_prevalence_total_genome_ge5_histogram.no1.pdf")
pie_out <- file.path(folder, "HGT_individual_gene_prevalence_total_genome_ge5_eq1_pie.pdf")

dat <- read.csv(infile, stringsAsFactors = FALSE)
dat$total_genome_count <- suppressWarnings(as.numeric(dat$total_genome_count))
dat$individual_gene_prev <- suppressWarnings(as.numeric(dat$individual_gene_prev))

plot_dat <- dat[
  !is.na(dat$total_genome_count) &
    dat$total_genome_count >= 5 &
    !is.na(dat$individual_gene_prev) &
    dat$individual_gene_prev >= 0 &
    dat$individual_gene_prev <= 1,
]

no1_dat <- plot_dat[plot_dat$individual_gene_prev < 1, ]

p_hist <- ggplot(no1_dat, aes(x = individual_gene_prev)) +
  geom_histogram(
    binwidth = 0.02,
    boundary = 0,
    fill = "#BDBDBD",
    color = "#FFFFFF",
    linewidth = 0.2
  ) +
  scale_x_continuous(
    breaks = seq(0, 1, 0.1),
    labels = percent_format(accuracy = 1),
    expand = c(0, 0)
  ) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = "Individual-level gene prevalence (excluding 100%)", y = "Count") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#D9D9D9", linewidth = 0.3),
    axis.line.x = element_line(color = "black", linewidth = 0.4),
    axis.line.y = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.margin = margin(8, 12, 8, 8)
  )

ggsave(hist_out, plot = p_hist, width = 7, height = 5, device = cairo_pdf)

pie_dat <- data.frame(
  label = c("Prevalence = 100%", "Prevalence < 100%"),
  count = c(sum(plot_dat$individual_gene_prev == 1), sum(plot_dat$individual_gene_prev < 1))
)
pie_dat$fraction <- pie_dat$count / sum(pie_dat$count)
pie_dat$ypos <- cumsum(pie_dat$fraction) - 0.5 * pie_dat$fraction
pie_dat$label_text <- sprintf("%s\n%s", pie_dat$label, percent(pie_dat$fraction, accuracy = 0.1))

p_pie <- ggplot(pie_dat, aes(x = "", y = count, fill = label)) +
  geom_col(width = 1, color = "white", linewidth = 0.4) +
  coord_polar(theta = "y") +
  geom_text(aes(y = ypos * sum(count), label = label_text), color = "white", size = 4) +
  scale_fill_manual(values = c("Prevalence = 100%" = "#969696", "Prevalence < 100%" = "#D9D9D9")) +
  labs(fill = NULL) +
  theme_void(base_size = 12) +
  theme(legend.position = "none")

ggsave(pie_out, plot = p_pie, width = 5.5, height = 5, device = cairo_pdf)

cat("wrote:", hist_out, "\n")
cat("rows excluding 1:", nrow(no1_dat), "\n")
cat("min/median/max (excluding 1):", min(no1_dat$individual_gene_prev), median(no1_dat$individual_gene_prev), max(no1_dat$individual_gene_prev), "\n")
cat("wrote:", pie_out, "\n")
cat("rows total:", nrow(plot_dat), "\n")
cat("rows eq 1:", pie_dat$count[pie_dat$label == "Prevalence = 100%"], "\n")
cat("rows lt 1:", pie_dat$count[pie_dat$label == "Prevalence < 100%"], "\n")
