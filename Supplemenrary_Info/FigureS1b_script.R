#!/usr/bin/env Rscript

base_dir <- "/scratch/p312334/project/10--HGT_isolates/1--General_infomation/conserve_cutoff"
result_dir <- file.path(base_dir, "result")
input_tsv <- file.path(result_dir, "ALL_speciespair_ribosomal42_snp.tsv")
output_pdf <- file.path(result_dir, "ALL_speciespair_ribosomal42_snp_histogram_0_10_R.pdf")
output_png <- file.path(result_dir, "ALL_speciespair_ribosomal42_snp_histogram_0_10_R.png")

df <- read.delim(input_tsv, sep = "\t", header = TRUE, check.names = FALSE)
values <- df$snp_per_2000bp
values <- values[!is.na(values) & values >= 0 & values <= 10]

breaks <- seq(-0.5, 10.5, by = 1)
hist_obj <- hist(values, breaks = breaks, plot = FALSE, right = FALSE)
y_max <- max(hist_obj$counts)
y_top <- if (y_max > 0) ceiling(y_max * 1.12) else 1

draw_plot <- function() {
  par(
    family = "Helvetica",
    mar = c(4.6, 4.8, 2.2, 0.8),
    mgp = c(2.6, 0.7, 0),
    las = 1,
    bty = "l"
  )

  plot(
    NA,
    xlim = c(-0.5, 10.5),
    ylim = c(0, y_top),
    xaxs = "i",
    yaxs = "i",
    xaxt = "n",
    yaxt = "n",
    xlab = "Pairwise SNP per 2000 bp",
    ylab = "Number of Species Pairs"
  )

  y_ticks <- pretty(c(0, y_top))
  y_ticks <- y_ticks[y_ticks >= 0 & y_ticks <= y_top]
  abline(h = y_ticks, col = "#E3E3E3", lwd = 0.8)

  rect(
    xleft = hist_obj$breaks[-length(hist_obj$breaks)],
    ybottom = 0,
    xright = hist_obj$breaks[-1],
    ytop = hist_obj$counts,
    col = "#8A8A8A",
    border = "#6A6A6A",
    lwd = 0.8
  )

  axis(1, at = 0:10, labels = 0:10, col.axis = "#4F4F4F", col = "#777777")
  axis(2, at = y_ticks, labels = y_ticks, col.axis = "#4F4F4F", col = "#777777")
  box(bty = "l", col = "#777777")
}

pdf(output_pdf, width = 6.2, height = 5.5, useDingbats = FALSE)
draw_plot()
dev.off()

png(output_png, width = 6.2, height = 5.5, units = "in", res = 200, type = "cairo")
draw_plot()
dev.off()

cat(sprintf("[INFO] Wrote PDF: %s\n", output_pdf))
cat(sprintf("[INFO] Wrote PNG: %s\n", output_png))
