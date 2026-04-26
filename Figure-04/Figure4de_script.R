#!/usr/bin/env Rscript

base_dir <- "/scratch/p312334/project/10--HGT_isolates/5--Enriched_ME/Plasmid_model_Final/results"
out_dir <- file.path(base_dir, "depth_cohort_comparison", "R_plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

input <- file.path(base_dir, "moduleA_moduleB_vs_all_genomes_annotated.tsv")
df <- read.delim(input, sep = "\t", header = TRUE, quote = "", comment.char = "", check.names = FALSE)

num_cols <- c("qcov", "tcov", "pident", "bitscore", "contig_depth")
for (col in num_cols) {
  df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
}
df <- df[df$Cohort %in% c("DAG3", "RISEUP") & !is.na(df$contig_depth), ]

module_label <- function(x) {
  x <- gsub("01_Module_A_4896bp_n29", "Module A", x)
  x <- gsub("02_Module_B_3683bp_n26", "Module B", x)
  x
}

best_per_genome <- function(x) {
  ord <- order(
    x$qseqid,
    x$Genome_file,
    -x$bitscore,
    -x$pident,
    -x$qcov,
    -x$tcov
  )
  x <- x[ord, ]
  x[!duplicated(paste(x$qseqid, x$Genome_file, sep = "\t")), ]
}

p_label <- function(p) {
  if (is.na(p)) {
    return("p = NA")
  }
  if (p < 0.001) {
    return(sprintf("p = %.2e", p))
  }
  sprintf("p = %.3f", p)
}

summarise_depth <- function(x, filter_name) {
  modules <- sort(unique(x$qseqid))
  rows <- list()
  idx <- 1
  for (module in modules) {
    sub <- x[x$qseqid == module, ]
    vals <- list()
    for (cohort in c("DAG3", "RISEUP")) {
      z <- sub[sub$Cohort == cohort, ]
      d <- z$contig_depth
      vals[[cohort]] <- d
      rows[[idx]] <- data.frame(
        filter = filter_name,
        qseqid = module,
        module = module_label(module),
        cohort = cohort,
        n_contigs = length(d),
        n_genomes = length(unique(z$Genome_file)),
        n_persons = length(unique(z$person_id[!is.na(z$person_id) & z$person_id != ""])),
        n_plasmid_contigs = length(unique(z$subject_id[z$plasmid == 1])),
        median_depth = ifelse(length(d) > 0, median(d), NA),
        mean_depth = ifelse(length(d) > 0, mean(d), NA),
        q1_depth = ifelse(length(d) > 0, as.numeric(quantile(d, 0.25)), NA),
        q3_depth = ifelse(length(d) > 0, as.numeric(quantile(d, 0.75)), NA),
        min_depth = ifelse(length(d) > 0, min(d), NA),
        max_depth = ifelse(length(d) > 0, max(d), NA),
        dag3_median = NA,
        riseup_median = NA,
        median_diff_RISEUP_minus_DAG3 = NA,
        wilcox_w = NA,
        wilcox_p = NA
      )
      idx <- idx + 1
    }
    p <- NA
    w <- NA
    if (length(vals[["DAG3"]]) > 0 && length(vals[["RISEUP"]]) > 0) {
      wt <- suppressWarnings(wilcox.test(vals[["DAG3"]], vals[["RISEUP"]], exact = FALSE))
      p <- wt$p.value
      w <- unname(wt$statistic)
    }
    rows[[idx]] <- data.frame(
      filter = filter_name,
      qseqid = module,
      module = module_label(module),
      cohort = "DAG3_vs_RISEUP",
      n_contigs = length(vals[["DAG3"]]) + length(vals[["RISEUP"]]),
      n_genomes = length(unique(sub$Genome_file)),
      n_persons = length(unique(sub$person_id[!is.na(sub$person_id) & sub$person_id != ""])),
      n_plasmid_contigs = length(unique(sub$subject_id[sub$plasmid == 1])),
      median_depth = NA,
      mean_depth = NA,
      q1_depth = NA,
      q3_depth = NA,
      min_depth = NA,
      max_depth = NA,
      dag3_median = ifelse(length(vals[["DAG3"]]) > 0, median(vals[["DAG3"]]), NA),
      riseup_median = ifelse(length(vals[["RISEUP"]]) > 0, median(vals[["RISEUP"]]), NA),
      median_diff_RISEUP_minus_DAG3 = ifelse(
        length(vals[["DAG3"]]) > 0 && length(vals[["RISEUP"]]) > 0,
        median(vals[["RISEUP"]]) - median(vals[["DAG3"]]),
        NA
      ),
      wilcox_w = w,
      wilcox_p = p
    )
    idx <- idx + 1
  }
  do.call(rbind, rows)
}

plot_depth <- function(x, summary_df, title, out_prefix) {
  modules <- c("01_Module_A_4896bp_n29", "02_Module_B_3683bp_n26")
  colors <- c("DAG3" = "#5EC2D9", "RISEUP" = "#E95D4A")
  plot_one <- function() {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    par(mfrow = c(1, 2), mar = c(4.5, 4.5, 4.2, 1), oma = c(0, 0, 2.5, 0))
    for (module in modules) {
      sub <- x[x$qseqid == module, ]
      d1 <- sub$contig_depth[sub$Cohort == "DAG3"]
      d2 <- sub$contig_depth[sub$Cohort == "RISEUP"]
      y <- c(d1, d2)
      y_max <- ifelse(length(y) > 0, max(y, na.rm = TRUE), 1)
      y_lim <- c(0, y_max * 1.18)
      boxplot(
        list(DAG3 = d1, RISEUP = d2),
        outline = FALSE,
        boxwex = 0.38,
        col = adjustcolor(colors, alpha.f = 0.55),
        border = "#333333",
        lwd = 1.2,
        ylab = "Contig depth",
        main = module_label(module),
        ylim = y_lim
      )
      set.seed(7)
      if (length(d1) > 0) {
        points(jitter(rep(1, length(d1)), amount = 0.075), d1, pch = 16, cex = 0.45, col = adjustcolor(colors["DAG3"], alpha.f = 0.55))
      }
      if (length(d2) > 0) {
        points(jitter(rep(2, length(d2)), amount = 0.075), d2, pch = 16, cex = 0.45, col = adjustcolor(colors["RISEUP"], alpha.f = 0.55))
      }
      s <- summary_df[summary_df$qseqid == module & summary_df$cohort == "DAG3_vs_RISEUP", ]
      p_text <- if (nrow(s) > 0) p_label(s$wilcox_p[1]) else "p = NA"
      n_text <- sprintf("DAG3 n=%d; RISEUP n=%d", length(d1), length(d2))
      text(1.5, y_max * 1.12, p_text, cex = 0.9)
      text(1.5, y_max * 1.04, n_text, cex = 0.78)
      grid(nx = NA, ny = NULL, col = "#e6e6e6")
    }
    mtext(title, outer = TRUE, cex = 1.05, font = 2)
  }
  pdf(file.path(out_dir, paste0(out_prefix, ".pdf")), width = 10, height = 5.3)
  plot_one()
  dev.off()
  png(file.path(out_dir, paste0(out_prefix, ".png")), width = 3000, height = 1600, res = 300)
  plot_one()
  dev.off()
}

filters <- list(
  qcov_gt0.75_tcov_gt0.75_no_pident_filter = list(
    title = "qcov > 0.75 and tcov > 0.75, no pident filter",
    data = df[df$qcov > 0.75 & df$tcov > 0.75, ]
  ),
  qcov_gt0.9_tcov_gt0.9_no_pident_filter = list(
    title = "qcov > 0.9 and tcov > 0.9, no pident filter",
    data = df[df$qcov > 0.9 & df$tcov > 0.9, ]
  ),
  pident_eq100_qcov_gt0.75_tcov_gt0.75 = list(
    title = "pident = 100, qcov > 0.75 and tcov > 0.75",
    data = df[df$pident == 100 & df$qcov > 0.75 & df$tcov > 0.75, ]
  ),
  pident_eq100_qcov_gt0.75 = list(
    title = "pident = 100 and qcov > 0.75",
    data = df[df$pident == 100 & df$qcov > 0.75, ]
  )
)

all_summary <- list()
for (name in names(filters)) {
  x <- best_per_genome(filters[[name]]$data)
  # One row per module/genome is used for all summaries and plots below.
  write.table(x, file.path(out_dir, paste0(name, "_best_per_genome.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  s <- summarise_depth(x, name)
  write.table(s, file.path(out_dir, paste0(name, "_depth_summary.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  plot_depth(x, s, filters[[name]]$title, paste0(name, "_depth_boxplot_R"))
  all_summary[[name]] <- s
}

combined <- do.call(rbind, all_summary)
write.table(combined, file.path(out_dir, "combined_depth_summary_R.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
tests <- combined[combined$cohort == "DAG3_vs_RISEUP", ]
write.table(tests, file.path(out_dir, "cohort_depth_comparison_tests_R.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("Wrote R plots and summaries to:", out_dir, "\n")
print(tests[, c("filter", "module", "n_contigs", "n_genomes", "dag3_median", "riseup_median", "median_diff_RISEUP_minus_DAG3", "wilcox_p")])
