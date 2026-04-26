#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
})

base_dir <- "/scratch/p312334/project/10--HGT_isolates/3.2--Pangenome"
detail_path <- file.path(base_dir, "species_fdr05_gene_hgt_detail.tsv")
species_to_plot <- c(
  "Bacteroides_uniformis__48_2",
  "Bacteroides_fragilis__32_1",
  "Bifidobacterium_adolescentis__335_1"
)

cols <- c(
  `IBD sig (beta > 0)` = "#E64B35FF",
  `Health sig (beta < 0)` = "#4DBBD5FF",
  `Not FDR < 0.05` = "#C9C9C9"
)
hgt_fill <- "#F2C14E"

find_header_line <- function(path) {
  con <- file(path, open = "r")
  on.exit(close(con))
  i <- 0L
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) return(NA_integer_)
    i <- i + 1L
    if (startsWith(line, "variant\t")) return(i)
  }
}

drop_plot_beta_outliers <- function(dat) {
  if (nrow(dat) < 10) return(dat)
  absb <- abs(dat$beta)
  q99 <- suppressWarnings(as.numeric(quantile(absb, probs = 0.99, na.rm = TRUE, names = FALSE)))
  if (!is.finite(q99) || q99 <= 0) return(dat)
  cutoff <- max(50, 20 * q99)
  keep <- !(absb > cutoff & dat$`lrt-pvalue` >= 0.99)
  dat[keep, , drop = FALSE]
}

read_pyseer_cogs <- function(path) {
  header_line <- find_header_line(path)
  if (is.na(header_line)) stop("Could not find pyseer header in: ", path)
  dat <- read.delim(
    path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    skip = header_line - 1L,
    check.names = FALSE
  )
  required_cols <- c("variant", "af", "lrt-pvalue", "beta")
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    stop("Missing columns in ", path, ": ", paste(missing_cols, collapse = ", "))
  }
  dat$af <- suppressWarnings(as.numeric(dat$af))
  dat$`lrt-pvalue` <- suppressWarnings(as.numeric(dat$`lrt-pvalue`))
  dat$beta <- suppressWarnings(as.numeric(dat$beta))
  dat <- dat[
    is.finite(dat$af) &
      is.finite(dat$`lrt-pvalue`) &
      dat$`lrt-pvalue` > 0 &
      is.finite(dat$beta),
    ,
    drop = FALSE
  ]
  if ("notes" %in% names(dat)) {
    dat <- dat[is.na(dat$notes) | dat$notes == "", , drop = FALSE]
  }
  dat <- drop_plot_beta_outliers(dat)
  dat$logp <- -log10(dat$`lrt-pvalue`)
  dat
}

read_table_if_present <- function(path) {
  info <- file.info(path)
  if (is.na(info$size) || info$size == 0) return(NULL)
  read.delim(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
}

get_cutoff <- function(path, p_col, label) {
  dat <- read_table_if_present(path)
  if (is.null(dat) || !p_col %in% names(dat)) return(NULL)
  p <- suppressWarnings(as.numeric(dat[[p_col]]))
  p <- p[is.finite(p) & p > 0]
  if (length(p) == 0) return(NULL)
  cutoff_p <- max(p)
  data.frame(
    label = label,
    cutoff_p = cutoff_p,
    cutoff_logp = -log10(cutoff_p),
    stringsAsFactors = FALSE
  )
}

read_hgt_label_summary <- function(path) {
  dat <- read_table_if_present(path)
  if (is.null(dat) || !all(c("group", "group_annotation") %in% names(dat))) {
    return(data.frame(group = character(), label = character(), stringsAsFactors = FALSE))
  }
  dat$label <- ifelse(
    is.na(dat$group_annotation) | dat$group_annotation == "",
    dat$group,
    dat$group_annotation
  )
  dat[, c("group", "label"), drop = FALSE]
}

wrap_label <- function(x, width = 34) {
  vapply(
    x,
    function(s) paste(strwrap(s, width = width), collapse = "\n"),
    character(1)
  )
}

stagger_label_positions <- function(label_dat, min_gap = 0.46, x_pad = 0.45, y_pad = 0.08) {
  if (nrow(label_dat) == 0) return(label_dat)

  place_side <- function(df, side = c("right", "left")) {
    side <- match.arg(side)
    if (nrow(df) == 0) return(df)
    df <- df[order(df$logp, decreasing = TRUE), , drop = FALSE]
    target_y <- df$logp + y_pad
    gap <- min_gap
    if (nrow(df) > 1L) {
      available_y <- max(target_y, na.rm = TRUE) - 0.45
      if (is.finite(available_y) && available_y > 0) {
        gap <- min(min_gap, available_y / (nrow(df) - 1L))
      }
    }
    label_y <- numeric(nrow(df))
    for (i in seq_len(nrow(df))) {
      if (i == 1L) {
        label_y[i] <- target_y[i]
      } else {
        label_y[i] <- min(target_y[i], label_y[i - 1L] - gap)
      }
    }
    df$label_y <- pmax(label_y, 0.35)
    df$label_x <- if (side == "right") max(df$beta, na.rm = TRUE) + x_pad else min(df$beta, na.rm = TRUE) + x_pad
    df$hjust <- 0
    df
  }

  out <- rbind(
    place_side(label_dat[label_dat$beta >= 0, , drop = FALSE], "right"),
    place_side(label_dat[label_dat$beta < 0, , drop = FALSE], "left")
  )
  out[order(match(out$variant, label_dat$variant)), , drop = FALSE]
}

species_label <- function(species) {
  parts <- strsplit(species, "__", fixed = TRUE)[[1]]
  if (length(parts) != 2) return(gsub("_", " ", species))
  paste0(gsub("_", " ", parts[1]), " (", parts[2], ")")
}

short_species_label <- function(species) {
  parts <- strsplit(species, "__", fixed = TRUE)[[1]]
  if (length(parts) != 2) return(gsub("_", " ", species))
  taxa <- strsplit(parts[1], "_", fixed = TRUE)[[1]]
  if (length(taxa) >= 2) {
    paste0(substr(taxa[1], 1, 1), ". ", taxa[2], " (", parts[2], ")")
  } else {
    paste0(gsub("_", " ", parts[1]), " (", parts[2], ")")
  }
}

make_plot <- function(species, detail_dat) {
  species_dir <- file.path(base_dir, species)
  pyseer_path <- file.path(species_dir, "IBD_COGs.txt")
  fdr_path <- file.path(species_dir, "IBD_COGs_FDR_lt_0.05.tsv")
  top50_path <- file.path(species_dir, "IBD_COGs_FDR_lt_0.05_top50.tsv")
  top50_hgt_summary_path <- file.path(species_dir, "top50_FDR_lt_0.05_within_HGT_overlap_summary.tsv")
  out_pdf <- file.path(species_dir, "IBD_COGs_top50_FDR_HGT_highlight.pdf")
  out_png <- file.path(species_dir, "IBD_COGs_top50_FDR_HGT_highlight.png")
  out_portrait_pdf <- file.path(species_dir, "IBD_COGs_top50_FDR_HGT_highlight_portrait.pdf")
  out_portrait_png <- file.path(species_dir, "IBD_COGs_top50_FDR_HGT_highlight_portrait.png")
  out_4x6_pdf <- file.path(species_dir, "IBD_COGs_top50_FDR_HGT_highlight_4x6.pdf")
  out_4x6_png <- file.path(species_dir, "IBD_COGs_top50_FDR_HGT_highlight_4x6.png")
  out_2x3_bigfont_pdf <- file.path(species_dir, "IBD_COGs_top50_FDR_HGT_highlight_2x3_bigfont.pdf")
  out_2x3_bigfont_png <- file.path(species_dir, "IBD_COGs_top50_FDR_HGT_highlight_2x3_bigfont.png")

  dat <- read_pyseer_cogs(pyseer_path)
  sig_dat <- detail_dat[detail_dat$species == species, , drop = FALSE]
  sig_genes <- unique(sig_dat$gene)
  hgt_genes <- unique(sig_dat$gene[toupper(sig_dat$in_hgt_region) == "TRUE"])

  dat$is_sig <- dat$variant %in% sig_genes
  dat$is_hgt_fdr05 <- dat$variant %in% hgt_genes
  dat$point_class <- ifelse(
    dat$is_sig & dat$beta > 0,
    "IBD sig (beta > 0)",
    ifelse(
      dat$is_sig & dat$beta < 0,
      "Health sig (beta < 0)",
      "Not FDR < 0.05"
    )
  )
  dat$point_class <- factor(dat$point_class, levels = names(cols))

  hgt_dat <- dat[dat$is_hgt_fdr05, , drop = FALSE]
  hgt_dat <- merge(
    hgt_dat,
    sig_dat[, c("gene", "fdr", "direction", "max_overlap_bp")],
    by.x = "variant",
    by.y = "gene",
    all.x = TRUE,
    sort = FALSE
  )

  top50_cutoff <- get_cutoff(top50_path, "lrt-pvalue", "Top50 cutoff")
  fdr_cutoff <- get_cutoff(fdr_path, "lrt-pvalue", "FDR 0.05 cutoff")
  cutoff_dat <- rbind(top50_cutoff, fdr_cutoff)
  cutoff_dat$line_type <- factor(cutoff_dat$label, levels = c("Top50 cutoff", "FDR 0.05 cutoff"))

  label_summary <- read_hgt_label_summary(top50_hgt_summary_path)
  label_dat <- merge(
    dat,
    label_summary,
    by.x = "variant",
    by.y = "group",
    all = FALSE,
    sort = FALSE
  )
  label_dat$plot_label <- wrap_label(label_dat$label, width = 44)
  label_dat <- stagger_label_positions(label_dat)

  n_sig <- sum(dat$is_sig)
  n_hgt <- nrow(hgt_dat)
  n_top50_hgt <- sum(
    hgt_dat$variant %in%
      read_table_if_present(top50_path)$variant
  )

  p <- ggplot(dat, aes(x = beta, y = logp)) +
    geom_hline(
      data = cutoff_dat,
      inherit.aes = FALSE,
      aes(yintercept = cutoff_logp, linetype = line_type),
      colour = "#222222",
      linewidth = 0.42
    ) +
    geom_vline(xintercept = 0, colour = "#808080", linewidth = 0.35, linetype = "solid") +
    geom_point(
      data = hgt_dat,
      inherit.aes = FALSE,
      aes(x = beta, y = logp, fill = "HGT overlap"),
      shape = 21,
      size = 3.8,
      stroke = 0,
      colour = hgt_fill,
      alpha = 0.28,
      show.legend = TRUE
    ) +
    geom_point(aes(colour = point_class), alpha = 0.30, size = 1.35, stroke = 0) +
    geom_segment(
      data = label_dat,
      inherit.aes = FALSE,
      aes(x = beta, y = logp, xend = label_x, yend = label_y),
      colour = "#6E6E6E",
      linewidth = 0.22,
      alpha = 0.72,
      show.legend = FALSE
    ) +
    geom_text(
      data = label_dat,
      inherit.aes = FALSE,
      aes(x = label_x, y = label_y, label = plot_label, hjust = hjust),
      size = 1.75,
      colour = "#111111",
      lineheight = 0.9,
      show.legend = FALSE
    ) +
    geom_label(
      data = cutoff_dat,
      inherit.aes = FALSE,
      aes(x = -Inf, y = cutoff_logp, label = label),
      hjust = -0.05,
      vjust = c(-0.2, 1.2)[seq_len(nrow(cutoff_dat))],
      size = 2.45,
      linewidth = 0,
      label.padding = unit(0.12, "lines"),
      fill = "white",
      alpha = 0.88,
      colour = "#222222",
      show.legend = FALSE
    ) +
    scale_colour_manual(
      values = cols,
      name = NULL
    ) +
    scale_linetype_manual(
      values = c("Top50 cutoff" = "dashed", "FDR 0.05 cutoff" = "dotted"),
      name = NULL
    ) +
    scale_fill_manual(
      values = c("HGT overlap" = hgt_fill),
      name = NULL
    ) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
    scale_x_continuous(expand = expansion(mult = c(0.04, 0.45))) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 11) +
    theme(
      plot.margin = margin(8, 24, 8, 10),
      plot.title = element_text(size = 12.8, face = "bold", colour = "#111111"),
      plot.subtitle = element_text(size = 8.8, colour = "#333333", margin = margin(b = 5)),
      axis.title = element_text(size = 10.5, colour = "#111111"),
      axis.text = element_text(size = 8.8, colour = "#333333"),
      axis.line = element_line(colour = "#111111", linewidth = 0.5),
      axis.ticks = element_line(colour = "#111111", linewidth = 0.4),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_blank(),
      legend.text = element_text(size = 8, colour = "#333333"),
      legend.key.width = grid::unit(1.4, "lines")
    ) +
    guides(
      colour = guide_legend(override.aes = list(alpha = 0.65, size = 2.2), order = 1),
      fill = guide_legend(override.aes = list(alpha = 0.45, size = 3.3, shape = 21, colour = hgt_fill), order = 2),
      linetype = "none"
    ) +
    labs(
      title = paste0(species_label(species), " IBD-associated COGs"),
      subtitle = paste0(
        "Red = beta > 0, blue = beta < 0, gray = not FDR < 0.05; point alpha = 30%\n",
        "Gold halo marks HGT-overlap FDR < 0.05 genes (", n_hgt, "/", n_sig, "); labels show top50 HGT gene annotations"
      ),
      x = "Effect size beta (Health enriched < 0; IBD enriched > 0)",
      y = "-log10(lrt p-value)"
    )

  ggsave(out_pdf, p, width = 11, height = 5.3, device = cairo_pdf)
  ggsave(out_png, p, width = 11, height = 5.3, dpi = 300)

  portrait <- p +
    scale_x_continuous(expand = expansion(mult = c(0.04, 0.62))) +
    theme(
      plot.margin = margin(8, 18, 8, 8),
      plot.title = element_text(size = 11.4, face = "bold", colour = "#111111"),
      plot.subtitle = element_text(size = 7.8, colour = "#333333", margin = margin(b = 4)),
      axis.title = element_text(size = 9.1, colour = "#111111"),
      axis.text = element_text(size = 7.8, colour = "#333333"),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.text = element_text(size = 7.2, colour = "#333333"),
      legend.key.width = grid::unit(1.0, "lines"),
      legend.key.height = grid::unit(0.9, "lines")
    ) +
    guides(
      colour = guide_legend(override.aes = list(alpha = 0.65, size = 1.9), order = 1, nrow = 2, byrow = TRUE),
      fill = guide_legend(override.aes = list(alpha = 0.45, size = 2.6, shape = 21, colour = hgt_fill), order = 2),
      linetype = "none"
    )
  portrait$layers[[3]]$aes_params$size <- 3.1
  portrait$layers[[3]]$aes_params$alpha <- 0.25
  portrait$layers[[4]]$aes_params$size <- 1.05
  portrait$layers[[5]]$aes_params$linewidth <- 0.18
  portrait$layers[[6]]$aes_params$size <- 1.6
  portrait$layers[[7]]$aes_params$size <- 2.05
  ggsave(out_portrait_pdf, portrait, width = 7.5, height = 8.5, device = cairo_pdf)
  ggsave(out_portrait_png, portrait, width = 7.5, height = 8.5, dpi = 300)

  compact_4x6 <- p +
    scale_x_continuous(expand = expansion(mult = c(0.04, 0.78))) +
    theme(
      plot.margin = margin(5, 8, 5, 5),
      plot.title = element_text(size = 10.2, face = "bold", colour = "#111111"),
      plot.subtitle = element_text(size = 7.4, colour = "#333333", margin = margin(b = 3), lineheight = 0.9),
      axis.title = element_text(size = 9.3, colour = "#111111"),
      axis.text = element_text(size = 8.1, colour = "#333333"),
      legend.position = "none"
    ) +
    labs(
      title = paste0(short_species_label(species), " IBD-associated COGs"),
      subtitle = paste0("Red=IBD, blue=Health, gray=NS; gold=HGT (", n_hgt, "/", n_sig, ")"),
      x = "Effect size beta",
      y = "-log10(lrt p-value)"
    ) +
    guides(
      colour = guide_legend(override.aes = list(alpha = 0.7, size = 2.0), order = 1, nrow = 2, byrow = TRUE),
      fill = guide_legend(override.aes = list(alpha = 0.45, size = 2.6, shape = 21, colour = hgt_fill), order = 2),
      linetype = "none"
    )
  compact_4x6$layers[[3]]$aes_params$size <- 2.8
  compact_4x6$layers[[3]]$aes_params$alpha <- 0.25
  compact_4x6$layers[[4]]$aes_params$size <- 0.9
  compact_4x6$layers[[5]]$aes_params$linewidth <- 0.15
  compact_4x6$layers[[6]]$aes_params$size <- 1.55
  compact_4x6$layers[[7]]$aes_params$size <- 1.85
  ggsave(out_4x6_pdf, compact_4x6, width = 4, height = 6, device = cairo_pdf)
  ggsave(out_4x6_png, compact_4x6, width = 4, height = 6, dpi = 300)

  compact_2x3_bigfont <- p +
    scale_x_continuous(expand = expansion(mult = c(0.04, 1.15))) +
    theme(
      plot.margin = margin(3, 5, 3, 3),
      plot.title = element_text(size = 20.4, face = "bold", colour = "#111111"),
      plot.subtitle = element_text(size = 14.8, colour = "#333333", margin = margin(b = 2), lineheight = 0.82),
      axis.title = element_text(size = 18.6, colour = "#111111"),
      axis.text = element_text(size = 16.2, colour = "#333333"),
      legend.position = "none"
    ) +
    labs(
      title = short_species_label(species),
      subtitle = paste0("R=IBD B=Health\nG=NS Au=HGT ", n_hgt, "/", n_sig),
      x = "beta",
      y = "-log10(P)"
    ) +
    guides(
      colour = "none",
      fill = "none",
      linetype = "none"
    )
  compact_2x3_bigfont$layers[[3]]$aes_params$size <- 2.4
  compact_2x3_bigfont$layers[[3]]$aes_params$alpha <- 0.24
  compact_2x3_bigfont$layers[[4]]$aes_params$size <- 0.8
  compact_2x3_bigfont$layers[[5]]$aes_params$linewidth <- 0.12
  compact_2x3_bigfont$layers[[6]]$aes_params$size <- 3.1
  compact_2x3_bigfont$layers[[7]]$aes_params$size <- 3.7
  ggsave(out_2x3_bigfont_pdf, compact_2x3_bigfont, width = 2, height = 3, device = cairo_pdf)
  ggsave(out_2x3_bigfont_png, compact_2x3_bigfont, width = 2, height = 3, dpi = 300)
  message("Wrote ", out_pdf, ", ", out_png, ", ", out_portrait_pdf, ", ", out_portrait_png, ", ", out_4x6_pdf, ", ", out_4x6_png, ", ", out_2x3_bigfont_pdf, " and ", out_2x3_bigfont_png)
}

detail_dat <- read.delim(detail_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
for (species in species_to_plot) {
  make_plot(species, detail_dat)
}
