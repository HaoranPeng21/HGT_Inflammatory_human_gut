#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(forcats)
  library(ape)
  library(ggplot2)
  library(ggtree)
  library(ggtreeExtra)
  library(ggnewscale)
  library(ggsci)
})

args <- commandArgs(trailingOnly = TRUE)

default_tree <- "/scratch/hb-fu/haoran/HGT_tree/gtdbtk_all_v3_rooted/output_hpc/infer/gtdbtk.bac120.decorated.tree"
default_meta <- "/scratch/p312334/project/10--HGT_isolates/data/__binstoget_isolates_summary_v4.csv"
default_outdir <- "/scratch/p312334/project/10--HGT_isolates/1--General_infomation/phylotree/result"
default_prefix <- "gtdbtk_all_v3_rooted_Phylum_rings_no_cohort"

tree_path <- if (length(args) >= 1) args[1] else default_tree
meta_path <- if (length(args) >= 2) args[2] else default_meta
outdir <- if (length(args) >= 3) args[3] else default_outdir
prefix <- if (length(args) >= 4) args[4] else default_prefix

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(tree_path)) {
  stop("Tree file not found: ", tree_path)
}
if (!file.exists(meta_path)) {
  stop("Metadata file not found: ", meta_path)
}

clean_chr <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x[x %in% c("", "NA", "na", "NaN", "nan")] <- NA_character_
  x
}

ring_palette_binary <- function() {
  c("0" = "white", "1" = "#82A8D8")
}

phylum_palette <- function(levels_vec) {
  fixed <- c(
    "Bacillota" = "#c9b6ff",
    "Actinomycetota" = "#7fd3e6",
    "Bacteroidota" = "#f3c46b",
    "Pseudomonadota" = "#f28cb1",
    "Desulfobacterota" = "#8fd19e",
    "Fusobacteriota" = "#f6a66a"
  )
  base_cols <- unique(c(
    fixed,
    pal_npg("nrc")(10),
    pal_jco("default")(10),
    pal_d3("category20")(20),
    pal_igv("default")(51)
  ))
  if (length(levels_vec) <= length(base_cols)) {
    cols <- base_cols[seq_along(levels_vec)]
  } else {
    cols <- grDevices::colorRampPalette(base_cols)(length(levels_vec))
  }
  names(cols) <- levels_vec
  cols[names(fixed)] <- fixed[names(fixed)]
  cols
}

tree <- read.tree(tree_path)
tree$tip.label <- tree$tip.label %>%
  str_remove("^'") %>%
  str_remove("'$") %>%
  str_remove("\\.fa$")

metadata <- read_csv(meta_path, show_col_types = FALSE) %>%
  transmute(
    label = clean_chr(Genome_file) %>% str_remove("\\.fa$"),
    Phylum = clean_chr(Phylum),
    plasmid = clean_chr(plasmid),
    HGT_within = clean_chr(HGT_within)
  ) %>%
  distinct(label, .keep_all = TRUE)

matched_labels <- intersect(tree$tip.label, metadata$label)
if (length(matched_labels) < 10) {
  stop("Too few tree tips matched metadata: ", length(matched_labels))
}

tree <- keep.tip(tree, matched_labels)

tip_df <- tibble(label = tree$tip.label) %>%
  left_join(metadata, by = "label") %>%
  mutate(
    Phylum = fct_infreq(Phylum),
    plasmid = case_when(
      is.na(plasmid) ~ NA_character_,
      suppressWarnings(as.numeric(plasmid)) > 0 ~ "1",
      TRUE ~ "0"
    ),
    plasmid = factor(plasmid, levels = c("0", "1")),
    HGT_within = factor(HGT_within, levels = c("0", "1"))
  )

matched_n <- sum(!is.na(tip_df$Phylum))

phylum_levels <- levels(tip_df$Phylum)
plasmid_levels <- levels(tip_df$plasmid)

ring_width <- 0.018
ring_tile_width <- 0.10
ring_offset_first <- 0.060
ring_offset_step <- 0.040
plasmid_cols <- c("0" = "white", "1" = "#8FB389")
hgt_within_cols <- ring_palette_binary()
phylum_cols <- phylum_palette(phylum_levels)

phylum_genomes_list <- tip_df %>%
  filter(!is.na(Phylum)) %>%
  group_by(Phylum) %>%
  summarise(genomes = list(label), .groups = "drop") %>%
  tibble::deframe()

tree_grouped <- groupOTU(tree, phylum_genomes_list)

ring_plasmid <- tip_df %>% transmute(label, ring = "plasmid", value = plasmid)
ring_hgt_within <- tip_df %>% transmute(label, ring = "HGT_within", value = HGT_within)

p <- ggtree(
    tree_grouped,
    aes(color = group),
    layout = "fan",
    open.angle = 24,
    branch.length = "branch.length",
    linewidth = 0.28
  ) %<+% tip_df +
  scale_color_manual(values = phylum_cols, na.value = "grey70", name = "Phylum") +
  geom_treescale(x = 0.02, y = 0.02, width = 0.10, fontsize = 3.2, linesize = 0.8, color = "grey30") +
  geom_fruit(
    data = ring_plasmid,
    geom = geom_tile,
    mapping = aes(y = label, x = ring, fill = value),
    color = NA,
    width = ring_tile_width,
    offset = ring_offset_step,
    pwidth = ring_width
  ) +
  scale_fill_manual(values = plasmid_cols, na.value = "white", name = "plasmid") +
  new_scale_fill() +
  geom_fruit(
    data = ring_hgt_within,
    geom = geom_tile,
    mapping = aes(y = label, x = ring, fill = value),
    color = NA,
    width = ring_tile_width,
    offset = ring_offset_step,
    pwidth = ring_width
  ) +
  scale_fill_manual(values = hgt_within_cols, na.value = "white", name = "HGT_within") +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 9),
    plot.margin = margin(8, 8, 8, 8),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(linewidth = 1.2, size = 3, alpha = 1)),
    fill = guide_legend(order = 2)
  )

pdf_path <- file.path(outdir, paste0(prefix, ".pdf"))
png_path <- file.path(outdir, paste0(prefix, ".png"))

ggsave(pdf_path, p, width = 20, height = 20, units = "in", device = cairo_pdf, limitsize = FALSE)
ggsave(png_path, p, width = 20, height = 20, units = "in", dpi = 300, limitsize = FALSE)

message("Matched tips with Phylum annotation: ", matched_n, " / ", length(tree$tip.label))
message("PDF written to: ", pdf_path)
message("PNG written to: ", png_path)
