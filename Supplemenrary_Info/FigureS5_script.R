#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

out_dir <- "/scratch/p312334/project/10--HGT_isolates/5--Enriched_ME/Plasmid_model_Final/results/secondary_cluster_network_qcov90_tcov90"

nodes_path <- file.path(out_dir, "qcov90_tcov90_secondary_cluster_nodes.tsv")
edges_path <- file.path(out_dir, "qcov90_tcov90_secondary_cluster_edges.tsv")
layout_path <- file.path(out_dir, "qcov90_tcov90_secondary_cluster_networkx_layout.tsv")
pdf_path <- file.path(out_dir, "secondary_cluster_network_qcov90_tcov90_R_editable_text.pdf")
plain_pdf_path <- file.path(out_dir, "secondary_cluster_network_qcov90_tcov90_R_ggsave_device_pdf.pdf")
png_path <- file.path(out_dir, "secondary_cluster_network_qcov90_tcov90_R_editable_text.png")

edge_colors <- c(
  "Module A" = "#d73027",
  "Module B" = "#4575b4",
  "Shared" = "#6a4c93"
)

phylum_colors <- c(
  "Actinomycetota" = "#e76f51",
  "Actinobacteriota" = "#e76f51",
  "Bacteroidota" = "#2a9d8f",
  "Bacillota" = "#f4a261",
  "Firmicutes" = "#f4a261",
  "Pseudomonadota" = "#457b9d"
)

read_tsv <- function(path) {
  read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

nodes <- read_tsv(nodes_path)
edges <- read_tsv(edges_path)
layout <- read_tsv(layout_path)

nodes$secondary_cluster <- as.character(nodes$secondary_cluster)
edges$secondary_cluster_1 <- as.character(edges$secondary_cluster_1)
edges$secondary_cluster_2 <- as.character(edges$secondary_cluster_2)
layout$secondary_cluster <- as.character(layout$secondary_cluster)

nodes_plot <- merge(nodes, layout, by = "secondary_cluster", all.x = TRUE, sort = FALSE)
nodes_plot$label <- ifelse(
  is.na(nodes_plot$species) | nodes_plot$species == "",
  nodes_plot$secondary_cluster,
  paste(nodes_plot$secondary_cluster, nodes_plot$species, sep = "\n")
)
nodes_plot$fill_plot <- ifelse(nodes_plot$phylum %in% names(phylum_colors), nodes_plot$phylum, "Other")

edge_plot <- merge(
  edges,
  layout,
  by.x = "secondary_cluster_1",
  by.y = "secondary_cluster",
  all.x = TRUE,
  sort = FALSE
)
names(edge_plot)[names(edge_plot) == "x"] <- "x"
names(edge_plot)[names(edge_plot) == "y"] <- "y"
edge_plot <- merge(
  edge_plot,
  layout,
  by.x = "secondary_cluster_2",
  by.y = "secondary_cluster",
  all.x = TRUE,
  sort = FALSE,
  suffixes = c("", "end")
)
names(edge_plot)[names(edge_plot) == "xend"] <- "xend"
names(edge_plot)[names(edge_plot) == "yend"] <- "yend"
edge_plot$edge_width <- ifelse(edge_plot$shared_modules > 1, 2.6, 1.2)
edge_plot$edge_alpha <- ifelse(edge_plot$edge_class == "Shared", 0.68, 0.45)

p <- ggplot() +
  geom_segment(
    data = edge_plot,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend,
      color = edge_class,
      linewidth = edge_width,
      alpha = edge_alpha
    ),
    lineend = "round",
    show.legend = c(color = TRUE, linewidth = FALSE, alpha = FALSE)
  ) +
  geom_point(
    data = nodes_plot,
    aes(x = x, y = y, fill = fill_plot),
    shape = 21,
    size = 8.8,
    color = "#222222",
    stroke = 0.9
  ) +
  geom_text(
    data = nodes_plot,
    aes(x = x, y = y, label = label),
    size = 2.46,
    lineheight = 0.86,
    family = "sans",
    color = "#111111"
  ) +
  scale_color_manual(values = edge_colors, breaks = names(edge_colors), name = "Edges") +
  scale_fill_manual(values = c(phylum_colors, "Other" = "#bdbdbd"), name = "Phylum") +
  scale_linewidth_identity() +
  scale_alpha_identity() +
  coord_equal(clip = "off") +
  labs(title = "Secondary cluster network: qcov > 0.9 and tcov > 0.9") +
  theme_void(base_family = "sans") +
  theme(
    plot.title = element_text(size = 13, hjust = 0.5, margin = margin(b = 8)),
    legend.position = "right",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    plot.margin = margin(12, 18, 12, 12)
  )

ggsave(
  filename = plain_pdf_path,
  plot = p,
  width = 12.5,
  height = 9.5,
  device = "pdf",
  useDingbats = FALSE
)

ggsave(
  filename = pdf_path,
  plot = p,
  width = 12.5,
  height = 9.5,
  device = cairo_pdf,
  family = "sans"
)

ggsave(
  filename = png_path,
  plot = p,
  width = 12.5,
  height = 9.5,
  dpi = 300,
  device = "png"
)
