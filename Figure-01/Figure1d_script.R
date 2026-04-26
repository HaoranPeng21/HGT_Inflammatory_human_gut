freq_file <- Sys.getenv("FREQ_FILE", "/scratch/p312334/project/10--HGT_isolates/1--General_infomation/function/result/COG24_CATEGORY_species_pair_frequency.csv")
annotated_freq_file <- Sys.getenv("ANNOTATED_FREQ_FILE", "/scratch/p312334/project/10--HGT_isolates/1--General_infomation/function/result/COG24_CATEGORY_species_pair_frequency_annotated.csv")
plot_file <- Sys.getenv("PLOT_FILE", "/scratch/p312334/project/10--HGT_isolates/1--General_infomation/function/result/COG24_CATEGORY_species_pair_frequency_vertical_lollipop_colors.pdf")

suppressPackageStartupMessages(library(ggplot2))

cat_file <- "/scratch/p312334/project/10--HGT_isolates/1--General_infomation/function/anvio9_pipeline/databases/COG/COG24/CATEGORIES.txt"

wrap_label <- function(x, width = 42) {
  vapply(strwrap(x, width = width, simplify = FALSE), paste, collapse = "\n", character(1))
}

cat_colors <- c(
  "C" = "#D73027",
  "B" = "#4575B4",
  "H" = "#1A9850",
  "R" = "#984EA3",
  "K" = "#FF7F00",
  "U" = "#A65628",
  "L" = "#4DAF4A"
)

freq_df <- read.csv(freq_file, stringsAsFactors = FALSE, check.names = FALSE)

cat_lines <- readLines(cat_file, warn = FALSE)
cat_lines <- cat_lines[grepl("^[A-Z]\t", cat_lines)]
cat_split <- strsplit(cat_lines, "\t", fixed = FALSE)
cat_df <- data.frame(
  COG24_CATEGORY = vapply(cat_split, `[[`, character(1), 1),
  category_name = vapply(cat_split, `[[`, character(1), 2),
  stringsAsFactors = FALSE
)

plot_df <- merge(freq_df, cat_df, by = "COG24_CATEGORY", all.x = TRUE, sort = FALSE)
plot_df$category_name[is.na(plot_df$category_name)] <- plot_df$COG24_CATEGORY[is.na(plot_df$category_name)]
plot_df$label <- paste0(plot_df$category_name, " (", plot_df$COG24_CATEGORY, ")")
plot_df$label_wrapped <- wrap_label(plot_df$label, width = 42)
plot_df$fill_group <- ifelse(plot_df$COG24_CATEGORY %in% names(cat_colors), plot_df$COG24_CATEGORY, "other")

ord <- order(-plot_df$frequency, plot_df$COG24_CATEGORY)
plot_df <- plot_df[ord, , drop = FALSE]

# Reverse factor levels so the flipped chart reads top-to-bottom in descending order.
plot_df$label_wrapped <- factor(plot_df$label_wrapped, levels = rev(plot_df$label_wrapped))

write.csv(plot_df, annotated_freq_file, row.names = FALSE, quote = TRUE)

fill_values <- c(cat_colors, "other" = "#CFCFCF")

p <- ggplot(plot_df, aes(x = label_wrapped, y = frequency, fill = fill_group)) +
  geom_col(width = 0.78, color = NA) +
  coord_flip() +
  scale_fill_manual(values = fill_values, guide = "none") +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)),
    labels = function(x) sprintf("%.0f%%", x * 100)
  ) +
  labs(
    x = NULL,
    y = "Frequency"
  ) +
  theme_minimal(base_size = 10.5, base_family = "Helvetica") +
  theme(
    axis.text.y = element_text(color = "black", size = 8.8, lineheight = 1.0),
    axis.text.x = element_text(color = "black", size = 9),
    axis.title.y = element_blank(),
    axis.title.x = element_text(face = "bold", size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "#D9D9D9", linewidth = 0.32),
    axis.line.x = element_line(color = "black", linewidth = 0.28),
    plot.margin = margin(5, 7, 5, 5)
  )

# Use the standard PDF device so text stays as editable/searchable text in the output PDF.
ggsave(plot_file, p, width = 4.3, height = 6.2, device = grDevices::pdf)
