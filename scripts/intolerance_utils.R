#!/usr/bin/env Rscript


make_bar <- function(df, title, fill_color) {
  ggplot( df ) +
    geom_col(aes(x = pos, y = distinct_alt_n), fill = fill_color, width = 1) +
    theme_minimal() +
    labs(title = title, x = NULL, y = NULL) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.title.y = element_blank(),
      plot.title = element_text(size = 10),
      plot.margin = margin(2, 5, 2, 5))
}


make_heatmap <- function(df, title, show_legend = FALSE) {
  ggplot(df, aes(x = pos, y = alt, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_manual(
      values = c("0" = "white", "1" = "#cce5ff", "2" = "#99c2e6", "3+" = "#6699cc"),
      name = "Allele Count"
    ) +
    theme_minimal() +
    labs(title = title, x = NULL, y = NULL) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.title.y = element_blank(),
      plot.title = element_text(size = 10),
      plot.margin = margin(2, 5, 2, 5),
      legend.position = if (show_legend) "right" else "none"
    )
}
