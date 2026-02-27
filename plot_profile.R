r
library(data.table)
library(dplyr)
library(ggplot2)

preprocess_data <- function(input_file) {
  df <- fread(input_file)
  df[strand == "-", distance := distance * -1]
  df <- df %>%
    group_by(distance) %>%
    summarise(mean_value = mean(mp_value, na.rm = TRUE), .groups = "drop") %>%
    arrange(distance)
  colnames(df) <- c("distance_bp", "mean_mp_percent")
  return(df)
}

create_tss_plot <- function(df) {
  x_limits <- c(-10000, 10000)
  y_limits <- c(0, 100)
  x_breaks <- c(-10000, 0, 10000)
  x_labels <- c("-10", "0", "10")
  y_breaks <- c(0, 25, 50, 75, 100)
  
  ggplot(df, aes(x = distance_bp, y = mean_mp_percent)) +
    geom_line(linewidth = 0.1, alpha = 1) +
    scale_x_continuous(
      limits = x_limits, 
      breaks = x_breaks, 
      labels = x_labels,
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = y_limits,
      breaks = y_breaks,
      expand = c(0, 0)
    ) +
    labs(
      x = "Distance from TSS (bp)",
      y = "Mean MP (%)"
    ) +
    theme(
      panel.background = element_blank(),
      plot.background = element_blank(),
      axis.line.y.left = element_line(linewidth = 0.25, lineend = "square"),
      axis.line.x.bottom = element_line(linewidth = 0.25, lineend = "square"),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_text(size = 5),
      axis.text = element_text(size = 5, color = "black"),
      legend.position = "none"
    )
}

# Usage
processed_data <- preprocess_data("MP_uniqueTSS_GRCh38.p14.tsv")
tss_plot <- create_tss_plot(processed_data)
print(tss_plot)
