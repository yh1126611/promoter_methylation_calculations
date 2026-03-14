library(data.table)
library(ggplot2)
library(ggdensity)
library(cowplot)
library(scales)

create_genomic_plots <- function(
  file_path,
  output_prefix = "genomic_plots",
  class_col = "Class",
  genome_size_col = "Genome size (bp)",
  promoter_size_col = "Promoter size (bp)",
  core_promoter_col = "Core promoter size (bp)",
  outlier_col = "Outlier",
  class_levels,
  class_colors,
  alpha_value = 0.14,
  plot_width_mm = 180
) {
  # Read and preprocess data
  df <- fread(file_path, dec = ".", fill = TRUE)
  
  # Forward-fill empty class values
  for (i in 1:nrow(df)) {
    if (df[[class_col]][i] == "") {
      df[[class_col]][i] <- df[[class_col]][i - 1]
    }
  }
  
  # Clean numeric columns: remove commas, convert to numeric
  df[[genome_size_col]] <- as.numeric(gsub(",", "", df[[genome_size_col]]))
  df[[promoter_size_col]] <- as.numeric(gsub(",", "", df[[promoter_size_col]]))
  df[[core_promoter_col]] <- as.numeric(df[[core_promoter_col]])
  
  # Filter data
  df <- df[df[[promoter_size_col]] < 10000, ]
  df <- df[!is.na(df[[core_promoter_col]]), ]
  df <- df[df[[outlier_col]] != "Yes", ]
  
  # Ensure class is factor with specified levels
  df[[class_col]] <- factor(df[[class_col]], levels = class_levels)
  
  # Common theme
  common_theme <- theme(
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.line.y.left = element_line(linewidth = 0.1, color = "black"),
    axis.line.x.bottom = element_line(linewidth = 0.1, color = "black"),
    axis.title.y = element_text(size = 5, color = "black"),
    axis.text.y = element_text(size = 5, color = "black"),
    axis.title.x = element_text(size = 5, color = "black"),
    axis.text.x = element_text(size = 5, color = "black"),
    legend.position = "none"
  )
  
  # Function to create density + points plot
  create_density_plot <- function(df, x_col, y_col, x_label, y_label, x_scale = 1) {
    p <- ggplot(df) +
      coord_cartesian(clip = "off")
    
    # Add density contours for each class
    for (i in seq_along(class_levels)) {
      class_val <- class_levels[i]
      subset_df <- df[df[[class_col]] == class_val, ]
      color_val <- class_colors[class_val]
      
      p <- p +
        stat_density_2d(
          data = subset_df,
          mapping = aes(x = .data[[x_col]], y = .data[[y_col]], fill = after_stat(level)),
          geom = "polygon",
          alpha = alpha_value
        ) +
        scale_fill_gradient(low = "white", high = color_val) +
        new_scale_fill()
    }
    
    p <- p +
      geom_point(
        mapping = aes(x = .data[[x_col]], y = .data[[y_col]], color = .data[[class_col]]),
        size = 0.3
      ) +
      scale_color_manual(values = class_colors, na.value = "pink") +
      scale_y_continuous(name = y_label, expand = c(0, 0), limits = c(0, max(df[[y_col]], na.rm = TRUE))) +
      scale_x_continuous(
        name = x_label,
        labels = label_number(scale = x_scale),
        expand = c(0, 0),
        limits = c(0, max(df[[x_col]], na.rm = TRUE))
      ) +
      common_theme
    
    return(p)
  }
  
  # Create three plots
  pg <- create_density_plot(
    df, genome_size_col, promoter_size_col,
    "Genome size (Gb)", "Promoter size (kb)", 1e-9
  )
  
  cg <- create_density_plot(
    df, genome_size_col, core_promoter_col,
    "Genome size (Gb)", "Core promoter size (bp)", 1
  )
  
  cp <- create_density_plot(
    df, promoter_size_col, core_promoter_col,
    "Promoter size (kb)", "Core promoter size (bp)", 1e-3
  )
  
  # Align plots
  genome_size_aligned <- align_plots(
    pg + theme(plot.margin = margin(5, 0, 0, 0)),
    cg + theme(plot.margin = margin(5, 0, 0, 0)),
    align = "v", axis = "lr"
  )
  
  core_aligned <- align_plots(
    cp + theme(plot.margin = margin(5, 0, 0, 0)),
    cg + theme(plot.margin = margin(5, 0, 0, 0)),
    align = "h", axis = "tb"
  )
  
  # Create legend data frame
  legend_df <- data.frame(
    Class = names(class_colors),
    X = seq_along(class_colors)
  )
  legend_df$Class <- factor(legend_df$Class, levels = names(class_colors))
  
  # Create legend
  legend <- as_ggplot(get_legend(
    ggplot(legend_df) +
      geom_point(mapping = aes(color = Class, x = X, y = X)) +
      scale_color_manual(name = "Class:", values = class_colors) +
      guides(color = guide_legend(ncol = 1, override.aes = list(size = 2))) +
      theme(
        plot.background = element_blank(),
        panel.background = element_blank(),
        legend.background = element_blank(),
        legend.text = element_text(color = "black", size = 5.5),
        legend.title = element_text(color = "black", size = 6, face = "bold"),
        legend.key.size = unit(0.75, "lines")
      )
  ))
  
  # Combine into lengths plot
  lengths_plot <- plot_grid(
    genome_size_aligned[[1]], core_aligned[[2]], core_aligned[[1]], legend,
    ncol = 4, labels = c("h", "i", "j"), label_size = 8
  )
  
  # Generate output filenames from prefix
  png_filename <- paste0(output_prefix, "_lengths.png")
  pdf_filename <- paste0(output_prefix, "_lengths.pdf")
  
  # Save plots
  ggsave(png_filename, lengths_plot, width = plot_width_mm, units = "mm", height = 12, dpi = 600)
  ggsave(pdf_filename, lengths_plot, width = plot_width_mm, units = "mm", height = 12)
  
  cat("Plots saved:\n", png_filename, "\n", pdf_filename, "\n")
  
  invisible(list(png = png_filename, pdf = pdf_filename, plot = lengths_plot))
}

# Example usage
class_levels <- c("Mammalia", "Aves", "Reptilia", "Amphibia", "Actinopterygii", "Chondrichthyes")
class_colors <- c(Mammalia = "#E69F00", Aves = "#00796B", Reptilia = "#A6D609", 
                  Amphibia = "#984EA3", Actinopterygii = "#56B4E9", Chondrichthyes = "#0072B2")

results <- create_genomic_plots(
  file_path = "your_file.tsv",
  output_prefix = "my_analysis",
  class_levels = class_levels,
  class_colors = class_colors
)

print(results$plot)
