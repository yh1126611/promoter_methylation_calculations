plot_gene_mp <- function(MATRIX, COLUMN, TYPE, TEXT, ALPHA) {
  class_colors=c("Amphibia"="#984EA3", "Aves"="#00796B", "Reptilia"="#A6D609", "Actinopterygii"="#56B4E9", "Sarcopterygii"="#A6761D", "Chondrichthyes"="#0072B2", "Mammalia"="#E69F00")
  ggplot(MATRIX[MATRIX[[COLUMN]] %in% TYPE, ],aes(x = Distance, y = MP)) +
    #annotation_raster(PIC, xmin = -5000, xmax = Inf, ymin = 25, ymax = Inf) +
    geom_point(aes(color=`Class (강)`), alpha = ALPHA, stroke = NA, shape = 19, size=0.25) +
    scale_color_manual(values=class_colors, na.value="black")+
    theme(
      plot.background = element_blank(),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.line.y.left = element_line(linewidth = 0.1),
      axis.line.x.bottom = element_line(linewidth = 0.1),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 5.5, color="black"),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 5.5, color="black"),
      legend.position = "none"
    ) +
    scale_y_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
    scale_x_continuous(labels = c("-10", "0", "10"), limits = c(-10000, 10000), breaks = c(-10000, 0, 10000)) +
    # annotate("text", y=25, x=10000, vjust=1, hjust="inward", label = TEXT, color = "black", size = 1.94) +
    coord_cartesian(clip = "off")
}

# Usage
plot_gene_mp(data_matrix_ACTB_merged, "Class", c("Mammalia"), "nGenome = 31\nnCpG = 19,432", 1)->mammal_gene
