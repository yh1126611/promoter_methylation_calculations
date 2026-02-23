library(ggplot2)
library(scales)

# Nature figure theme (linewidth=0.25 per Nature guidelines)
nature_theme <- theme(
  panel.background = element_blank(),
  plot.background = element_blank(),
  panel.grid = element_blank(),
  axis.line.y.left = element_line(linewidth = 0.25, color = "black"),
  axis.line.x.bottom = element_line(linewidth = 0.25, color = "black"),
  axis.ticks = element_blank(),
  axis.title = element_text(size = 6, color = "black"),
  axis.text = element_text(size = 5, color = "black"),
  legend.position = "none"
)

# Color schemes
gc_colors <- c("GC" = "#E7298A", "None" = "#BBBBBB")
mp_gradient_colors <- c(low = "#9EC9E2", high = "#0D4A70")

# Plot dimensions
axis_line_width <- 0.25
font_title_size <- 6
font_text_size <- 5
line_dotted_width <- 0.25

# Extended Data Fig. 1a
ggplot() +
  geom_bar(data = df, aes(x = genome_id, y = proportion, fill = feature_type), 
           stat = "identity", position = "fill") +
  scale_fill_manual(values = gc_colors, name = "Feature Type") + 
  labs(y = "Proportion", x = "Genome") +
  nature_theme

# Extended Data Fig. 1b
ggplot(df, aes(x = genome_id, y = proportion, fill = feature_color, group = feature_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_identity() +
  labs(y = "Proportion", x = "Genome") +
  nature_theme +
  theme(axis.line = element_line(linewidth = axis_line_width))

# Extended Data Fig. 1c
ggplot(df, aes(x = dinucleotide, y = total_length_mbp, fill = feature_color, group = dinucleotide)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_identity() +
  labs(y = "Length (Mbp)", x = "Dinucleotide") +
  scale_y_continuous(labels = label_number(scale = 1e-6)) +
  nature_theme +
  geom_hline(yintercept = 389659438, linetype = "dotted", 
             linewidth = line_dotted_width, color = "red")

# Extended Data Fig. 1d
ggplot(df, aes(x = modification_probability)) + 
  geom_histogram(aes(fill = after_stat(x)), binwidth = 0.1) +
  scale_fill_gradient(low = mp_gradient_colors["low"], 
                      high = mp_gradient_colors["high"], 
                      guide = "none") +
  labs(y = "Frequency", x = "Modification Probability (%)") +
  scale_y_continuous(labels = comma) +
  nature_theme

# Extended Data Fig. 1e
ggplot(stack_dataframe_sorted[stack_dataframe_sorted$chromosome != "MT", ], 
       aes(x = chromosome, y = cpg_count, fill = modification_probability)) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_gradient(low = mp_gradient_colors["low"], high = mp_gradient_colors["high"]) + 
  labs(y = "Number of CpGs", x = "Chromosome", fill = "Modification Probability") +
  scale_y_continuous(labels = comma) +
  nature_theme +
  theme(legend.position = "right",
        legend.key.size = unit(0.25, "cm"),
        legend.title = element_text(size = 8),
        legend.text = element_blank(),
        legend.ticks = element_blank())

# Extended Data Fig. 1f
chromosomal_proportion_plot <- ggplot(
  stack_dataframe_sorted[stack_dataframe_sorted$chromosome != "MT", ], 
  aes(x = chromosome, y = cpg_count, fill = modification_probability)
) + 
  geom_bar(stat = "identity", position = "fill") + 
  scale_fill_gradient(low = mp_gradient_colors["low"], high = mp_gradient_colors["high"]) +  
  labs(y = "Proportion of CpGs", x = "Chromosome", fill = "Modification Probability") +
  scale_y_continuous(labels = comma) +
  nature_theme +
  theme(legend.position = "none")
